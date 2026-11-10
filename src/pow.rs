/* origin: FreeBSD /usr/src/lib/msun/src/e_pow.c */
/*
 * ====================================================
 * Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

// |x|^y = 2^(y*log2(|x|)). High/low parts retain precision during reduction.
// Handle exceptional values and integer exponents before evaluating the approximation.

use core::f64::consts::{LN_2, LOG2_E};

use crate::{
    exp::exp_correction, fabs, get_high_word, get_low_word, log_constants, scalbn::scalbn, sqrt,
    with_set_high_word, with_set_low_word,
};

const LOG_BASES: [f64; 2] = [1.0, 3.0 / 2.0];
const LOG2_BASE_HI: [f64; 2] = [0.0, log_constants::LOG2_THREE_HALVES.hi];
const LOG2_BASE_LO: [f64; 2] = [0.0, log_constants::LOG2_THREE_HALVES.lo];
const SUBNORMAL_SCALE: f64 = crate::expi(2.0, f64::MANTISSA_DIGITS as isize);
const OVERFLOW_FACTOR: f64 = f64::MAX;
const UNDERFLOW_FACTOR: f64 = f64::MIN_POSITIVE;

// Fitted coefficients for the logarithm correction.
const LOG_L1: f64 = 5.99999999999994648725e-01;
const LOG_L2: f64 = 4.28571428578550184252e-01;
const LOG_L3: f64 = 3.33333329818377432918e-01;
const LOG_L4: f64 = 2.72728123808534006489e-01;
const LOG_L5: f64 = 2.30660745775561754067e-01;
const LOG_L6: f64 = 2.06975017800338417784e-01;

const LN_2_HI: f64 = log_constants::LN_TWO.hi;
const LN_2_LO: f64 = log_constants::LN_TWO.lo;
const LOG_SCALE: f64 = 2.0 * LOG2_E / 3.0;
const LOG_SCALE_HI: f64 = log_constants::LOG_SCALE.hi;
const LOG_SCALE_LO: f64 = log_constants::LOG_SCALE.lo;
const LOG2_E_HI: f64 = log_constants::LOG2_E.hi;
const LOG2_E_LO: f64 = log_constants::LOG2_E.lo;
// Correction at the overflow rounding boundary.
const OVERFLOW_CORRECTION: f64 = LOG2_E * crate::expi(2.0, -(f64::MANTISSA_DIGITS as isize + 1));

// IEEE-754 fields in the high 32-bit word.
const FRACTION_BITS: i32 = f64::MANTISSA_DIGITS as i32 - 1;
const HIGH_FRACTION_BITS: i32 = FRACTION_BITS - u32::BITS as i32;
const EXPONENT_BIAS: i32 = f64::MAX_EXP - 1;
const EXPONENT_UNIT: i32 = 1 << HIGH_FRACTION_BITS;
const HIGH_FRACTION_MASK: i32 = EXPONENT_UNIT - 1;
const ABS_HIGH_MASK: i32 = i32::MAX;
const HALF_EXPONENT_BIAS_WORD: u32 = (((EXPONENT_BIAS + 1) / 2) as u32) << HIGH_FRACTION_BITS;
const HALF_EXPONENT_UNIT_WORD: u32 = (EXPONENT_UNIT as u32) >> 1;
const ONE_HIGH_WORD: i32 = get_high_word(1.0) as i32;
const TWO_HIGH_WORD: i32 = get_high_word(2.0) as i32;
const HALF_HIGH_WORD: i32 = get_high_word(0.5) as i32;
const INFINITY_HIGH_WORD: i32 = get_high_word(f64::INFINITY) as i32;
const MIN_NORMAL_HIGH_WORD: i32 = get_high_word(f64::MIN_POSITIVE) as i32;
const EVEN_INTEGER_HIGH_WORD: i32 = get_high_word(SUBNORMAL_SCALE) as i32;
const LARGE_POWER_HIGH_WORD: i32 = get_high_word(crate::expi(2.0, (i32::BITS - 1) as isize)) as i32;
const HUGE_POWER_HIGH_WORD: i32 = get_high_word(crate::expi(2.0, u64::BITS as isize)) as i32;
const SQRT_THREE_HALVES_FRACTION: i32 = get_high_word(sqrt(3.0 / 2.0)) as i32 & HIGH_FRACTION_MASK;
const SQRT_THREE_FRACTION: i32 = get_high_word(sqrt(3.0)) as i32 & HIGH_FRACTION_MASK;
const EXP2_OVERFLOW_HIGH_WORD: i32 = get_high_word(f64::MAX_EXP as f64) as i32;
const ZERO_ROUNDING_EXPONENT: i32 = f64::MIN_EXP - f64::MANTISSA_DIGITS as i32 - 1;
const EXP2_ZERO_HIGH_WORD: i32 = get_high_word(ZERO_ROUNDING_EXPONENT as f64) as i32;

/// Power function
///
/// ```
/// # use trig_const::pow;
/// let _a = pow(10.0, 2.5);
/// assert_eq!(pow(10.0, 2.0), 100.0);
/// ```
#[must_use]
pub const fn pow(x: f64, y: f64) -> f64 {
    #[cfg(feature = "nightly")]
    {
        #[cfg(feature = "std")]
        {
            std::intrinsics::const_eval_select((x, y), pow_inner, f64::powf)
        }
        #[cfg(not(feature = "std"))]
        {
            core::intrinsics::const_eval_select((x, y), pow_inner, libm::pow)
        }
    }
    #[cfg(not(feature = "nightly"))]
    {
        pow_inner(x, y)
    }
}

const fn pow_inner(x: f64, y: f64) -> f64 {
    let t1: f64;

    let (hx, lx): (i32, u32) = (get_high_word(x) as i32, get_low_word(x));
    let (hy, ly): (i32, u32) = (get_high_word(y) as i32, get_low_word(y));

    let mut ix: i32 = hx & ABS_HIGH_MASK;
    let iy: i32 = hy & ABS_HIGH_MASK;

    /* x**0 = 1, even if x is NaN */
    if ((iy as u32) | ly) == 0 {
        return 1.0;
    }

    /* 1**y = 1, even if y is NaN */
    if hx == ONE_HIGH_WORD && lx == 0 {
        return 1.0;
    }

    /* NaN if either arg is NaN */
    if ix > INFINITY_HIGH_WORD
        || (ix == INFINITY_HIGH_WORD && lx != 0)
        || iy > INFINITY_HIGH_WORD
        || (iy == INFINITY_HIGH_WORD && ly != 0)
    {
        return x + y;
    }

    /* determine if y is an odd int when x < 0
     * yisint = 0       ... y is not an integer
     * yisint = 1       ... y is an odd int
     * yisint = 2       ... y is an even int
     */
    let mut yisint: i32 = 0;
    let mut k: i32;
    let mut j: i32;
    if hx < 0 {
        if iy >= EVEN_INTEGER_HIGH_WORD {
            yisint = 2; /* even integer y */
        } else if iy >= ONE_HIGH_WORD {
            k = (iy >> HIGH_FRACTION_BITS) - EXPONENT_BIAS; /* exponent */

            if k > HIGH_FRACTION_BITS {
                j = (ly >> (FRACTION_BITS - k)) as i32;

                if (j << (FRACTION_BITS - k)) == (ly as i32) {
                    yisint = 2 - (j & 1);
                }
            } else if ly == 0 {
                j = iy >> (HIGH_FRACTION_BITS - k);

                if (j << (HIGH_FRACTION_BITS - k)) == iy {
                    yisint = 2 - (j & 1);
                }
            }
        }
    }

    if ly == 0 {
        /* special value of y */
        if iy == INFINITY_HIGH_WORD {
            /* y is +-inf */

            return if ((ix - ONE_HIGH_WORD) | (lx as i32)) == 0 {
                /* (-1)**+-inf is 1 */
                1.0
            } else if ix >= ONE_HIGH_WORD {
                /* (|x|>1)**+-inf = inf,0 */
                if hy >= 0 {
                    y
                } else {
                    0.0
                }
            } else {
                /* (|x|<1)**+-inf = 0,inf */
                if hy >= 0 {
                    0.0
                } else {
                    -y
                }
            };
        }

        if iy == ONE_HIGH_WORD {
            /* y is +-1 */
            return if hy >= 0 { x } else { 1.0 / x };
        }

        if hy == TWO_HIGH_WORD {
            /* y is 2 */
            return x * x;
        }

        if hy == HALF_HIGH_WORD {
            /* y is 0.5 */
            if hx >= 0 {
                /* x >= +0 */
                return sqrt(x);
            }
        }
    }

    let mut ax: f64 = fabs(x);
    if lx == 0 {
        /* special value of x */
        if ix == INFINITY_HIGH_WORD || ix == 0 || ix == ONE_HIGH_WORD {
            /* x is +-0,+-inf,+-1 */
            let mut z: f64 = ax;

            if hy < 0 {
                /* z = (1/|x|) */
                z = 1.0 / z;
            }

            if hx < 0 {
                if ((ix - ONE_HIGH_WORD) | yisint) == 0 {
                    z = f64::NAN; /* (-1)**non-int is NaN */
                } else if yisint == 1 {
                    z = -z; /* (x<0)**odd = -(|x|**odd) */
                }
            }

            return z;
        }
    }

    let mut s: f64 = 1.0; /* sign of result */
    if hx < 0 {
        if yisint == 0 {
            /* (x<0)**(non-int) is NaN */
            return f64::NAN;
        }

        if yisint == 1 {
            /* (x<0)**(odd int) */
            s = -1.0;
        }
    }

    // Large powers use a separate expansion near x=1.
    let t2: f64 = if iy > LARGE_POWER_HIGH_WORD {
        /* if |y| > 2**31 */
        if iy > HUGE_POWER_HIGH_WORD {
            /* if |y| > 2**64, must o/uflow */
            if ix <= (ONE_HIGH_WORD - 1) {
                return if hy < 0 {
                    OVERFLOW_FACTOR * OVERFLOW_FACTOR
                } else {
                    UNDERFLOW_FACTOR * UNDERFLOW_FACTOR
                };
            }

            if ix >= ONE_HIGH_WORD {
                return if hy > 0 {
                    OVERFLOW_FACTOR * OVERFLOW_FACTOR
                } else {
                    UNDERFLOW_FACTOR * UNDERFLOW_FACTOR
                };
            }
        }

        /* over/underflow if x is not close to one */
        if ix < (ONE_HIGH_WORD - 1) {
            return if hy < 0 {
                s * OVERFLOW_FACTOR * OVERFLOW_FACTOR
            } else {
                s * UNDERFLOW_FACTOR * UNDERFLOW_FACTOR
            };
        }
        if ix > ONE_HIGH_WORD {
            return if hy > 0 {
                s * OVERFLOW_FACTOR * OVERFLOW_FACTOR
            } else {
                s * UNDERFLOW_FACTOR * UNDERFLOW_FACTOR
            };
        }

        // Near x=1, use log(1+t) = t-t^2/2+t^3/3-t^4/4.
        let t: f64 = ax - 1.0; /* t has 20 trailing zeros */
        let w: f64 = (t * t) * (0.5 - t * (1.0 / 3.0 - t * 0.25));
        let u: f64 = LOG2_E_HI * t;
        let v: f64 = t * LOG2_E_LO - w * LOG2_E;
        t1 = with_set_low_word(u + v, 0);
        v - (t1 - u)
    } else {
        let mut n: i32 = 0;

        if ix < MIN_NORMAL_HIGH_WORD {
            /* take care subnormal number */
            ax *= SUBNORMAL_SCALE;
            n -= f64::MANTISSA_DIGITS as i32;
            ix = get_high_word(ax) as i32;
        }

        n += (ix >> HIGH_FRACTION_BITS) - EXPONENT_BIAS;
        j = ix & HIGH_FRACTION_MASK;

        /* determine interval */
        let k: i32;
        ix = j | ONE_HIGH_WORD; /* normalize ix */
        if j <= SQRT_THREE_HALVES_FRACTION {
            /* |x|<sqrt(3/2) */
            k = 0;
        } else if j < SQRT_THREE_FRACTION {
            /* |x|<sqrt(3)   */
            k = 1;
        } else {
            k = 0;
            n += 1;
            ix -= EXPONENT_UNIT;
        }
        ax = with_set_high_word(ax, ix as u32);

        /* compute ss = s_h+s_l = (x-1)/(x+1) or (x-1.5)/(x+1.5) */
        let u: f64 = ax - LOG_BASES[k as usize];
        let v: f64 = 1.0 / (ax + LOG_BASES[k as usize]);
        let ss: f64 = u * v;
        let s_h = with_set_low_word(ss, 0);

        // Build ax+base's high word without rounding its low word.
        let t_h: f64 = with_set_high_word(
            0.0,
            ((ix as u32 >> 1) | HALF_EXPONENT_BIAS_WORD)
                + HALF_EXPONENT_UNIT_WORD
                + ((k as u32) << (HIGH_FRACTION_BITS - 2)),
        );
        let t_l: f64 = ax - (t_h - LOG_BASES[k as usize]);
        let s_l: f64 = v * ((u - s_h * t_h) - s_h * t_l);

        /* compute log(ax) */
        let s2: f64 = ss * ss;
        let mut r: f64 = s2
            * s2
            * (LOG_L1
                + s2 * (LOG_L2 + s2 * (LOG_L3 + s2 * (LOG_L4 + s2 * (LOG_L5 + s2 * LOG_L6)))));
        r += s_l * (s_h + ss);
        let s2: f64 = s_h * s_h;
        let t_h: f64 = with_set_low_word(3.0 + s2 + r, 0);
        let t_l: f64 = r - ((t_h - 3.0) - s2);

        /* u+v = ss*(1+...) */
        let u: f64 = s_h * t_h;
        let v: f64 = s_l * t_h + t_l * ss;

        /* 2/(3log2)*(ss+...) */
        let p_h: f64 = with_set_low_word(u + v, 0);
        let p_l = v - (p_h - u);
        let z_h: f64 = LOG_SCALE_HI * p_h;
        let z_l: f64 = LOG_SCALE_LO * p_h + p_l * LOG_SCALE + LOG2_BASE_LO[k as usize];

        // Combine the exponent, base logarithm, and correction.
        let t: f64 = n as f64;
        t1 = with_set_low_word(((z_h + z_l) + LOG2_BASE_HI[k as usize]) + t, 0);
        z_l - (((t1 - t) - LOG2_BASE_HI[k as usize]) - z_h)
    };

    /* split up y into y1+y2 and compute (y1+y2)*(t1+t2) */
    let y1: f64 = with_set_low_word(y, 0);
    let p_l: f64 = (y - y1) * t1 + y * t2;
    let mut p_h: f64 = y1 * t1;
    let z: f64 = p_l + p_h;
    let mut j: i32 = get_high_word(z) as i32;
    let i: i32 = get_low_word(z) as i32;

    if j >= EXP2_OVERFLOW_HIGH_WORD {
        /* z >= 1024 */
        if (j - EXP2_OVERFLOW_HIGH_WORD) | i != 0 {
            /* if z > 1024 */
            return s * OVERFLOW_FACTOR * OVERFLOW_FACTOR; /* overflow */
        }

        if p_l + OVERFLOW_CORRECTION > z - p_h {
            return s * OVERFLOW_FACTOR * OVERFLOW_FACTOR; /* overflow */
        }
    } else if (j & ABS_HIGH_MASK) >= (EXP2_ZERO_HIGH_WORD & ABS_HIGH_MASK) {
        /* z <= -1075 */
        // FIXME: instead of abs(j) use unsigned j

        if (((j as u32) - EXP2_ZERO_HIGH_WORD as u32) | (i as u32)) != 0 {
            /* z < -1075 */
            return s * UNDERFLOW_FACTOR * UNDERFLOW_FACTOR; /* underflow */
        }

        if p_l <= z - p_h {
            return s * UNDERFLOW_FACTOR * UNDERFLOW_FACTOR; /* underflow */
        }
    }

    /* compute 2**(p_h+p_l) */
    let i: i32 = j & ABS_HIGH_MASK;
    k = (i >> HIGH_FRACTION_BITS) - EXPONENT_BIAS;
    let mut n: i32 = 0;

    if i > HALF_HIGH_WORD {
        /* if |z| > 0.5, set n = [z+0.5] */
        n = j + (EXPONENT_UNIT >> (k + 1));
        k = ((n & ABS_HIGH_MASK) >> HIGH_FRACTION_BITS) - EXPONENT_BIAS; /* new k for n */
        let t: f64 = with_set_high_word(0.0, (n & !(HIGH_FRACTION_MASK >> k)) as u32);
        n = ((n & HIGH_FRACTION_MASK) | EXPONENT_UNIT) >> (HIGH_FRACTION_BITS - k);
        if j < 0 {
            n = -n;
        }
        p_h -= t;
    }

    let t: f64 = with_set_low_word(p_l + p_h, 0);
    let u: f64 = t * LN_2_HI;
    let v: f64 = (p_l - (t - p_h)) * LN_2 + t * LN_2_LO;
    let mut z: f64 = u + v;
    let w: f64 = v - (z - u);
    let t1 = exp_correction(z);
    let r: f64 = (z * t1) / (t1 - 2.0) - (w + z * w);
    z = 1.0 - (r - z);
    j = get_high_word(z) as i32;
    j += n << HIGH_FRACTION_BITS;

    if (j >> HIGH_FRACTION_BITS) <= 0 {
        /* subnormal output */
        z = scalbn(z, n);
    } else {
        z = with_set_high_word(z, j as u32);
    }

    s * z
}

#[cfg(test)]
mod tests {
    use super::pow;

    #[test]
    fn exponent_limits_const() {
        const MIN_SUBNORMAL: f64 = f64::MIN_POSITIVE * f64::EPSILON;
        const MAX_POWER_OF_TWO: f64 = crate::expi(2.0, (f64::MAX_EXP - 1) as isize);
        const LAST_ODD_INTEGER: f64 = crate::expi(2.0, f64::MANTISSA_DIGITS as isize) - 1.0;
        assert_eq!((const { pow(2.0, 1023.0) }), MAX_POWER_OF_TWO);
        assert!((const { pow(2.0, 1024.0) }).is_infinite());
        assert_eq!((const { pow(2.0, -1074.0) }), MIN_SUBNORMAL);
        assert_eq!((const { pow(2.0, -1075.0) }).to_bits(), 0);
        assert_eq!(
            (const { pow(-2.0, -1075.0) }).to_bits(),
            (-0.0_f64).to_bits()
        );
        assert_eq!((const { pow(-1.0, LAST_ODD_INTEGER) }), -1.0);
        assert_eq!((const { pow(-1.0, LAST_ODD_INTEGER + 1.0) }), 1.0);
    }
}
