#![doc = include_str!("../README.md")]
#![cfg_attr(not(feature = "std"), no_std)]
#![forbid(unsafe_code)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]
#![allow(internal_features)]
#![cfg_attr(feature = "nightly", feature(core_intrinsics, const_eval_select))]

#[macro_use]
pub(crate) mod macros;

mod acos;
mod acosh;
mod asin;
mod asinh;
mod atan;
mod atan2;
mod atanh;
mod cos;
mod exp;
mod expm1;
mod floor;
mod k_cos;
mod k_sin;
pub(crate) mod k_tan;
mod ln;
pub(crate) mod log1p;
mod pow;
mod rem_pio2;
mod rem_pio2_large;
pub(crate) mod scalbn;
mod sin;
mod tan;
pub use acos::acos;
pub use acosh::acosh;
pub use asin::asin;
pub use asinh::asinh;
pub use atan::atan;
pub use atan2::atan2;
pub use atanh::atanh;
pub use cos::cos;
pub use exp::exp;
pub use floor::floor;
pub use ln::ln;
pub use pow::pow;
pub use sin::sin;
pub use tan::tan;

/// Cotangent
///
/// ```
/// # use trig_const::cot;
/// # use core::f64::consts::PI;
/// # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
/// const COT_PI_4: f64 = cot(PI / 4.0);
/// float_eq(COT_PI_4, 1.0);
/// ```
pub const fn cot(x: f64) -> f64 {
    let sin_calc = sin(x);
    if sin_calc == 0.0 {
        f64::INFINITY
    } else {
        cos(x) / sin_calc
    }
}

/// Secant
///
/// ```
/// # use trig_const::sec;
/// # use core::f64::consts::PI;
/// const SEC_PI: f64 = sec(PI);
/// assert_eq!(SEC_PI, -1.0);
/// ```
pub const fn sec(x: f64) -> f64 {
    let cos_calc = cos(x);
    if cos_calc == 0.0 {
        f64::INFINITY
    } else {
        1.0 / cos_calc
    }
}

/// Cosecant
///
/// ```
/// # use trig_const::csc;
/// # use core::f64::consts::PI;
/// const CSC_PI_2: f64 = csc(PI / 2.0);
/// assert_eq!(CSC_PI_2, 1.0);
/// ```
pub const fn csc(x: f64) -> f64 {
    let sin_calc = sin(x);
    if sin_calc == 0.0 {
        f64::INFINITY
    } else {
        1.0 / sin_calc
    }
}

/// Hyperbolic Sine
///
/// ```
/// # use trig_const::sinh;
/// const SINH_0: f64 = sinh(0.0);
/// assert_eq!(SINH_0, 0.0);
/// ```
pub const fn sinh(x: f64) -> f64 {
    nightly_exp!(sinh, sinh_inner, x)
}

const fn sinh_inner(x: f64) -> f64 {
    let magnitude = x.abs();
    let result = if magnitude < 709.782712893384 {
        // expm1 retains the low bits that exp(x) - 1 loses near zero.
        if magnitude < 1.4901161193847656e-8 {
            return x;
        }
        let t = expm1::expm1(magnitude);
        if magnitude < 1.0 {
            0.5 * (2.0 * t - t * t / (t + 1.0))
        } else {
            0.5 * (t + t / (t + 1.0))
        }
    } else {
        // Scale before multiplication so exp(x) need not fit in an f64.
        let t = exp(0.5 * magnitude);
        (0.5 * t) * t
    };
    if x.is_sign_negative() {
        -result
    } else {
        result
    }
}

/// Hyperbolic Cosine
///
/// ```
/// # use trig_const::cosh;
/// const COSH_0: f64 = cosh(0.0);
/// assert_eq!(COSH_0, 1.0);
/// ```
pub const fn cosh(x: f64) -> f64 {
    nightly_exp!(cosh, cosh_inner, x)
}

const fn cosh_inner(x: f64) -> f64 {
    (exp(x) + exp(-x)) / 2.0
}

/// x^pow
pub const fn expi(x: f64, mut pow: isize) -> f64 {
    let mut o = 1.0;

    while pow > 0 {
        o *= x;
        pow -= 1;
    }
    while pow < 0 {
        o /= x;
        pow += 1;
    }

    o
}

/// Factorial (x!)
pub const fn factorial(mut x: f64) -> f64 {
    if x == 0.0 || x == 1.0 {
        1.0
    } else {
        let mut s = 1.0;
        while x > 1.0 {
            s *= x;
            x -= 1.0;
        }
        s
    }
}

/// Square root, rounded to the nearest representable value.
///
/// Preserves signed zero and positive infinity. Negative inputs return NaN.
pub const fn sqrt(x: f64) -> f64 {
    nightly_exp!(sqrt, sqrt_inner, x)
}

const fn sqrt_inner(x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    } else if x.is_infinite() || x == 0.0 {
        return x;
    }

    // Normalize x = significand * 2^(exponent - 52), including subnormals.
    let bits = x.to_bits();
    let mut significand = bits & 0x000f_ffff_ffff_ffff;
    let mut exponent = ((bits >> 52) as i32) - 1023;
    if exponent == -1023 {
        let shift = significand.leading_zeros() - 11;
        significand <<= shift;
        exponent = -1022 - shift as i32;
    } else {
        significand |= 1_u64 << 52;
    }

    // Make the exponent even so its square root is an exact power of two.
    if exponent & 1 != 0 {
        significand <<= 1;
        exponent -= 1;
    }
    let radicand = (significand as u128) << 52;

    // Integer square root by long division, consuming two bits per step.
    // After 53 steps, root = floor(sqrt(radicand)) and
    // remainder = radicand - root^2. All arithmetic here is exact.
    let mut root = 0_u64;
    let mut remainder = 0_u64;
    let mut pairs = 53;
    while pairs > 0 {
        pairs -= 1;
        remainder = (remainder << 2) | ((radicand >> (pairs * 2)) & 3) as u64;
        root <<= 1;
        let trial = (root << 1) | 1;
        if remainder >= trial {
            remainder -= trial;
            root += 1;
        }
    }

    // Round up precisely when sqrt(radicand) > root + 1/2.
    // A tie would require the integer radicand to equal root^2 + root + 1/4.
    if remainder > root {
        root += 1;
    }

    // Every nonzero f64 square root is normal. Addition also handles a rounded
    // significand carrying into the next exponent.
    let result_exponent = (exponent / 2 + 1023) as u64;
    f64::from_bits((result_exponent << 52) + (root - (1_u64 << 52)))
}

pub const fn fabs(x: f64) -> f64 {
    x.abs()
}

const fn with_set_high_word(f: f64, hi: u32) -> f64 {
    let mut tmp = f.to_bits();
    tmp &= 0x00000000_ffffffff;
    tmp |= (hi as u64) << 32;
    f64::from_bits(tmp)
}
const fn with_set_low_word(f: f64, lo: u32) -> f64 {
    let mut tmp = f.to_bits();
    tmp &= 0xffffffff_00000000;
    tmp |= lo as u64;
    f64::from_bits(tmp)
}
const fn get_high_word(x: f64) -> u32 {
    (x.to_bits() >> 32) as u32
}

const fn get_low_word(x: f64) -> u32 {
    x.to_bits() as u32
}

#[cfg(test)]
mod tests {
    use core::f64::consts::{E, PI};

    use crate::{cos, cosh, exp, expi, factorial, ln, sin, sinh, sqrt};

    macro_rules! float_eq {
        ($lhs:expr, $rhs:expr) => {
            assert!(($lhs - $rhs).abs() < 0.0001, "lhs: {}, rhs: {}", $lhs, $rhs);
        };
    }

    #[test]
    fn test_factorial() {
        assert_eq!(factorial(0.0), 1.0);
        assert_eq!(factorial(1.0), 1.0);
        assert_eq!(factorial(2.0), 2.0);
        assert_eq!(factorial(3.0), 6.0);
        assert_eq!(factorial(4.0), 24.0);
        assert_eq!(factorial(5.0), 120.0);
    }

    #[test]
    fn test_expi() {
        assert_eq!(expi(2.0, 0), 1.0);
        assert_eq!(expi(2.0, 4), 16.0);
        assert_eq!(expi(2.0, 5), 32.0);
        assert_eq!(expi(3.0, 3), 27.0);
    }

    #[test]
    fn test_exp() {
        float_eq!(exp(0.0), 1.0);
        float_eq!(exp(1.0), E);
    }

    #[test]
    fn test_sqrt() {
        float_eq!(sqrt(0.0), 0.0);
        float_eq!(sqrt(1.0), 1.0);
        float_eq!(sqrt(4.0), 2.0);
        float_eq!(sqrt(9.0), 3.0);
        float_eq!(sqrt(16.0), 4.0);
        float_eq!(sqrt(25.0), 5.0);
    }

    #[test]
    fn test_sqrt_extreme_magnitudes_const() {
        const MIN_SUBNORMAL: f64 = f64::from_bits(1);
        const MAX_SUBNORMAL: f64 = f64::from_bits(0x000f_ffff_ffff_ffff);

        assert_eq!((const { sqrt(1e-20) }), 1e-10);
        assert_eq!((const { sqrt(1e20) }), 1e10);
        assert_eq!((const { sqrt(MIN_SUBNORMAL) }), MIN_SUBNORMAL.sqrt());
        assert_eq!((const { sqrt(MAX_SUBNORMAL) }), MAX_SUBNORMAL.sqrt());
        assert_eq!(
            (const { sqrt(f64::MIN_POSITIVE) }),
            f64::MIN_POSITIVE.sqrt()
        );
        assert_eq!((const { sqrt(f64::MAX) }), f64::MAX.sqrt());
    }

    #[test]
    fn test_sqrt_exponents_and_rounding() {
        // Exercise every normal exponent and significands near rounding boundaries.
        for exponent in 1..0x7ff_u64 {
            for fraction in [0, 1, 0x0007_ffff_ffff_ffff, 0x000f_ffff_ffff_ffff] {
                let x = f64::from_bits((exponent << 52) | fraction);
                // Check the const implementation even when nightly forwards runtime calls.
                assert_eq!(super::sqrt_inner(x).to_bits(), x.sqrt().to_bits());
            }
        }
    }

    #[test]
    fn test_sqrt_special_values() {
        assert_eq!(super::sqrt_inner(0.0).to_bits(), 0.0_f64.to_bits());
        assert_eq!(super::sqrt_inner(-0.0).to_bits(), (-0.0_f64).to_bits());
        assert_eq!(super::sqrt_inner(f64::INFINITY), f64::INFINITY);
        assert!(super::sqrt_inner(-1.0).is_nan());
        assert!(super::sqrt_inner(f64::NEG_INFINITY).is_nan());
        assert!(super::sqrt_inner(f64::NAN).is_nan());
    }

    #[test]
    fn test_sqrt_sampled_significands() {
        // Fixed seed keeps this coverage deterministic, including subnormal inputs.
        let mut bits = 0x1234_5678_9abc_def0_u64;
        for _ in 0..10_000 {
            bits = bits.wrapping_mul(6364136223846793005).wrapping_add(1);
            for sample in [bits & 0x7fef_ffff_ffff_ffff, bits & 0x000f_ffff_ffff_ffff] {
                let x = f64::from_bits(sample);
                assert_eq!(super::sqrt_inner(x).to_bits(), x.sqrt().to_bits());
            }
        }
    }

    #[test]
    fn test_cos() {
        float_eq!(cos(0.0), 0.0_f64.cos());
        float_eq!(cos(1.0), 1.0_f64.cos());
        float_eq!(cos(PI), PI.cos());
        float_eq!(cos(PI * 8.0), (PI * 8.0).cos());
    }

    #[test]
    fn test_sin() {
        float_eq!(sin(0.0), 0.0_f64.sin());
        float_eq!(sin(1.0), 1.0_f64.sin());
        float_eq!(sin(PI), PI.sin());
        float_eq!(sin(PI * 8.0), (PI * 8.0).sin());
    }

    #[test]
    fn test_sinh() {
        for x in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5] {
            float_eq!(sinh(x), x.sinh());
        }
    }

    #[test]
    fn test_sinh_near_zero_const() {
        assert_eq!((const { sinh(1e-20) }), 1e-20);
        assert_eq!((const { sinh(-1e-20) }), -1e-20);
        assert_eq!((const { sinh(f64::from_bits(1)) }), f64::from_bits(1));
        assert_eq!((const { sinh(-0.0) }).to_bits(), (-0.0_f64).to_bits());
        assert_eq!((const { sinh(1e-6) }), 1.0000000000001666e-6);
        assert_eq!((const { sinh(-1e-6) }), -1.0000000000001666e-6);
    }

    #[test]
    fn test_sinh_large_and_special_values() {
        assert!((const { sinh(710.0) } / 710.0_f64.sinh() - 1.0).abs() <= f64::EPSILON);
        assert!((const { sinh(-710.0) } / (-710.0_f64).sinh() - 1.0).abs() <= f64::EPSILON);
        assert_eq!((const { sinh(711.0) }), f64::INFINITY);
        assert_eq!((const { sinh(f64::INFINITY) }), f64::INFINITY);
        assert_eq!((const { sinh(f64::NEG_INFINITY) }), f64::NEG_INFINITY);
        assert!(const { sinh(f64::NAN) }.is_nan());
    }

    #[test]
    fn test_cosh() {
        for x in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5] {
            float_eq!(cosh(x), x.cosh());
        }
    }

    #[test]
    fn test_ln() {
        // float_eq!(ln(0.01), 0.01_f64.ln());
        // float_eq!(ln(0.5), 0.5_f64.ln());
        float_eq!(ln(1.0), 1.0_f64.ln());
        float_eq!(ln(2.0), 2.0_f64.ln());
        float_eq!(ln(10.0), 10.0_f64.ln());
        float_eq!(ln(1_000.0), 1_000.0_f64.ln());
    }
}
