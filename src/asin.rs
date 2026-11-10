/* origin: FreeBSD /usr/src/lib/msun/src/e_asin.c */
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */
/* asin(x)
 * Method :
 *      Since  asin(x) = x + x^3/6 + x^5*3/40 + x^7*15/336 + ...
 *      we approximate asin(x) on [0,0.5] by
 *              asin(x) = x + x*x^2*R(x^2)
 *      where
 *              R(x^2) is a rational approximation of (asin(x)-x)/x^3
 *      and its remez error is bounded by
 *              |(asin(x)-x)/x^3 - R(x^2)| < 2^(-58.75)
 *
 *      For x in [0.5,1]
 *              asin(x) = pi/2-2*asin(sqrt((1-x)/2))
 *      Let y = (1-x), z = y/2, s := sqrt(z), and pio2_hi+pio2_lo=pi/2;
 *      then for x>0.98
 *              asin(x) = pi/2 - 2*(s+s*z*R(z))
 *                      = pio2_hi - (2*(s+s*z*R(z)) - pio2_lo)
 *      For x<=0.98, let pio4_hi = pio2_hi/2, then
 *              f = hi part of s;
 *              c = sqrt(z) - f = (z-f*f)/(s+f)         ...f+c=sqrt(z)
 *      and
 *              asin(x) = pi/2 - 2*(s+s*z*R(z))
 *                      = pio4_hi+(pio4-2s)-(2s*z*R(z)-pio2_lo)
 *                      = pio4_hi+(pio4-2f)-(2s*z*R(z)-(pio2_lo+2c))
 *
 * Special cases:
 *      if x is NaN, return x itself;
 *      if |x|>1, return NaN with invalid signal.
 *
 */

use super::{fabs, get_high_word, get_low_word, sqrt, with_set_low_word, PIO2_HI, PIO2_LO};

// Fitted for asin and acos: P(z)/Q(z) approximates (asin(x) - x)/x^3,
// where z = x^2 lies in [0, 1/4].
// P_Sn and Q_Sn multiply z^n in P and Q, with Q(0) = 1.
const P_S0: f64 = 1.0 / 6.0;
const P_S1: f64 = -3.25565818622400915405e-01;
const P_S2: f64 = 2.01212532134862925881e-01;
const P_S3: f64 = -4.00555345006794114027e-02;
const P_S4: f64 = 7.91534994289814532176e-04;
const P_S5: f64 = 3.47933107596021167570e-05;
const Q_S1: f64 = -2.40339491173441421878e+00;
const Q_S2: f64 = 2.02094576023350569471e+00;
const Q_S3: f64 = -6.88283971605453293030e-01;
const Q_S4: f64 = 7.70381505559019352791e-02;

pub(crate) const fn asin_correction(z: f64) -> f64 {
    let p = z * (P_S0 + z * (P_S1 + z * (P_S2 + z * (P_S3 + z * (P_S4 + z * P_S5)))));
    let q = 1.0 + z * (Q_S1 + z * (Q_S2 + z * (Q_S3 + z * Q_S4)));
    p / q
}

/// Arcsine
///
/// ```
/// # use trig_const::asin;
/// const ASIN_PI: f64 = asin(0.0);
/// assert_eq!(ASIN_PI, 0.0);
/// ```
#[must_use]
pub const fn asin(x: f64) -> f64 {
    nightly_exp!(asin, asin_inner, x)
}

const fn asin_inner(mut x: f64) -> f64 {
    const TWO_POW_NEG_120: f64 = crate::expi(2.0, -120);
    let hx: u32 = get_high_word(x);
    let ix: u32 = hx & 0x7fffffff;
    /* |x| >= 1 or nan */
    if ix >= 0x3ff00000 {
        let lx: u32 = get_low_word(x);
        if ((ix - 0x3ff00000) | lx) == 0 {
            /* asin(1) = +-pi/2 with inexact */
            return x * PIO2_HI + TWO_POW_NEG_120;
        } else {
            return f64::INFINITY;
        }
    }
    /* |x| < 0.5 */
    if ix < 0x3fe00000 {
        /* if 0x1p-1022 <= |x| < 0x1p-26, avoid raising underflow */
        if ix >= 0x00100000 && ix < 0x3e500000 {
            return x;
        } else {
            return x + x * asin_correction(x * x);
        }
    }
    /* 1 > |x| >= 0.5 */
    let z: f64 = (1.0 - fabs(x)) * 0.5;
    let s: f64 = sqrt(z);
    let r: f64 = asin_correction(z);
    if ix >= 0x3fef3333 {
        /* if |x| > 0.975 */
        x = PIO2_HI - (2. * (s + s * r) - PIO2_LO);
    } else {
        /* sqrt_hi + sqrt_lo = sqrt(z) */
        let sqrt_hi: f64 = with_set_low_word(s, 0);
        let sqrt_lo: f64 = (z - sqrt_hi * sqrt_hi) / (s + sqrt_hi);
        x = 0.5 * PIO2_HI
            - (2.0 * s * r - (PIO2_LO - 2.0 * sqrt_lo) - (0.5 * PIO2_HI - 2.0 * sqrt_hi));
    }
    if hx >> 31 != 0 {
        -x
    } else {
        x
    }
}
