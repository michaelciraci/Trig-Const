/* origin: FreeBSD /usr/src/lib/msun/src/e_acos.c */
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
/* acos(x)
 * Method :
 *      acos(x)  = pi/2 - asin(x)
 *      acos(-x) = pi/2 + asin(x)
 * For |x|<=0.5
 *      acos(x) = pi/2 - (x + x*x^2*R(x^2))     (see asin.c)
 * For x>0.5
 *      acos(x) = pi/2 - (pi/2 - 2asin(sqrt((1-x)/2)))
 *              = 2asin(sqrt((1-x)/2))
 *              = 2s + 2s*z*R(z)        ...z=(1-x)/2, s=sqrt(z)
 *              = 2f + (2c + 2s*z*R(z))
 *     where f=hi part of s, and c = (z-f*f)/(s+f) is the correction term
 *     for f so that f+c ~ sqrt(z).
 * For x<-0.5
 *      acos(x) = pi - 2asin(sqrt((1-|x|)/2))
 *              = pi - 0.5*(s+s*z*R(z)), where z=(1-|x|)/2,s=sqrt(z)
 *
 * Special cases:
 *      if x is NaN, return x itself;
 *      if |x|>1, return NaN with invalid signal.
 *
 * Function needed: sqrt
 */

use super::{asin::asin_correction, sqrt, PIO2_HI, PIO2_LO};

/// Arccosine
///
/// ```
/// # use trig_const::acos;
/// # use core::f64::consts::PI;
/// const ACOS_1: f64 = acos(1.0);
/// assert_eq!(ACOS_1, 0.0);
/// ```
#[must_use]
pub const fn acos(x: f64) -> f64 {
    nightly_exp!(acos, acos_inner, x)
}

const fn acos_inner(x: f64) -> f64 {
    const TWO_POW_NEG_120: f64 = crate::expi(2.0, -120);
    let hx: u32 = (x.to_bits() >> 32) as u32;
    let ix: u32 = hx & 0x7fffffff;
    /* |x| >= 1 or nan */
    if ix >= 0x3ff00000 {
        let lx: u32 = x.to_bits() as u32;

        if ((ix - 0x3ff00000) | lx) == 0 {
            /* acos(1)=0, acos(-1)=pi */
            if (hx >> 31) != 0 {
                return 2.0 * PIO2_HI + TWO_POW_NEG_120;
            }
            return 0.0;
        }
        return f64::INFINITY;
    }
    /* |x| < 0.5 */
    if ix < 0x3fe00000 {
        if ix <= 0x3c600000 {
            /* |x| < 2**-57 */
            return PIO2_HI + TWO_POW_NEG_120;
        }
        return PIO2_HI - (x - (PIO2_LO - x * asin_correction(x * x)));
    }
    /* x < -0.5 */
    if (hx >> 31) != 0 {
        let z = (1.0 + x) * 0.5;
        let s = sqrt(z);
        let w = asin_correction(z) * s - PIO2_LO;
        return 2. * (PIO2_HI - (s + w));
    }
    /* x > 0.5 */
    let z = (1.0 - x) * 0.5;
    let s = sqrt(z);
    let df = crate::with_set_low_word(s, 0);

    let c: f64 = (z - df * df) / (s + df);
    let w = asin_correction(z) * s + c;
    2. * (df + w)
}
