// origin: FreeBSD /usr/src/lib/msun/src/s_cos.c */
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================

// cos(x)
// Return cosine function of x.
//
// kernel function:
//      k_sin           ... sine function on [-pi/4,pi/4]
//      k_cos           ... cosine function on [-pi/4,pi/4]
//      rem_pio2        ... argument reduction routine
//
// Method.
//      Let S,C and T denote the sin, cos and tan respectively on
//      [-PI/4, +PI/4]. Reduce the argument x to y1+y2 = x-k*pi/2
//      in [-pi/4 , +pi/4], and let n = k mod 4.
//      We have
//
//          n        sin(x)      cos(x)        tan(x)
//     ----------------------------------------------------------
//          0          S           C             T
//          1          C          -S            -1/T
//          2         -S          -C             T
//          3         -C           S            -1/T
//     ----------------------------------------------------------
//
// Special cases:
//      Let trig be any of sin, cos, or tan.
//      trig(+-INF)  is NaN, with signals;
//      trig(NaN)    is that NaN;
//
// Accuracy:
//      TRIG(x) returns trig(x) nearly rounded
//

use crate::{k_cos::k_cos, k_sin::k_sin, rem_pio2::rem_pio2};

/// Cosine
///
/// ```
/// # use trig_const::cos;
/// # use core::f64::consts::PI;
/// const COS_PI: f64 = cos(PI);
/// assert_eq!(COS_PI, -1.0);
/// ```
#[must_use]
pub const fn cos(x: f64) -> f64 {
    nightly_exp!(cos, cos_inner, x)
}

const fn cos_inner(x: f64) -> f64 {
    const SMALL_ANGLE: f64 = core::f64::consts::SQRT_2 * crate::expi(2.0, -27);
    let ix = crate::get_high_word(x.abs());

    /* |x| ~< pi/4 */
    if ix <= crate::get_high_word(core::f64::consts::FRAC_PI_4) {
        if ix < crate::get_high_word(SMALL_ANGLE) {
            /* if x < 2**-27 * sqrt(2) */
            /* raise inexact if x != 0 */
            if x as i32 == 0 {
                return 1.0;
            }
        }
        return k_cos(x, 0.0);
    }

    /* cos(Inf or NaN) is NaN */
    if ix >= crate::get_high_word(f64::INFINITY) {
        return f64::NAN;
    }

    /* argument reduction needed */
    let (n, y0, y1) = rem_pio2(x);
    match n & 3 {
        0 => k_cos(y0, y1),
        1 => -k_sin(y0, y1, 1),
        2 => -k_cos(y0, y1),
        _ => k_sin(y0, y1, 1),
    }
}
