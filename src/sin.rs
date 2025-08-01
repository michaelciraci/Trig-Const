// origin: FreeBSD /usr/src/lib/msun/src/s_sin.c */
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================

// sin(x)
// Return sine function of x.
//
// kernel function:
//      k_sin            ... sine function on [-pi/4,pi/4]
//      k_cos            ... cose function on [-pi/4,pi/4]
//      rem_pio2         ... argument reduction routine
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

use crate::{k_cos::k_cos, k_sin::k_sin, rem_pio2::rem_pio2};

/// Sine
///
/// ```
/// # use trig_const::sin;
/// # use core::f64::consts::PI;
/// # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
/// const SIN_PI: f64 = sin(PI);
/// float_eq(SIN_PI, 0.0);
/// ```
pub const fn sin(x: f64) -> f64 {
    /* High word of x. */
    let ix = (f64::to_bits(x) >> 32) as u32 & 0x7fffffff;

    /* |x| ~< pi/4 */
    if ix <= 0x3fe921fb {
        if ix < 0x3e500000 {
            /* |x| < 2**-26 */
            /* raise inexact if x != 0 and underflow if subnormal*/
            return x;
        }
        return k_sin(x, 0.0, 0);
    }

    /* sin(Inf or NaN) is NaN */
    if ix >= 0x7ff00000 {
        return f64::NAN;
    }

    /* argument reduction needed */
    let (n, y0, y1) = rem_pio2(x);
    match n & 3 {
        0 => k_sin(y0, y1, 1),
        1 => k_cos(y0, y1),
        2 => -k_sin(y0, y1, 1),
        _ => -k_cos(y0, y1),
    }
}
