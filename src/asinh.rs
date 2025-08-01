use crate::{ln, log1p::log1p, sqrt};

const LN2: f64 = 0.693147180559945309417232121458176568; /* 0x3fe62e42,  0xfefa39ef*/

/* asinh(x) = sign(x)*log(|x|+sqrt(x*x+1)) ~= x - x^3/6 + o(x^5) */
/// Inverse hyperbolic sine (f64)
///
/// Calculates the inverse hyperbolic sine of `x`.
/// Is defined as `sgn(x)*log(|x|+sqrt(x*x+1))`.
/// ```
/// # use trig_const::asinh;
/// const ASINH_1: f64 = asinh(0.0);
/// assert_eq!(ASINH_1, 0.0);
/// ```
pub const fn asinh(mut x: f64) -> f64 {
    let mut u = x.to_bits();
    let e = ((u >> 52) as usize) & 0x7ff;
    let sign = (u >> 63) != 0;

    /* |x| */
    u &= (!0) >> 1;
    x = f64::from_bits(u);

    if e >= 0x3ff + 26 {
        /* |x| >= 0x1p26 or inf or nan */
        x = ln(x) + LN2;
    } else if e > 0x3ff {
        /* |x| >= 2 */
        x = ln(2.0 * x + 1.0 / (sqrt(x * x + 1.0) + x));
    } else if e >= 0x3ff - 26 {
        /* |x| >= 0x1p-26, up to 1.6ulp error in [0.125,0.5] */
        x = log1p(x + x * x / (sqrt(x * x + 1.0) + 1.0));
    } else {
        /* |x| < 0x1p-26, raise inexact if x != 0 */
    }

    if sign {
        -x
    } else {
        x
    }
}
