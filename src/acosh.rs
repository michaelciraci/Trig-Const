use crate::ln;

use super::{log1p::log1p, sqrt};

const LN2: f64 = 0.693147180559945309417232121458176568; /* 0x3fe62e42,  0xfefa39ef*/

/// Inverse hyperbolic cosine (f64)
///
/// Calculates the inverse hyperbolic cosine of `x`.
/// Is defined as `log(x + sqrt(x*x-1))`.
/// `x` must be a number greater than or equal to 1.
///
/// ```
/// # use trig_const::acosh;
/// const ACOSH_1: f64 = acosh(1.0);
/// assert_eq!(ACOSH_1, 0.0);
/// ```
pub const fn acosh(x: f64) -> f64 {
    let u = x.to_bits();
    let e = ((u >> 52) as usize) & 0x7ff;

    /* x < 1 domain error is handled in the called functions */

    if e < 0x3ff + 1 {
        /* |x| < 2, up to 2ulp error in [1,1.125] */
        return log1p(x - 1.0 + sqrt((x - 1.0) * (x - 1.0) + 2.0 * (x - 1.0)));
    }
    if e < 0x3ff + 26 {
        /* |x| < 0x1p26 */
        return ln(2.0 * x - 1.0 / (x + sqrt(x * x - 1.0)));
    }
    /* |x| >= 0x1p26 or nan */
    ln(x) + LN2
}
