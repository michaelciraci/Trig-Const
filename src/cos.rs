/// Cosine of an angle in radians.
///
/// Const evaluation uses fixed-point intervals to determine the nearest `f64`.
/// With the `nightly` feature, runtime calls use the configured math library.
///
/// ```
/// # use trig_const::cos;
/// # use core::f64::consts::PI;
/// const COS_PI: f64 = cos(PI);
/// assert_eq!(COS_PI, -1.0);
/// ```
pub const fn cos(x: f64) -> f64 {
    nightly_exp!(cos, cos_inner, x)
}

const fn cos_inner(x: f64) -> f64 {
    crate::trig::sincos(x, false)
}
