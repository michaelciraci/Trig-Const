/// Sine of an angle in radians.
///
/// Const evaluation uses fixed-point intervals to determine the nearest `f64`.
/// With the `nightly` feature, runtime calls use the configured math library.
///
/// ```
/// # use trig_const::sin;
/// # use core::f64::consts::PI;
/// # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
/// const SIN_PI: f64 = sin(PI);
/// float_eq(SIN_PI, 0.0);
/// ```
pub const fn sin(x: f64) -> f64 {
    nightly_exp!(sin, sin_inner, x)
}

const fn sin_inner(x: f64) -> f64 {
    crate::trig::sincos(x, true)
}
