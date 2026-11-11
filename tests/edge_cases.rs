use core::f64::consts::{FRAC_PI_2, PI};
use trig_const::{asinh, atanh, cos, exp, ln, pow, sin, tan};

#[test]
fn inverse_hyperbolic_special_values() {
    assert_eq!((const { asinh(0.0) }).to_bits(), 0.0_f64.to_bits());
    assert_eq!((const { asinh(-0.0) }).to_bits(), (-0.0_f64).to_bits());
    assert_eq!((const { asinh(f64::INFINITY) }), f64::INFINITY);
    assert_eq!((const { asinh(f64::NEG_INFINITY) }), f64::NEG_INFINITY);
    assert!(const { asinh(f64::NAN) }.is_nan());

    assert_eq!((const { atanh(0.0) }).to_bits(), 0.0_f64.to_bits());
    assert_eq!((const { atanh(-0.0) }).to_bits(), (-0.0_f64).to_bits());
    assert_eq!((const { atanh(1.0) }), f64::INFINITY);
    assert_eq!((const { atanh(-1.0) }), f64::NEG_INFINITY);
    assert!(const { atanh(1.1) }.is_nan());
    assert!(const { atanh(-1.1) }.is_nan());
    assert!(const { atanh(f64::INFINITY) }.is_nan());
    assert!(const { atanh(f64::NAN) }.is_nan());
}

#[test]
fn exp_and_ln_special_values() {
    assert_eq!((const { exp(0.0) }), 1.0);
    assert_eq!((const { exp(-0.0) }), 1.0);
    assert_eq!((const { exp(f64::INFINITY) }), f64::INFINITY);
    assert_eq!((const { exp(f64::NEG_INFINITY) }), 0.0);
    assert!(const { exp(f64::NAN) }.is_nan());

    assert_eq!((const { ln(0.0) }), f64::NEG_INFINITY);
    assert_eq!((const { ln(-0.0) }), f64::NEG_INFINITY);
    assert_eq!((const { ln(f64::INFINITY) }), f64::INFINITY);
    assert!(const { ln(-1.0) }.is_nan());
    assert!(const { ln(f64::NEG_INFINITY) }.is_nan());
    assert!(const { ln(f64::NAN) }.is_nan());
}

#[test]
fn pow_special_values() {
    assert_eq!((const { pow(f64::NAN, 0.0) }), 1.0);
    assert_eq!((const { pow(1.0, f64::NAN) }), 1.0);
    assert!(const { pow(-2.0, 0.5) }.is_nan());
    assert_eq!((const { pow(-0.0, 3.0) }).to_bits(), (-0.0_f64).to_bits());
    assert_eq!((const { pow(-0.0, -3.0) }), f64::NEG_INFINITY);
    assert_eq!((const { pow(2.0, f64::INFINITY) }), f64::INFINITY);
    assert_eq!((const { pow(2.0, f64::NEG_INFINITY) }), 0.0);
}

#[test]
fn trig_argument_reduction_const() {
    fn assert_close(actual: f64, expected: f64) {
        assert!((actual - expected).abs() <= 2e-15 * expected.abs().max(1.0));
    }

    assert_close(const { sin(FRAC_PI_2 - 1e-12) }, (FRAC_PI_2 - 1e-12).sin());
    assert_close(const { cos(FRAC_PI_2 + 1e-12) }, (FRAC_PI_2 + 1e-12).cos());
    assert_close(const { tan(FRAC_PI_2 - 1e-12) }, (FRAC_PI_2 - 1e-12).tan());
    assert_close(const { sin(-PI - 1e-12) }, (-PI - 1e-12).sin());
    assert_close(const { sin(1e20) }, 1e20_f64.sin());
    assert_close(const { cos(-1e20) }, (-1e20_f64).cos());
    assert_close(const { tan(f64::MAX) }, f64::MAX.tan());
}
