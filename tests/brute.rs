use core::f64::consts::PI;

use trig_const::{acos, acosh, asin, asinh, atan, atan2, cos, ln, pow, sin, sqrt, tan};

fn float_loop(start: f64, stop: f64, step: f64) -> impl Iterator<Item = f64> {
    core::iter::successors(Some(start), move |prev| {
        let next = prev + step;
        (next < stop).then_some(next)
    })
}

macro_rules! float_eq {
    ($lhs:expr, $rhs:expr, $tol:expr) => {
        if !$lhs.is_nan() && !$rhs.is_nan() {
            assert!(($lhs - $rhs).abs() < $tol, "lhs: {}, rhs: {}", $lhs, $rhs);
        }
    };
    ($lhs:expr, $rhs:expr) => {
        float_eq!($lhs, $rhs, 10e-12);
    };
}

#[test]
fn test_sin() {
    for x in float_loop(-8.0 * PI, 8.0 * PI, 0.01) {
        float_eq!(sin(x), x.sin());
    }
}

#[test]
fn test_cos() {
    for x in float_loop(-8.0 * PI, 8.0 * PI, 0.01) {
        float_eq!(cos(x), x.cos());
    }
}

#[test]
fn test_tan() {
    for x in float_loop(-8.0 * PI, 8.0 * PI, 0.01) {
        float_eq!(tan(x), x.tan());
    }
}

#[test]
fn test_asin() {
    for x in float_loop(-1.0, 1.0, 0.01) {
        float_eq!(asin(x), x.asin());
    }
}

#[test]
fn test_acos() {
    for x in float_loop(-1.0, 1.0, 0.01) {
        float_eq!(acos(x), x.acos());
    }
}

#[test]
fn test_atan() {
    for x in float_loop(-2.0 * PI, 2.0 * PI, 0.01) {
        float_eq!(atan(x), x.atan());
    }
}

#[test]
fn test_atan2() {
    for x in float_loop(-2.0 * PI, 2.0 * PI, 0.01) {
        for y in float_loop(-2.0 * PI, 2.0 * PI, 0.01) {
            float_eq!(atan2(x, y), x.atan2(y));
        }
    }
}

#[test]
fn test_asinh() {
    for x in float_loop(0.0, 1.0, 0.01) {
        float_eq!(asinh(x), x.asinh());
    }
}

#[test]
fn test_acosh() {
    for x in float_loop(0.0, 1.0, 0.01) {
        float_eq!(acosh(x), x.acosh());
    }
}

#[test]
fn test_ln() {
    for x in float_loop(1.0, 10.0, 0.01) {
        float_eq!(ln(x), x.ln());
    }
}

#[test]
fn test_sqrt() {
    for x in float_loop(0.0, 10.0, 0.01) {
        float_eq!(sqrt(x), x.sqrt());
    }
}

#[test]
fn test_pow() {
    for x in float_loop(-10.0, 10.0, 1.0) {
        for y in float_loop(-2.0, 3.0, 0.5) {
            if x != 0.0 {
                float_eq!(pow(x, y), x.powf(y));
            }
        }
    }
}
