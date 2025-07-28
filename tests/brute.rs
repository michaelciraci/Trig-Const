use core::f64::consts::PI;

use trig_const::{
    acos, acosh, asin, asinh, atan, atan2, atanh, cos, cosh, cot, csc, ln, sec, sin, sinh, tan,
};

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
        float_eq!($lhs, $rhs, 0.0000000001);
    };
}

#[test]
fn test_sin() {
    for x in float_loop(-8.0 * PI, 8.0 * PI, 0.1) {
        float_eq!(sin(x), x.sin());
    }
}

#[test]
fn test_cos() {
    for x in float_loop(-8.0 * PI, 8.0 * PI, 0.1) {
        float_eq!(cos(x), x.cos());
    }
}

#[test]
fn test_tan() {
    for x in float_loop(-8.0 * PI, 8.0 * PI, 0.1) {
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
fn test_ln() {
    for x in float_loop(0.01, 1000.0, 0.1) {
        float_eq!(ln(x), x.ln());
    }
}

#[test]
fn test_atan() {
    for x in float_loop(-2.0 * PI, 2.0 * PI, 0.1) {
        float_eq!(atan(x), x.atan());
    }
}

#[test]
fn test_atan2() {
    for x in float_loop(-2.0 * PI, 2.0 * PI, 0.1) {
        for y in float_loop(-2.0 * PI, 2.0 * PI, 0.1) {
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
fn test_atanh() {
    for x in float_loop(-0.99, 0.99, 0.01) {
        float_eq!(atanh(x), x.atanh());
    }
}

#[test]
fn test_sinh() {
    // Test over a moderate range to avoid huge numbers.
    for x in float_loop(-5.0, 5.0, 0.1) {
        float_eq!(sinh(x), x.sinh());
    }
}

#[test]
fn test_cosh() {
    // Test over a moderate range to avoid huge numbers.
    for x in float_loop(-5.0, 5.0, 0.1) {
        float_eq!(cosh(x), x.cosh());
    }
}

#[test]
fn test_sec() {
    for x in float_loop(-2.0 * PI, 2.0 * PI, 0.1) {
        float_eq!(sec(x), 1.0 / x.cos());
    }
}

#[test]
fn test_cot() {
    for x in float_loop(-2.0 * PI, 2.0 * PI, 0.1) {
        let our_cot = cot(x);
        // Test the stable identity: cot(x) * sin(x) = cos(x)
        if our_cot.is_finite() {
            let std_sin = x.sin();
            let std_cos = x.cos();
            float_eq!(our_cot * std_sin, std_cos);
        }
    }
}

#[test]
fn test_csc() {
    for x in float_loop(-2.0 * PI, 2.0 * PI, 0.1) {
        let our_csc = csc(x);
        // Test the stable identity: csc(x) * sin(x) = 1
        if our_csc.is_finite() {
            let std_sin = x.sin();
            float_eq!(our_csc * std_sin, 1.0_f64);
        }
    }
}
