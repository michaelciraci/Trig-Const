use core::f64::consts::PI;

use trig_const::{acos, asin, cos, sin, tan};

fn float_loop(start: f64, stop: f64, step: f64) -> impl Iterator<Item = f64> {
    core::iter::successors(Some(start), move |prev| {
        let next = prev + step;
        (next < stop).then_some(next)
    })
}

macro_rules! float_eq {
    ($lhs:expr, $rhs:expr) => {
        assert!(
            ($lhs - $rhs).abs() < 0.0000000001,
            "lhs: {}, rhs: {}",
            $lhs,
            $rhs
        );
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
        dbg!(x);
        float_eq!(asin(x), x.asin());
    }
}

#[test]
fn test_acos() {
    for x in float_loop(-1.0, 1.0, 0.01) {
        dbg!(x);
        float_eq!(acos(x), x.acos());
    }
}
