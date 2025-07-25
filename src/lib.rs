//! ## trig-const
//!
//! Rust implementation of const trig functions.
//!
//! This is implemented using a 16-term Taylor series approximation of cosine.
//! Correctness is favored over speed, especially considering the main use case for
//! this crate is to expose trigonometric functions for compile time.
//!
//! The implemntation was largely inspired by the work of Dr. Austin Henley and Dr. Stephen Marz:
//!   - GitHub Repo: <https://github.com/AZHenley/cosine>
//!   - Article: <https://austinhenley.com/blog/cosine.html>
//!
//! The implementation carries forward the original MIT license contained in the GitHub repo above.
//!
//! ## Requirements
//!
//! This crate supports any compiler version back to rustc 1.85
//!
//! ```toml
//! [dependencies]
//! trig-const = "0"
//! ```
//!
//! ## Example
//!
//! ```
//! # use trig_const::cos;
//! # use core::f64::consts::PI;
//! # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
//! const COS_PI: f64 = cos(PI);
//! float_eq(COS_PI, -1.0);
//! ```

#![no_std]
use core::f64::consts::PI;

/// Cosine
///
/// ```
/// # use trig_const::cos;
/// # use core::f64::consts::PI;
/// # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
/// const COS_PI: f64 = cos(PI);
/// float_eq(COS_PI, -1.0);
/// ```
pub const fn cos(mut x: f64) -> f64 {
    // If value is large, fold into smaller value
    while x < -0.1 {
        x += 2.0 * PI;
    }
    while x > 2.0 * PI + 0.1 {
        x -= 2.0 * PI;
    }
    let div = (x / PI) as u32;
    x = x - (div as f64 * PI);
    let sign = if div % 2 != 0 { -1.0 } else { 1.0 };

    let mut result = 1.0;
    let mut inter = 1.0;
    let num = x * x;

    let mut i = 1;
    while i <= 16 {
        let comp = 2.0 * i as f64;
        let den = comp * (comp - 1.0);
        inter *= num / den;
        if i % 2 == 0 {
            result += inter;
        } else {
            result -= inter;
        }
        i += 1;
    }

    sign * result
}

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
    cos(x - PI / 2.0)
}

/// Tangent
///
/// ```
/// # use trig_const::tan;
/// # use core::f64::consts::PI;
/// # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
/// const TAN_PI_4: f64 = tan(PI / 4.0);
/// float_eq(TAN_PI_4, 1.0);
/// ```
pub const fn tan(x: f64) -> f64 {
    sin(x) / cos(x)
}

/// Cotangent
///
/// ```
/// # use trig_const::cot;
/// # use core::f64::consts::PI;
/// # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
/// const COT_PI_4: f64 = cot(PI / 4.0);
/// float_eq(COT_PI_4, 1.0);
/// ```
pub const fn cot(x: f64) -> f64 {
    cos(x) / sin(x)
}

/// Secant
///
/// ```
/// # use trig_const::sec;
/// # use core::f64::consts::PI;
/// # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
/// let SEC_PI: f64 = sec(PI);
/// float_eq(SEC_PI, -1.0);
/// ```
pub const fn sec(x: f64) -> f64 {
    1.0 / cos(x)
}

/// Cosecant
///
/// ```
/// # use trig_const::csc;
/// # use core::f64::consts::PI;
/// # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
/// const CSC_PI_2: f64 = csc(PI / 2.0);
/// float_eq(CSC_PI_2, 1.0);
/// ```
pub const fn csc(x: f64) -> f64 {
    1.0 / sin(x)
}

#[cfg(test)]
mod tests {
    use core::f64::consts::PI;

    use crate::{cos, sin};

    macro_rules! float_eq {
        ($lhs:expr, $rhs:expr) => {
            assert!(($lhs - $rhs).abs() < 0.0001, "lhs: {}, rhs: {}", $lhs, $rhs);
        };
    }

    #[test]
    fn test_cos() {
        float_eq!(cos(0.0), 0.0_f64.cos());
        float_eq!(cos(1.0), 1.0_f64.cos());
        float_eq!(cos(PI), PI.cos());
        float_eq!(cos(PI * 8.0), (PI * 8.0).cos());
    }

    #[test]
    fn test_sin() {
        float_eq!(sin(0.0), 0.0_f64.sin());
        float_eq!(sin(1.0), 1.0_f64.sin());
        float_eq!(sin(PI), PI.sin());
        float_eq!(sin(PI * 8.0), (PI * 8.0).sin());
    }
}
