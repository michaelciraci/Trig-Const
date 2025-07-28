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
//! ## Examples
//!
//! ```
//! # use trig_const::cos;
//! # use core::f64::consts::PI;
//! # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
//! const COS_PI: f64 = cos(PI);
//! float_eq(COS_PI, -1.0);
//! ```
//!
//! ```
//! use std::f64::consts::PI;
//! use trig_const::{atan2, sin};
//!
//! const SIN_PI_4: f64 = sin(PI / 2.0);
//! const ATAN2_0_0: f64 = atan2(0.0, 0.0);
//!
//! fn main() {
//!     println!("{SIN_PI_4}\n{ATAN2_0_0}");
//! }
//! ```
//!

#![no_std]
#![forbid(unsafe_code)]

use core::f64::{
    self,
    consts::{FRAC_PI_2, PI},
};

// A const-compatible floor function.
const fn floor(x: f64) -> f64 {
    let i = x as i64;
    let f = i as f64;
    if f == x || x >= 0.0 {
        f
    } else {
        // For negative x, casting truncates towards zero, so we subtract 1
        // to get the correct floor.
        f - 1.0
    }
}

/// Number of sum iterations for Taylor series
const TAYLOR_SERIES_SUMS: usize = 16;

/// Cosine
///
/// ```
/// # use trig_const::cos;
/// # use core::f64::consts::PI;
/// # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
/// const COS_PI: f64 = cos(PI);
/// float_eq(COS_PI, -1.0);
/// ```
pub const fn cos(x: f64) -> f64 {
    // Handle edge cases. cos(inf) is undefined.
    if x.is_nan() || x.is_infinite() {
        return f64::NAN;
    }

    const fn taylor_cos(x: f64) -> f64 {
        let x_squared = x * x;
        let mut s = 1.0;
        let mut term = 1.0;
        let mut n = 1;
        while n <= TAYLOR_SERIES_SUMS {
            let comp = 2.0 * n as f64;
            term *= -x_squared / (comp * (comp - 1.0));
            s += term;
            n += 1;
        }
        s
    }

    const FRAC_PI_2: f64 = core::f64::consts::FRAC_PI_2;
    const PI: f64 = core::f64::consts::PI;
    const TWO_PI: f64 = 2.0 * PI;

    let x = x - floor(x / TWO_PI) * TWO_PI;

    // Corrected Quadrant Reduction
    if x < FRAC_PI_2 {
        // Quadrant I: [0, PI/2)
        taylor_cos(x)
    } else if x < PI {
        // Quadrant II: [PI/2, PI)
        -taylor_cos(PI - x)
    } else if x < 1.5 * PI {
        // Quadrant III: [PI, 1.5*PI)
        -taylor_cos(x - PI)
    } else {
        // Quadrant IV: [1.5*PI, 2*PI)
        taylor_cos(TWO_PI - x)
    }
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
    if x.is_nan() || x.is_infinite() {
        return f64::NAN;
    }

    const fn floor(x: f64) -> f64 {
        let i = x as i64;
        let f = i as f64;
        if f == x || x >= 0.0 {
            f
        } else {
            f - 1.0
        }
    }

    // Taylor series for sin(x) = x - x^3/3! + x^5/5! - ...
    // This is only called with a small input in [0, PI/2].
    const fn taylor_sin(x: f64) -> f64 {
        let x_squared = x * x;
        let mut s = x;
        let mut term = x;
        let mut n = 1;
        while n <= TAYLOR_SERIES_SUMS {
            let comp = 2.0 * n as f64 + 1.0;
            term *= -x_squared / (comp * (comp - 1.0));
            s += term;
            n += 1;
        }
        s
    }

    const PI: f64 = core::f64::consts::PI;
    const TWO_PI: f64 = 2.0 * PI;
    const FRAC_PI_2: f64 = core::f64::consts::FRAC_PI_2;

    // Perform a correct Euclidean-style remainder to map x to [0, 2*PI).
    let x = x - floor(x / TWO_PI) * TWO_PI;

    // Reduce to [0, PI/2] and apply the correct sign based on quadrant.
    if x < FRAC_PI_2 {
        // Quadrant I: [0, PI/2)
        taylor_sin(x)
    } else if x < PI {
        // Quadrant II: [PI/2, PI) -> sin(x) = sin(PI - x)
        taylor_sin(PI - x)
    } else if x < 1.5 * PI {
        // Quadrant III: [PI, 1.5*PI) -> sin(x) = -sin(x - PI)
        -taylor_sin(x - PI)
    } else {
        // Quadrant IV: [1.5*PI, 2*PI) -> sin(x) = -sin(2*PI - x)
        -taylor_sin(TWO_PI - x)
    }
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
    let sin_calc = sin(x);
    if sin_calc == 0.0 {
        f64::INFINITY
    } else {
        cos(x) / sin_calc
    }
}

/// Secant
///
/// ```
/// # use trig_const::sec;
/// # use core::f64::consts::PI;
/// # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
/// const SEC_PI: f64 = sec(PI);
/// float_eq(SEC_PI, -1.0);
/// ```
pub const fn sec(x: f64) -> f64 {
    let cos_calc = cos(x);
    if cos_calc == 0.0 {
        f64::INFINITY
    } else {
        1.0 / cos_calc
    }
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
    let sin_calc = sin(x);
    if sin_calc == 0.0 {
        f64::INFINITY
    } else {
        1.0 / sin_calc
    }
}

/// Hyperbolic Sine
///
/// ```
/// # use trig_const::sinh;
/// const SINH_0: f64 = sinh(0.0);
/// assert_eq!(SINH_0, 0.0);
/// ```
pub const fn sinh(x: f64) -> f64 {
    (exp(x) - exp(-x)) / 2.0
}

/// Hyperbolic Cosine
///
/// ```
/// # use trig_const::cosh;
/// const COSH_0: f64 = cosh(0.0);
/// assert_eq!(COSH_0, 1.0);
/// ```
pub const fn cosh(x: f64) -> f64 {
    (exp(x) + exp(-x)) / 2.0
}

/// Arcsine
///
/// ```
/// # use trig_const::asin;
/// const ASIN_PI: f64 = asin(0.0);
/// assert_eq!(ASIN_PI, 0.0);
/// ```
pub const fn asin(x: f64) -> f64 {
    if x.is_infinite() || x.abs() > 1.0 {
        return f64::NAN;
    } else if x == 1.0 {
        return f64::consts::FRAC_PI_2;
    } else if x == -1.0 {
        return -f64::consts::FRAC_PI_2;
    }

    // As we start to get past 0.8, the number of summations needed for an accurate
    // Taylor series approximation starts to get unweidy. We can use the property
    // that arcsin(x) = pi/2 - 2*arcsin(sqrt((1 - x) / 2)) to reduce
    const RANGE_REDUCTION_THRESHOLD: f64 = 0.5;
    if x.abs() > RANGE_REDUCTION_THRESHOLD {
        let sign = x.signum();
        let abs_x = x.abs();

        let y = sqrt((1.0 - abs_x) / 2.0);
        return sign * (f64::consts::FRAC_PI_2 - 2.0 * asin(y));
    }

    let mut n = 1;
    let mut s = x;

    while n < TAYLOR_SERIES_SUMS {
        let numer1 = factorial(2.0 * n as f64);
        let numer2 = expi(x, 2 * n + 1);

        // Calculate all denom terms;
        let denom1 = expi(4.0, n);
        let denom2 = factorial(n as f64) * factorial(n as f64);
        let denom3 = 2.0 * n as f64 + 1.0;

        // Try to match terms to divide to stop number getting too large
        let f1 = numer1 / denom2;
        let f2 = numer2 / denom1;

        s += f1 * f2 / denom3;

        n += 1;
    }

    s
}

/// Arccosine
///
/// ```
/// # use trig_const::acos;
/// # use core::f64::consts::PI;
/// const ACOS_1: f64 = acos(1.0);
/// assert_eq!(ACOS_1, 0.0);
/// ```
pub const fn acos(x: f64) -> f64 {
    if x.is_infinite() || x.abs() > 1.0 {
        f64::NAN
    } else {
        f64::consts::FRAC_PI_2 - asin(x)
    }
}

/// Arctangent
///
/// ```
/// # use trig_const::atan;
/// # use core::f64::consts::PI;
/// # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
/// const ATAN_1: f64 = atan(1.0);
/// float_eq(ATAN_1, PI / 4.0);
/// ```
pub const fn atan(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    } else if x.is_infinite() {
        if x > 0.0 {
            return FRAC_PI_2;
        } else {
            return -FRAC_PI_2;
        }
    } else if x == 0.0 {
        return 0.0;
    }

    const fn atan_taylor_series(x: f64) -> f64 {
        let mut s = 0.0;
        let mut term = x;
        let mut sign = 1.0;
        let x_squared = x * x;

        let mut n = 0;
        while n < TAYLOR_SERIES_SUMS {
            let denom = (2 * n + 1) as f64;
            s += sign * term / denom;
            term *= x_squared;
            if sign == 1.0 {
                sign = -1.0;
            } else {
                sign = 1.0;
            }
            n += 1;
        }

        s
    }
    let sign = x.signum();
    let mut val = x.abs();
    if val > 1.0 {
        sign * (FRAC_PI_2 - atan(1.0 / val))
    } else {
        // Range reduction using the identity
        // atan(x) = 2 atan(x / (1 + sqrt(1 + x^2)))
        // until |x| \in [0, 0.125]
        let mut multiplicand = 1.0;
        while val > 0.125 {
            val /= 1.0 + sqrt(1.0 + val * val);
            multiplicand *= 2.0;
        }
        sign * multiplicand * atan_taylor_series(val)
    }
}

/// Arctan2
///
/// ```
/// # use trig_const::atan2;
/// # use core::f64::consts::PI;
/// # fn float_eq(lhs: f64, rhs: f64) { assert!((lhs - rhs).abs() < 0.0001, "lhs: {}, rhs: {}", lhs, rhs); }
/// const ATAN2_0_1: f64 = atan2(0.0, 1.0);
/// float_eq(ATAN2_0_1, 0.0);
/// ```
pub const fn atan2(y: f64, x: f64) -> f64 {
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }

    if x == 0.0 {
        if y > 0.0 {
            FRAC_PI_2
        } else if y < 0.0 {
            -FRAC_PI_2
        } else {
            0.0
        }
    } else if x > 0.0 {
        atan(y / x) // Quadrant I or IV
    } else if y >= 0.0 {
        atan(y / x) + PI // Quadrant II
    } else {
        atan(y / x) - PI // Quadrant III
    }
}

/// Inverse hyperbolic sine
///
/// ```
/// # use trig_const::asinh;
/// const ASINH_1: f64 = asinh(0.0);
/// assert_eq!(ASINH_1, 0.0);
/// ```
pub const fn asinh(x: f64) -> f64 {
    if x.is_nan() {
        f64::NAN
    } else if x.is_infinite() {
        x
    } else {
        ln(x + sqrt(x * x + 1.0))
    }
}

/// Inverse hyperbolic cosine
///
/// ```
/// # use trig_const::acosh;
/// const ACOSH_1: f64 = acosh(1.0);
/// assert_eq!(ACOSH_1, 0.0);
/// ```
pub const fn acosh(x: f64) -> f64 {
    if x.is_nan() {
        f64::NAN
    } else if x.is_infinite() {
        x
    } else if x < 1.0 {
        f64::NAN
    } else {
        ln(x + sqrt(x * x - 1.0))
    }
}

/// Inverse hyperbolic tangent using range reduction.
pub const fn atanh(x: f64) -> f64 {
    if x.is_nan() || x.abs() > 1.0 {
        return f64::NAN;
    } else if x == 1.0 {
        return f64::INFINITY;
    } else if x == -1.0 {
        return f64::NEG_INFINITY;
    } else if x == 0.0 {
        return 0.0;
    }

    let sign = x.signum();
    let mut y = x.abs();

    // Iteratively reduce `|x|` until it's in `[0,0.125]`
    // using the identity arctanh(x) = 2arctanh(x / (1 + sqrt(1 - x^2)))
    let mut multiplicand = 1.0;
    while y > 0.125 {
        y /= 1.0 + sqrt(1.0 - y * y);
        multiplicand *= 2.0;
    }

    // atanh(y) = y + y^3/3 + y^5/5 + ...
    let y_squared = y * y;
    let mut term = y;
    let mut s = y;
    let mut n = 1;
    while n < TAYLOR_SERIES_SUMS {
        term *= y_squared;
        s += term / ((2 * n + 1) as f64);
        n += 1;
    }
    sign * multiplicand * s
}

/// e^x using range reduction for high precision.
const fn exp(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    } else if x == f64::INFINITY {
        return f64::INFINITY;
    } else if x == f64::NEG_INFINITY {
        return 0.0;
    }

    const LN2: f64 = core::f64::consts::LN_2;
    let k = (x / LN2 + 0.5) as i64;

    // Add checks to handle overflow before bit manipulation.
    // The maximum exponent for a normal f64 corresponds to k=1023.
    // k=1024 results in INFINITY.
    if k >= 1024 {
        return f64::INFINITY;
    }
    // The minimum exponent for a normal f64 corresponds to roughly k=-1022.
    // Anything smaller will underflow to zero. -1075 is a safe threshold.
    if k <= -1075 {
        return 0.0;
    }

    let r = x - k as f64 * LN2;

    // Taylor series for e^r
    let mut s = 1.0;
    let mut term = 1.0;
    let mut n = 1;
    while n < 20 {
        term *= r / n as f64;
        s += term;
        n += 1;
    }

    // Calculate 2^k by directly constructing the f64's bit representation.
    // An f64's exponent is stored with a bias of 1023. We shift the
    // biased exponent left by 52 bits to position it correctly, leaving
    // the mantissa as zero, effectively creating the number 1.0 * 2^k.
    let two_k_bits = ((1023 + k) as u64) << 52;
    let scale = f64::from_bits(two_k_bits);

    scale * s
}

/// x^pow
const fn expi(x: f64, mut pow: usize) -> f64 {
    let mut o = 1.0;

    while pow > 0 {
        o *= x;
        pow -= 1;
    }

    o
}

/// Factorial (x!)
const fn factorial(mut x: f64) -> f64 {
    if x == 0.0 || x == 1.0 {
        1.0
    } else {
        let mut s = 1.0;
        while x > 1.0 {
            s *= x;
            x -= 1.0;
        }
        s
    }
}

/// Const sqrt function using Newton's method
const fn sqrt(x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    } else if x.is_infinite() || x == 0.0 {
        return x;
    }

    // Use Newton's method for sqrt calculation
    let mut current_guess = 1.0;

    let mut i = 0;
    while i < TAYLOR_SERIES_SUMS {
        current_guess = 0.5 * (current_guess + x / current_guess);
        i += 1;
    }

    current_guess
}

/// Computes natural log using the identity
/// `ln(x) = 2 arctanh((x - 1) / (x+1))`
pub const fn ln(x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    } else if x == 0.0 {
        return f64::NEG_INFINITY;
    } else if x == 1.0 {
        return 0.0;
    } else if x.is_infinite() {
        return f64::INFINITY;
    }
    2.0 * atanh((x - 1.0) / (x + 1.0))
}

#[cfg(test)]
mod tests {
    use core::f64::consts::{E, PI};

    use crate::{cos, cosh, exp, expi, factorial, ln, sin, sinh, sqrt};

    macro_rules! float_eq {
        ($lhs:expr, $rhs:expr) => {
            assert!(($lhs - $rhs).abs() < 0.0001, "lhs: {}, rhs: {}", $lhs, $rhs);
        };
    }

    #[test]
    fn test_factorial() {
        assert_eq!(factorial(0.0), 1.0);
        assert_eq!(factorial(1.0), 1.0);
        assert_eq!(factorial(2.0), 2.0);
        assert_eq!(factorial(3.0), 6.0);
        assert_eq!(factorial(4.0), 24.0);
        assert_eq!(factorial(5.0), 120.0);
    }

    #[test]
    fn test_expi() {
        assert_eq!(expi(2.0, 0), 1.0);
        assert_eq!(expi(2.0, 4), 16.0);
        assert_eq!(expi(2.0, 5), 32.0);
        assert_eq!(expi(3.0, 3), 27.0);
    }

    #[test]
    fn test_exp() {
        float_eq!(exp(0.0), 1.0);
        float_eq!(exp(1.0), E);
    }

    #[test]
    fn test_sqrt() {
        float_eq!(sqrt(0.0), 0.0);
        float_eq!(sqrt(1.0), 1.0);
        float_eq!(sqrt(4.0), 2.0);
        float_eq!(sqrt(9.0), 3.0);
        float_eq!(sqrt(16.0), 4.0);
        float_eq!(sqrt(25.0), 5.0);
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

    #[test]
    fn test_sinh() {
        for x in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5] {
            float_eq!(sinh(x), x.sinh());
        }
    }

    #[test]
    fn test_cosh() {
        for x in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5] {
            float_eq!(cosh(x), x.cosh());
        }
    }

    #[test]
    fn test_ln() {
        float_eq!(ln(0.01), 0.01_f64.ln());
        float_eq!(ln(0.5), 0.5_f64.ln());
        float_eq!(ln(1.0), 1.0_f64.ln());
        float_eq!(ln(2.0), 2.0_f64.ln());
        float_eq!(ln(10.0), 10.0_f64.ln());
        float_eq!(ln(1_000.0), 1_000.0_f64.ln());
    }
}
