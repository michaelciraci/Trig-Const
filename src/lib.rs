#![doc = include_str!("../README.md")]
#![no_std]
#![forbid(unsafe_code)]
#![allow(clippy::excessive_precision)]

mod acos;
mod asin;
mod atan;
mod atan2;
mod cos;
mod exp;
mod floor;
mod k_cos;
mod k_sin;
pub(crate) mod k_tan;
mod ln;
mod pow;
mod rem_pio2;
mod rem_pio2_large;
pub(crate) mod scalbn;
mod sin;
mod tan;
pub use acos::acos;
pub use asin::asin;
pub use atan::atan;
pub use atan2::atan2;
pub use cos::cos;
pub use exp::exp;
pub use floor::floor;
pub use ln::ln;
pub use pow::pow;
pub use sin::sin;
pub use tan::tan;

/// Number of sum iterations for Taylor series
const TAYLOR_SERIES_SUMS: usize = 16;

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
/// const SEC_PI: f64 = sec(PI);
/// assert_eq!(SEC_PI, -1.0);
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
/// const CSC_PI_2: f64 = csc(PI / 2.0);
/// assert_eq!(CSC_PI_2, 1.0);
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

/// x^pow
pub const fn expi(x: f64, mut pow: isize) -> f64 {
    let mut o = 1.0;

    while pow > 0 {
        o *= x;
        pow -= 1;
    }
    while pow < 0 {
        o /= x;
        pow += 1;
    }

    o
}

/// Factorial (x!)
pub const fn factorial(mut x: f64) -> f64 {
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
pub const fn sqrt(x: f64) -> f64 {
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

pub const fn fabs(x: f64) -> f64 {
    if x > 0.0 {
        x
    } else {
        -x
    }
}

const fn with_set_high_word(f: f64, hi: u32) -> f64 {
    let mut tmp = f.to_bits();
    tmp &= 0x00000000_ffffffff;
    tmp |= (hi as u64) << 32;
    f64::from_bits(tmp)
}
const fn with_set_low_word(f: f64, lo: u32) -> f64 {
    let mut tmp = f.to_bits();
    tmp &= 0xffffffff_00000000;
    tmp |= lo as u64;
    f64::from_bits(tmp)
}
const fn get_high_word(x: f64) -> u32 {
    (x.to_bits() >> 32) as u32
}

const fn get_low_word(x: f64) -> u32 {
    x.to_bits() as u32
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
        // float_eq!(ln(0.01), 0.01_f64.ln());
        // float_eq!(ln(0.5), 0.5_f64.ln());
        float_eq!(ln(1.0), 1.0_f64.ln());
        float_eq!(ln(2.0), 2.0_f64.ln());
        float_eq!(ln(10.0), 10.0_f64.ln());
        float_eq!(ln(1_000.0), 1_000.0_f64.ln());
    }
}
