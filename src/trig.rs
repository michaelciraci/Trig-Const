//! Accuracy-first sine/cosine using integer fixed-point intervals.
//!
//! Both argument reduction and evaluation retain lower and upper bounds.
//! Precision increases until the bounds agree on rounding, up to 512 bits.
//! An unresolved final interval supplies a best estimate instead of panicking.
//! See docs/trigonometry.md for the arithmetic and error accounting.

use crate::trig_coefficients::{COS_COEFFS, SIN_COEFFS};
use crate::trig_constants::{INV_PIO2, PIO2};
use crate::uint::Uint;

// N-1 fractional limbs and one integer limb, little-endian. The working
// precisions are 128, 256, and 512 fractional bits (N = 3, 5, and 9).
#[derive(Clone, Copy)]
struct Fixed<const N: usize> {
    raw: Uint<N>,
}

impl<const N: usize> Fixed<N> {
    const ZERO: Self = Self { raw: Uint::ZERO };

    const fn from_limbs(limbs: [u64; N]) -> Self {
        Self {
            raw: Uint::from_limbs(limbs),
        }
    }

    const fn one() -> Self {
        Self {
            raw: Uint::power_of_two((N - 1) * 64),
        }
    }

    const fn units(n: u64) -> Self {
        Self {
            raw: Uint::from_u64(n),
        }
    }

    const fn less(self, rhs: Self) -> bool {
        self.raw.less(rhs.raw)
    }

    const fn add(self, rhs: Self) -> Self {
        Self {
            raw: self.raw.add(rhs.raw),
        }
    }

    const fn sub(self, rhs: Self) -> Self {
        Self {
            raw: self.raw.sub(rhs.raw),
        }
    }

    const fn saturating_sub(self, rhs: Self) -> Self {
        if self.less(rhs) {
            Self::ZERO
        } else {
            self.sub(rhs)
        }
    }

    // Full integer product, then discard N-1 fractional limbs. The discarded
    // part is less than one fixed-point unit; callers round upper bounds up.
    const fn mul(self, rhs: Self) -> Self {
        if N == 3 && self.raw.limbs[2] == 0 && rhs.raw.limbs[2] == 0 {
            // Common-path operands are fractions below one. Four 64-bit
            // products suffice for the high half of their 128x128 product.
            // Retain the low-product carry: omitting it would invalidate the
            // outward-rounding bound.
            let low = self.raw.limbs[0] as u128 * rhs.raw.limbs[0] as u128;
            let cross = self.raw.limbs[0] as u128 * rhs.raw.limbs[1] as u128 + (low >> 64);
            let other = self.raw.limbs[1] as u128 * rhs.raw.limbs[0] as u128 + cross as u64 as u128;
            let high = self.raw.limbs[1] as u128 * rhs.raw.limbs[1] as u128
                + (cross >> 64)
                + (other >> 64);
            let mut result = Self::ZERO;
            result.raw.limbs[0] = high as u64;
            result.raw.limbs[1] = (high >> 64) as u64;
            return result;
        }
        let mut product = [0_u64; 18];
        let mut i = 0;
        while i < N {
            let mut carry = 0_u128;
            let mut j = 0;
            while j < N {
                let p = self.raw.limbs[i] as u128 * rhs.raw.limbs[j] as u128
                    + product[i + j] as u128
                    + carry;
                product[i + j] = p as u64;
                carry = p >> 64;
                j += 1;
            }
            product[i + N] = carry as u64;
            i += 1;
        }
        assert!(product[2 * N - 1] == 0);
        let mut result = Self::ZERO;
        i = 0;
        while i < N {
            result.raw.limbs[i] = product[i + N - 1];
            i += 1;
        }
        result
    }

    const fn div(self, divisor: u64) -> Self {
        assert!(divisor <= u32::MAX as u64);
        Self {
            raw: self.raw.div_small(divisor as u32),
        }
    }

    const fn bit(self, bit: usize) -> bool {
        self.raw.bit(bit)
    }

    // One round-to-nearest, ties-to-even conversion, using only integer bits.
    // Nonzero values in our fixed-point formats are normal binary64 numbers.
    const fn rounded_bits(self) -> u64 {
        let mut words = N;
        while words > 0 && self.raw.limbs[words - 1] == 0 {
            words -= 1;
        }
        if words == 0 {
            return 0;
        }
        let top = words * 64 - self.raw.limbs[words - 1].leading_zeros() as usize;
        let exponent = top as i32 - 1 - ((N - 1) * 64) as i32;
        let mut significand = if top <= 53 {
            self.raw.limbs[0] << (53 - top)
        } else {
            self.raw.window(top - 53)
        };
        if top > 53 {
            let round_bit = top - 54;
            let word = round_bit / 64;
            let mask = (1_u64 << (round_bit % 64)) - 1;
            let mut sticky = self.raw.limbs[word] & mask != 0;
            let mut i = 0;
            while i < word {
                sticky |= self.raw.limbs[i] != 0;
                i += 1;
            }
            if self.bit(round_bit) && (sticky || significand & 1 != 0) {
                significand += 1;
            }
        }
        // Addition deliberately handles a significand carry into the exponent.
        (((exponent + 1023) as u64) << 52) + significand - (1_u64 << 52)
    }
}

#[derive(Clone, Copy)]
struct Interval<const N: usize> {
    lo: Fixed<N>,
    hi: Fixed<N>,
}

impl<const N: usize> Interval<N> {
    const fn exact(value: Fixed<N>) -> Self {
        Self {
            lo: value,
            hi: value,
        }
    }

    const fn mul(self, rhs: Self) -> Self {
        Self {
            lo: self.lo.mul(rhs.lo),
            hi: self.hi.mul(rhs.hi).add(Fixed::units(1)),
        }
    }

    const fn div(self, divisor: u64) -> Self {
        Self {
            lo: self.lo.div(divisor),
            hi: self.hi.div(divisor).add(Fixed::units(1)),
        }
    }
}

// Return |r|, the quadrant, and the sign of r in x = quadrant*pi/2 + r (mod 2*pi).
// If the interval crosses a quadrant boundary, return that boundary instead.
const fn reduce<const N: usize>(bits: u64) -> Result<(Interval<N>, u64, bool), u64> {
    let exponent = ((bits >> 52) & 0x7ff) as i32 - 1023;
    let mantissa = (bits & 0x000f_ffff_ffff_ffff) | (1_u64 << 52);
    let precision = (N - 1) * 64;
    if bits <= 0x3fe0_0000_0000_0000 {
        // x <= 1/2, exact conversion
        let shift = (precision as i32 + exponent - 52) as usize;
        let mut value = Fixed::ZERO;
        value.raw.limbs[shift / 64] = mantissa << (shift % 64);
        if shift % 64 != 0 {
            value.raw.limbs[shift / 64 + 1] = mantissa >> (64 - shift % 64);
        }
        return Ok((Interval::exact(value), 0, false));
    }

    let mut product = Uint::<37>::ZERO;
    let mut carry = 0_u128;
    let mut i = 0;
    while i < INV_PIO2.len() {
        carry += INV_PIO2[i] as u128 * mantissa as u128;
        product.limbs[i] = carry as u64;
        carry >>= 64;
        i += 1;
    }
    product.limbs[36] = carry as u64;
    let binary_point = (2304 + 52 - exponent) as usize;
    let mut quadrant = product.window(binary_point) & 3;
    let mut lo = Fixed::ZERO;
    i = 0;
    while i < N - 1 {
        lo.raw.limbs[i] = product.window(binary_point - precision + 64 * i);
        i += 1;
    }
    // Constant truncation contributes < 2^(1024-2304), and extracting the
    // fraction loses < one working unit. Two units enclose their sum at all
    // working precisions. Retry if the interval could cross a quadrant boundary.
    let hi = lo.add(Fixed::units(2));
    if hi.raw.limbs[N - 1] != 0 {
        return Err((quadrant + 1) & 3);
    }
    let reflected = lo.bit(precision - 1);
    let fraction = if reflected {
        quadrant = (quadrant + 1) & 3;
        Interval {
            lo: Fixed::one().sub(hi),
            hi: Fixed::one().sub(lo),
        }
    } else {
        Interval { lo, hi }
    };
    let mut pi_lo = Fixed::ZERO;
    i = 0;
    while i < N {
        pi_lo.raw.limbs[i] = PIO2[i + PIO2.len() - N];
        i += 1;
    }
    let pi = Interval {
        lo: pi_lo,
        hi: pi_lo.add(Fixed::units(1)),
    };
    Ok((fraction.mul(pi), quadrant, reflected))
}

const fn series<const N: usize>(argument: Interval<N>, sine: bool) -> Interval<N> {
    let square = argument.mul(argument);
    let mut term = if sine {
        argument
    } else {
        Interval::exact(Fixed::one())
    };
    let mut sum = term;
    let mut degree = if sine { 1 } else { 0 };
    let mut subtract = true;
    loop {
        degree += 2;
        term = term.mul(square).div(degree * (degree - 1));
        if !Fixed::units(2).less(term.hi) {
            // The alternating-series remainder is bounded by the next term.
            return Interval {
                // Both kernels are nonnegative on the reduced domain.
                lo: sum.lo.saturating_sub(term.hi),
                hi: sum.hi.add(term.hi),
            };
        }
        sum = if subtract {
            Interval {
                lo: sum.lo.saturating_sub(term.hi),
                hi: sum.hi.sub(term.lo),
            }
        } else {
            Interval {
                lo: sum.lo.add(term.lo),
                hi: sum.hi.add(term.hi),
            }
        };
        subtract = !subtract;
    }
}

// Division-free first pass. Reciprocal factorials are stored as certified
// integer bounds rather than rounded floats. Horner evaluation uses the same
// outward-rounded interval operations as the wider Taylor-series fallback.
const fn polynomial(argument: Interval<3>, sine: bool) -> Option<Interval<3>> {
    const fn coefficient(value: u128) -> Interval<3> {
        let lo = Fixed::from_limbs([value as u64, (value >> 64) as u64, 0]);
        Interval {
            lo,
            hi: lo.add(Fixed::units(1)),
        }
    }
    let coefficients = if sine { &SIN_COEFFS } else { &COS_COEFFS };
    let square = argument.mul(argument);
    let mut i = coefficients.len() - 1;
    let mut tail = coefficient(coefficients[i]);
    while i > 0 {
        i -= 1;
        let c = coefficient(coefficients[i]);
        let product = square.mul(tail);
        if c.lo.less(product.hi) {
            return None;
        }
        tail = Interval {
            lo: c.lo.sub(product.hi),
            hi: c.hi.sub(product.lo),
        };
    }
    let mut correction = square.mul(tail);
    let base = if sine {
        correction = argument.mul(correction);
        argument
    } else {
        Interval::exact(Fixed::one())
    };
    // For r < 0.79, the omitted r^35/35! or r^34/34! is < 2^-128.
    // trig_coefficients verifies this bound with compile-time integer assertions.
    let upper_correction = correction.hi.add(Fixed::units(1));
    if base.lo.less(upper_correction) {
        return None;
    }
    Some(Interval {
        lo: base.lo.sub(upper_correction),
        hi: base.hi.sub(correction.lo).add(Fixed::units(1)),
    })
}

const fn with_sign(magnitude: u64, bits: u64, sine: bool, quadrant: u64, reflected: bool) -> u64 {
    let kernel_sine = sine ^ (quadrant & 1 != 0);
    let negative = if sine {
        (quadrant >= 2) ^ (kernel_sine && reflected) ^ (bits >> 63 != 0)
    } else {
        (quadrant == 1 || quadrant == 2) ^ (kernel_sine && reflected)
    };
    magnitude | ((negative as u64) << 63)
}

// Ok contains a settled rounding; Err contains the best estimate to use if
// no higher precision is available. All arithmetic stays in fixed point.
const fn finish<const N: usize>(
    result: Interval<N>,
    bits: u64,
    sine: bool,
    quadrant: u64,
    reflected: bool,
) -> Result<u64, u64> {
    let lo = result.lo.rounded_bits();
    let hi = result.hi.rounded_bits();
    if lo == hi {
        Ok(with_sign(lo, bits, sine, quadrant, reflected))
    } else {
        let midpoint = result.lo.add(result.hi.sub(result.lo).div(2));
        Err(with_sign(
            midpoint.rounded_bits(),
            bits,
            sine,
            quadrant,
            reflected,
        ))
    }
}

const fn fast(bits: u64, sine: bool) -> Option<u64> {
    let (argument, quadrant, reflected) = match reduce::<3>(bits & 0x7fff_ffff_ffff_ffff) {
        Ok(value) => value,
        Err(_) => return None,
    };
    match polynomial(argument, sine ^ (quadrant & 1 != 0)) {
        Some(result) => match finish(result, bits, sine, quadrant, reflected) {
            Ok(value) => Some(value),
            Err(_) => None,
        },
        None => None,
    }
}

const fn evaluate<const N: usize>(bits: u64, sine: bool) -> Result<u64, u64> {
    let (argument, quadrant, reflected) = match reduce::<N>(bits & 0x7fff_ffff_ffff_ffff) {
        Ok(value) => value,
        Err(quadrant) => {
            // The fraction of pi/2 is within two working units of this boundary.
            // Retry first; at the precision limit use its exact sine/cosine.
            let magnitude = if sine ^ (quadrant & 1 != 0) {
                0
            } else {
                1.0_f64.to_bits()
            };
            return Err(with_sign(magnitude, bits, sine, quadrant, false));
        }
    };
    let kernel_sine = sine ^ (quadrant & 1 != 0);
    let result = series(argument, kernel_sine);
    finish(result, bits, sine, quadrant, reflected)
}

pub(crate) const fn sincos(x: f64, sine: bool) -> f64 {
    let bits = x.to_bits();
    let magnitude = bits & 0x7fff_ffff_ffff_ffff;
    if magnitude >= 0x7ff0_0000_0000_0000 {
        return f64::NAN;
    }
    if magnitude < 0x3e40_0000_0000_0000 {
        // |x| < 2^-27
        return if sine { x } else { 1.0 };
    }
    if let Some(result) = fast(bits, sine) {
        return f64::from_bits(result);
    }
    match evaluate::<5>(bits, sine) {
        Ok(result) => f64::from_bits(result),
        Err(_) => match evaluate::<9>(bits, sine) {
            Ok(result) | Err(result) => f64::from_bits(result),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::{evaluate, fast, finish, series, sincos, Fixed, Interval};
    use crate::{cos, sin};

    #[test]
    fn difficult_rounding_const() {
        // Hard rounding regressions. The latter two also require accurate
        // reduction of the original input, beyond a double-double kernel alone.
        assert_eq!((const { sin(0.27772473379937157) }), 0.2741682781713897);
        assert_eq!((const { sin(1.0306029048601308) }), 0.8576092201418881);
        assert_eq!((const { cos(1.0240812040791702) }), 0.5198839755256663);
        assert_eq!((const { sin(0.37239483228894443) }), 0.36384716049001037);
        assert_eq!((const { cos(0.6608777347288374) }), 0.7894537733004547);
        assert_eq!((const { sin(1e20) }), -0.6452512852657808);
        assert_eq!((const { cos(1e20) }), 0.7639704044417283);
        assert_eq!((const { sin(f64::MAX) }), 0.004961954789184062);
        assert_eq!((const { cos(f64::MAX) }), -0.9999876894265599);
    }

    #[test]
    fn special_values_const() {
        assert_eq!((const { sin(0.0).to_bits() }), 0);
        assert_eq!((const { sin(-0.0).to_bits() }), 1_u64 << 63);
        assert_eq!((const { cos(0.0) }), 1.0);
        assert_eq!((const { cos(-0.0) }), 1.0);
        assert_eq!((const { sin(f64::from_bits(1)).to_bits() }), 1);
        assert_eq!((const { cos(f64::from_bits(1)) }), 1.0);
        assert!((const { sin(f64::INFINITY).is_nan() }));
        assert!((const { cos(f64::NEG_INFINITY).is_nan() }));
        assert!((const { sin(f64::NAN).is_nan() }));
        assert!((const { cos(f64::NAN).is_nan() }));
        for x in [f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
            assert!(sincos(x, true).is_nan());
            assert!(sincos(x, false).is_nan());
        }
    }

    #[test]
    fn independent_reference_values() {
        use core::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

        // Previously verified with GNU MPFR, rounding to nearest-even. These
        // inline regressions are independent of this crate's own arithmetic.
        // The final four inputs are hard cases from Tables V/VI of Lefevre,
        // Ly, and Zimmermann (ARITH 2026), linked in docs/trigonometry.md.
        let cases: &[(f64, f64, f64)] = &[
            (0.27772473379937157, 0.2741682781713897, 0.9616817328225256),
            (1.0306029048601308, 0.8576092201418881, 0.5143018817072541),
            (1.0240812040791702, 0.8542368828326417, 0.5198839755256663),
            (0.37239483228894443, 0.36384716049001037, 0.9314586645704449),
            (0.6608777347288374, 0.6138100193233035, 0.7894537733004547),
            (1e6, -0.34999350217129294, 0.9367521275331447),
            (1e20, -0.6452512852657808, 0.7639704044417283),
            (1e+100, -0.3806377310050287, 0.9247242387519338),
            (f64::MAX, 0.004961954789184062, -0.9999876894265599),
            (FRAC_PI_4, 0.7071067811865475, 0.7071067811865476),
            (FRAC_PI_2, 1.0, 6.123233995736766e-17),
            (PI, 1.2246467991473532e-16, -1.0),
            (
                f64::from_bits(FRAC_PI_2.to_bits() - 1),
                1.0,
                2.83276944882399e-16,
            ),
            (
                f64::from_bits(FRAC_PI_2.to_bits() + 1),
                1.0,
                -1.6081226496766366e-16,
            ),
            (
                9.90254958294253e+40,
                0.26993362518960334,
                -0.9628789321571009,
            ),
            (
                9.294860169592323e+302,
                0.39741478245503786,
                0.9176390851997396,
            ),
            (
                1.9042023571864292e+276,
                -0.9949191153293739,
                -0.1006774748998799,
            ),
            (
                6.985204031187422e+85,
                0.9744229903242662,
                0.22472168548565838,
            ),
        ];
        for &(x, sine, cosine) in cases {
            // Call the private implementation even with nightly forwarding enabled.
            for (x, sine) in [(x, sine), (-x, -sine)] {
                assert_eq!(sincos(x, true), sine, "sin({x:?})");
                assert_eq!(sincos(x, false), cosine, "cos({x:?})");
            }
            assert_eq!(evaluate::<9>(x.to_bits(), true), Ok(sine.to_bits()));
            assert_eq!(evaluate::<9>(x.to_bits(), false), Ok(cosine.to_bits()));
        }
    }

    #[test]
    fn precisions_agree_at_boundaries() {
        use core::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

        // A consistency check; the reference table above is the accuracy oracle.
        for x in [7.450580596923828e-9, 0.5, FRAC_PI_4, FRAC_PI_2, PI] {
            for bits in (x.to_bits() - 1)..=(x.to_bits() + 1) {
                for sine in [true, false] {
                    assert_eq!(
                        Ok(sincos(f64::from_bits(bits), sine).to_bits()),
                        evaluate::<9>(bits, sine),
                        "input bits {bits:#x}, sine={sine}"
                    );
                }
            }
        }
    }

    #[test]
    fn unresolved_rounding_uses_the_midpoint() {
        let one = Fixed::<3>::one();
        let half_ulp = Fixed::from_limbs([0, 1 << 11, 0]); // 2^-53
        let midpoint = one.add(half_ulp);
        let interval = Interval {
            lo: midpoint.sub(Fixed::units(1)),
            hi: midpoint.add(Fixed::units(3)),
        };
        // The endpoints disagree; the interval's center rounds up.
        assert_eq!(
            finish(interval, 0, true, 0, false),
            Err(1.0_f64.to_bits() + 1)
        );
        assert_eq!(
            finish(interval, 1 << 63, true, 0, false),
            Err((-1.0_f64).to_bits() + 1)
        );
        assert_eq!(
            finish(Interval::exact(one), 0, true, 0, false),
            Ok(1.0_f64.to_bits())
        );
    }

    #[test]
    fn unresolved_reduction_uses_the_quadrant_boundary() {
        // Zero fractional bits deliberately force an unresolved reduction.
        assert_eq!(
            (const { evaluate::<1>(1.0_f64.to_bits(), true) }),
            Err(1.0_f64.to_bits())
        );
        assert_eq!(
            (const { evaluate::<1>(4.0_f64.to_bits(), true) }),
            Err((-1.0_f64).to_bits())
        );
        assert_eq!(
            (const { evaluate::<1>(1.0_f64.to_bits(), false) }),
            Err((-0.0_f64).to_bits())
        );
    }

    #[test]
    fn series_bounds_near_zero_do_not_underflow() {
        let result = series(
            Interval::<9> {
                lo: Fixed::ZERO,
                hi: Fixed::units(2),
            },
            true,
        );
        assert!(result.lo.raw.is_zero());
        assert!(!result.hi.less(Fixed::units(2)));
    }

    #[test]
    fn fast_path_retries_hard_inputs() {
        assert!(fast(0.5_f64.to_bits(), true).is_some());
        assert!(fast(0.5_f64.to_bits(), false).is_some());
        // MPFR references for actual hard cases that the 128-bit polynomial
        // cannot settle. The complete function must recompute at 256 bits.
        for (x, expected) in [
            (7.096235686398445e-8_f64, 7.096235686398438e-8_f64),
            (2.149119332890821e-8, 2.149119332890821e-8),
        ] {
            assert_eq!(fast(x.to_bits(), true), None);
            assert_eq!(sincos(x, true).to_bits(), expected.to_bits());
            assert_eq!(evaluate::<5>(x.to_bits(), true), Ok(expected.to_bits()));
        }
        assert_eq!((const { sin(7.096235686398445e-8) }), 7.096235686398438e-8);
    }

    #[test]
    fn fractional_product_preserves_low_carries() {
        let max = Fixed::from_limbs([u64::MAX, u64::MAX, 0]);
        assert_eq!(max.mul(max).raw.limbs, [u64::MAX - 1, u64::MAX, 0]);
        let half = Fixed::from_limbs([0, 1 << 63, 0]);
        assert_eq!(max.mul(half).raw.limbs, [u64::MAX, u64::MAX >> 1, 0]);
        let small = Fixed::from_limbs([0, 1, 0]); // 2^-64
        assert_eq!(small.mul(small).raw.limbs, [1, 0, 0]);
    }

    #[test]
    fn integer_conversion_rounds_ties_to_even() {
        let one = Fixed::<5>::one();
        let mut half_ulp = Fixed::<5>::ZERO;
        half_ulp.raw.limbs[3] = 1 << 11; // 2^-53 at scale 2^256
        let midpoint = one.add(half_ulp);
        assert_eq!(midpoint.rounded_bits(), 1.0_f64.to_bits());
        assert_eq!(
            midpoint.sub(Fixed::units(1)).rounded_bits(),
            1.0_f64.to_bits()
        );
        assert_eq!(
            midpoint.add(Fixed::units(1)).rounded_bits(),
            1.0_f64.to_bits() + 1
        );
        assert_eq!(
            midpoint.add(half_ulp).add(half_ulp).rounded_bits(),
            1.0_f64.to_bits() + 2
        );
        // Carry from an all-ones significand into the exponent at 1.0.
        half_ulp.raw.limbs[3] >>= 1;
        assert_eq!(one.sub(half_ulp).rounded_bits(), 1.0_f64.to_bits());
        assert_eq!(
            one.sub(half_ulp).sub(Fixed::units(1)).rounded_bits(),
            1.0_f64.to_bits() - 1
        );
    }
}
