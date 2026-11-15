//! High/low logarithm constants, computed at compile time.

pub(crate) const LN_TWO: SplitLog = LN_TWO_BOUNDS.split_at_f32();
// Retain 32 significant bits for exact multiplication by small integer exponents.
pub(crate) const LN_TWO_TRUNCATED: SplitLog = LN_TWO_BOUNDS.split_truncated(32);
pub(crate) const LOG2_E: SplitLog = LOG2_E_BOUNDS.split_at_f32();
pub(crate) const LOG_SCALE: SplitLog = (Bounds {
    lo: 2 * (LOG2_E_BOUNDS.lo / 3),
    hi: 2 * (LOG2_E_BOUNDS.hi / 3 + 1),
})
.split_at_f32();
pub(crate) const LOG2_THREE_HALVES: SplitLog = (Bounds {
    lo: divide_scaled(LN_THREE_HALVES_BOUNDS.lo, LN_TWO_BOUNDS.hi),
    hi: divide_scaled(LN_THREE_HALVES_BOUNDS.hi, LN_TWO_BOUNDS.lo) + 1,
})
.split_at_f32();

const LN_TWO_BOUNDS: Bounds = Bounds::ln_ratio(3);
const LN_THREE_HALVES_BOUNDS: Bounds = Bounds::ln_ratio(5);
const LOG2_E_BOUNDS: Bounds = Bounds {
    lo: divide_scaled(SCALE, LN_TWO_BOUNDS.hi),
    hi: divide_scaled(SCALE, LN_TWO_BOUNDS.lo) + 1,
};

const SCALE: u128 = 1_u128 << (u128::BITS - 1);

/// Inclusive bounds on a value, in units of 2^-127.
#[derive(Clone, Copy)]
struct Bounds {
    lo: u128,
    hi: u128,
}

/// A logarithm constant represented as `hi + lo` for higher-precision arithmetic.
/// `hi` has reduced precision; `lo` corrects the rounding or truncation.
pub(crate) struct SplitLog {
    pub hi: f64,
    pub lo: f64,
}

// floor(numerator*SCALE/denominator), without a 256-bit product.
const fn divide_scaled(numerator: u128, denominator: u128) -> u128 {
    assert!(numerator <= SCALE && denominator > SCALE / 2 && denominator < SCALE);
    let mut quotient = (numerator / denominator) * SCALE;
    let mut remainder = numerator % denominator;
    let mut bit = u128::BITS - 1;
    // Long division produces one fractional quotient bit per iteration.
    while bit > 0 {
        bit -= 1;
        remainder *= 2;
        if remainder >= denominator {
            remainder -= denominator;
            quotient |= 1_u128 << bit;
        }
    }
    quotient
}

impl Bounds {
    // ln((q+1)/(q-1)) = 2*(1/q + 1/(3q^3) + ...).
    const fn ln_ratio(q: u128) -> Self {
        assert!(q >= 3 && q <= 5);
        let mut power = SCALE / q;
        let mut lo = 0;
        let mut hi = 0;
        let mut denominator = 1;
        loop {
            let term = power / denominator;
            if term == 0 {
                // The remaining positive series is below two fixed-point units.
                return Bounds {
                    lo: 2 * lo,
                    hi: 2 * (hi + 2),
                };
            }
            lo += term;
            hi += term + 1;
            power /= q * q;
            denominator += 2;
        }
    }

    const fn rounded(self) -> f64 {
        let result = self.lo as f64 / SCALE as f64;
        // Both endpoints must round alike to determine a unique result.
        assert!(result.to_bits() == (self.hi as f64 / SCALE as f64).to_bits());
        result
    }

    const fn split_at_f32(self) -> SplitLog {
        let hi = self.rounded() as f32 as f64;
        self.split_with_hi(hi)
    }

    const fn split_truncated(self, significant_bits: u32) -> SplitLog {
        assert!(significant_bits > 0 && significant_bits <= f64::MANTISSA_DIGITS);
        let clear_bits = f64::MANTISSA_DIGITS - significant_bits;
        let hi = f64::from_bits(self.rounded().to_bits() & (u64::MAX << clear_bits));
        self.split_with_hi(hi)
    }

    const fn split_with_hi(self, hi: f64) -> SplitLog {
        let scaled = (hi * SCALE as f64) as u128;
        let (lo, upper, negative) = if scaled <= self.lo {
            (self.lo - scaled, self.hi - scaled, false)
        } else {
            // An interval straddling hi would leave the correction's sign unknown.
            assert!(scaled >= self.hi);
            (scaled - self.hi, scaled - self.lo, true)
        };
        let tail = (Bounds { lo, hi: upper }).rounded();
        SplitLog {
            hi,
            lo: if negative { -tail } else { tail },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        divide_scaled, Bounds, LN_TWO, LN_TWO_TRUNCATED, LOG2_E, LOG2_THREE_HALVES, LOG_SCALE,
        SCALE,
    };

    #[test]
    fn bounds_for_ln_two() {
        // ln_ratio(3) calculates ln((3 + 1) / (3 - 1)) = ln(2).
        let bounds = Bounds::ln_ratio(3);

        // ln(2) lies between 1/2 and 3/4; bounds are scaled by SCALE.
        assert!(bounds.lo > SCALE / 2);
        assert!(bounds.lo < bounds.hi);
        assert!(bounds.hi < (SCALE / 4) * 3);
    }

    #[test]
    fn ln_ratio_five_returns_bounds_for_ln_three_halves() {
        // (q + 1)/(q - 1) = 6/4 = 3/2 when q = 5.
        // ln(3/2) = 2 * (1/5 + 1/(3*5^3) + ...).
        let bounds = Bounds::ln_ratio(5);

        // All terms are positive. A geometric upper bound is 5/12 < 1/2.
        // The interval is expressed in units of 1/SCALE.
        assert!(bounds.lo > 2 * (SCALE / 5));
        assert!(bounds.lo < bounds.hi);
        assert!(bounds.hi < SCALE / 2);
    }

    #[test]
    fn logarithm_bounds_enclose_high_precision_references() {
        // floor(ln((q+1)/(q-1)) * 2^127), checked with exact rational series bounds.
        for (q, reference) in [
            (3, 117_932_881_612_756_647_068_972_071_382_077_242_199_u128),
            (5, 68_986_313_345_450_186_732_775_020_948_158_169_859),
        ] {
            let bounds = Bounds::ln_ratio(q);
            assert!(bounds.lo <= reference);
            assert!(bounds.hi > reference);
        }
    }

    #[test]
    fn divide_scaled_keeps_the_result_in_fixed_point_units() {
        // (1/2) / (3/4) = 2/3. Inputs and output are multiplied by SCALE.
        let half = SCALE / 2;
        let three_quarters = 3 * (SCALE / 4);
        // SCALE leaves remainder 2 when divided by 3.
        let two_thirds_rounded_down = 2 * (SCALE / 3) + 1;
        assert_eq!(divide_scaled(half, three_quarters), two_thirds_rounded_down);
    }

    #[test]
    fn divide_scaled_handles_exact_and_repeating_quotients() {
        // With denominator 3/4, the last quotient is 4/3, rounded down.
        let denominator = 3 * (SCALE / 4);
        assert_eq!(divide_scaled(0, denominator), 0);
        assert_eq!(divide_scaled(denominator / 2, denominator), SCALE / 2);
        assert_eq!(divide_scaled(denominator, denominator), SCALE);
        assert_eq!(divide_scaled(SCALE, denominator), SCALE + SCALE / 3);
    }

    #[test]
    fn split_at_f32_preserves_positive_and_negative_corrections() {
        // Adjacent f32 values at 1 are 1 and 1 + 2^-23.
        // Choose values one quarter and three quarters of the way between them.
        const UNIT: u128 = SCALE / (1_u128 << 25);
        const BELOW_MIDPOINT: u128 = SCALE + UNIT;
        const ABOVE_MIDPOINT: u128 = SCALE + 3 * UNIT;
        const POSITIVE: super::SplitLog = Bounds {
            lo: BELOW_MIDPOINT,
            hi: BELOW_MIDPOINT,
        }
        .split_at_f32();
        const NEGATIVE: super::SplitLog = Bounds {
            lo: ABOVE_MIDPOINT,
            hi: ABOVE_MIDPOINT,
        }
        .split_at_f32();
        // The first rounds down, so its correction is positive.
        assert_eq!(POSITIVE.hi, 1.0);
        assert_eq!(POSITIVE.lo, 1.0 / (1_u64 << 25) as f64);
        assert_eq!(POSITIVE.hi + POSITIVE.lo, 1.0 + 1.0 / (1_u64 << 25) as f64);

        // The second rounds up, so its correction is negative.
        assert_eq!(NEGATIVE.hi, 1.0 + 1.0 / (1_u64 << 23) as f64);
        assert_eq!(NEGATIVE.lo, -POSITIVE.lo);
        assert_eq!(NEGATIVE.hi + NEGATIVE.lo, 1.0 + 3.0 / (1_u64 << 25) as f64);
    }

    #[test]
    fn split_truncated_retains_discarded_bits() {
        // 1 + 3/16 is 1.0011 in binary. Keeping three significant bits
        // gives hi = 1.00; lo keeps the discarded 0.0011 = 3/16.
        const VALUE: u128 = SCALE + 3 * (SCALE / 16);
        let bounds = Bounds {
            lo: VALUE,
            hi: VALUE,
        };
        let split = bounds.split_truncated(3);
        assert_eq!(split.hi, 1.0);
        assert_eq!(split.lo, 3.0 / 16.0);
        let full = bounds.split_truncated(f64::MANTISSA_DIGITS);
        assert_eq!(full.hi, 1.0 + 3.0 / 16.0);
        assert_eq!(full.lo, 0.0);
    }

    #[test]
    fn rounded_accepts_distinct_bounds_with_the_same_rounding() {
        // Both endpoints are close enough to 1 to round to 1 as f64.
        let bounds = Bounds {
            lo: SCALE + (1_u128 << 72),
            hi: SCALE + (1_u128 << 73),
        };
        assert_eq!(bounds.rounded(), 1.0);
    }

    #[test]
    #[should_panic]
    fn split_with_hi_rejects_bounds_straddling_hi() {
        // The correction's sign is unknown when the bounds straddle hi.
        Bounds {
            lo: SCALE - 1,
            hi: SCALE + 1,
        }
        .split_with_hi(1.0);
    }

    #[test]
    #[should_panic]
    fn rounded_rejects_ambiguous_bounds() {
        Bounds {
            lo: SCALE / 2,
            hi: SCALE,
        }
        .rounded();
    }

    #[test]
    #[should_panic]
    fn split_truncated_rejects_zero_precision() {
        Bounds {
            lo: SCALE,
            hi: SCALE,
        }
        .split_truncated(0);
    }

    #[test]
    fn logarithm_splits_preserve_original_bits() {
        assert_eq!(LN_TWO.hi.to_bits(), 6.93147182464599609375e-1_f64.to_bits());
        assert_eq!(
            LN_TWO.lo.to_bits(),
            (-1.90465429995776804525e-9_f64).to_bits()
        );

        assert_eq!(
            LN_TWO_TRUNCATED.hi.to_bits(),
            6.93147180369123816490e-1_f64.to_bits()
        );
        assert_eq!(
            LN_TWO_TRUNCATED.lo.to_bits(),
            1.90821492927058770002e-10_f64.to_bits()
        );

        assert_eq!(LOG2_E.hi.to_bits(), 1.44269502162933349609_f64.to_bits());
        assert_eq!(LOG2_E.lo.to_bits(), 1.92596299112661746887e-8_f64.to_bits());

        assert_eq!(
            LOG_SCALE.hi.to_bits(),
            9.61796700954437255859e-1_f64.to_bits()
        );
        assert_eq!(
            LOG_SCALE.lo.to_bits(),
            (-7.02846165095275826516e-9_f64).to_bits()
        );

        assert_eq!(
            LOG2_THREE_HALVES.hi.to_bits(),
            5.84962487220764160156e-1_f64.to_bits()
        );
        assert_eq!(
            LOG2_THREE_HALVES.lo.to_bits(),
            1.35003920212974897128e-8_f64.to_bits()
        );
    }
}
