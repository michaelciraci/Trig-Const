//! High/low logarithm constants, computed at compile time.

pub(crate) const LN_TWO: SplitLog = split_at_f32(LN_TWO_BOUNDS);
// Retain 32 significant bits for exact multiplication by small integer exponents.
pub(crate) const LN_TWO_TRUNCATED: SplitLog = split_truncated(LN_TWO_BOUNDS, 32);
pub(crate) const LOG2_E: SplitLog = split_at_f32(LOG2_E_BOUNDS);
pub(crate) const LOG_SCALE: SplitLog = split_at_f32(Bounds {
    lo: 2 * (LOG2_E_BOUNDS.lo / 3),
    hi: 2 * (LOG2_E_BOUNDS.hi / 3 + 1),
});
pub(crate) const LOG2_THREE_HALVES: SplitLog = split_at_f32(Bounds {
    lo: divide_scaled(LN_THREE_HALVES_BOUNDS.lo, LN_TWO_BOUNDS.hi),
    hi: divide_scaled(LN_THREE_HALVES_BOUNDS.hi, LN_TWO_BOUNDS.lo) + 1,
});

const LN_TWO_BOUNDS: Bounds = ln_ratio(3);
const LN_THREE_HALVES_BOUNDS: Bounds = ln_ratio(5);
const LOG2_E_BOUNDS: Bounds = Bounds {
    lo: divide_scaled(SCALE, LN_TWO_BOUNDS.hi),
    hi: divide_scaled(SCALE, LN_TWO_BOUNDS.lo) + 1,
};

const SCALE: u128 = 1_u128 << (u128::BITS - 1);

#[derive(Clone, Copy)]
struct Bounds {
    lo: u128,
    hi: u128,
}

pub(crate) struct SplitLog {
    pub hi: f64,
    pub lo: f64,
}

// ln((q+1)/(q-1)) = 2*(1/q + 1/(3q^3) + ...).
const fn ln_ratio(q: u128) -> Bounds {
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

// floor(numerator*SCALE/denominator), without a 256-bit product.
const fn divide_scaled(numerator: u128, denominator: u128) -> u128 {
    assert!(numerator <= SCALE && denominator > SCALE / 2 && denominator < SCALE);
    let mut quotient = (numerator / denominator) * SCALE;
    let mut remainder = numerator % denominator;
    let mut bit = u128::BITS - 1;
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

const fn rounded(value: Bounds) -> f64 {
    let result = value.lo as f64 / SCALE as f64;
    assert!(result.to_bits() == (value.hi as f64 / SCALE as f64).to_bits());
    result
}

const fn split_at_f32(value: Bounds) -> SplitLog {
    let hi = rounded(value) as f32 as f64;
    split_with_hi(value, hi)
}

const fn split_truncated(value: Bounds, significant_bits: u32) -> SplitLog {
    assert!(significant_bits > 0 && significant_bits <= f64::MANTISSA_DIGITS);
    let clear_bits = f64::MANTISSA_DIGITS - significant_bits;
    let hi = f64::from_bits(rounded(value).to_bits() & (u64::MAX << clear_bits));
    split_with_hi(value, hi)
}

const fn split_with_hi(value: Bounds, hi: f64) -> SplitLog {
    let scaled = (hi * SCALE as f64) as u128;
    let (lo, upper, negative) = if scaled <= value.lo {
        (value.lo - scaled, value.hi - scaled, false)
    } else {
        assert!(scaled >= value.hi);
        (scaled - value.hi, scaled - value.lo, true)
    };
    let tail = rounded(Bounds { lo, hi: upper });
    SplitLog {
        hi,
        lo: if negative { -tail } else { tail },
    }
}

#[cfg(test)]
mod tests {
    use super::{ln_ratio, LN_TWO, LN_TWO_TRUNCATED, LOG2_E, LOG2_THREE_HALVES, LOG_SCALE, SCALE};

    #[test]
    fn bounds_for_ln_two() {
        // ln_ratio(3) calculates ln((3 + 1) / (3 - 1)) = ln(2).
        let bounds = ln_ratio(3);

        // ln(2) lies between 1/2 and 3/4; bounds are scaled by SCALE.
        assert!(bounds.lo > SCALE / 2);
        assert!(bounds.lo < bounds.hi);
        assert!(bounds.hi < (SCALE / 4) * 3);
    }

    #[test]
    fn logarithm_splits_preserve_original_bits() {
        for (value, hi, lo) in [
            (
                LN_TWO,
                6.93147182464599609375e-1_f64,
                -1.90465429995776804525e-9_f64,
            ),
            (
                LN_TWO_TRUNCATED,
                6.93147180369123816490e-1_f64,
                1.90821492927058770002e-10_f64,
            ),
            (LOG2_E, 1.44269502162933349609, 1.92596299112661746887e-8),
            (
                LOG_SCALE,
                9.61796700954437255859e-1,
                -7.02846165095275826516e-9,
            ),
            (
                LOG2_THREE_HALVES,
                5.84962487220764160156e-1,
                1.35003920212974897128e-8,
            ),
        ] {
            assert_eq!(value.hi.to_bits(), hi.to_bits());
            assert_eq!(value.lo.to_bits(), lo.to_bits());
        }
    }
}
