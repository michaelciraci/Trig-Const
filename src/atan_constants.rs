//! Arctangent reference angles, computed with integer bounds at compile time.

pub(crate) const ATAN_HALF: SplitAngle = Bounds::atan_reciprocal(2).split();
pub(crate) const ATAN_ONE: SplitAngle = PI_OVER_FOUR.split();
// atan(3/2) = pi/4 + atan(1/5).
pub(crate) const ATAN_THREE_HALVES: SplitAngle = (Bounds {
    lo: 5 * ATAN_FIFTH.lo - ATAN_239TH.hi,
    hi: 5 * ATAN_FIFTH.hi - ATAN_239TH.lo,
})
.split();
pub(crate) const ATAN_INFINITY: SplitAngle = (Bounds {
    lo: 2 * PI_OVER_FOUR.lo,
    hi: 2 * PI_OVER_FOUR.hi,
})
.split();

const ATAN_FIFTH: Bounds = Bounds::atan_reciprocal(5);
const ATAN_239TH: Bounds = Bounds::atan_reciprocal(239);
// Machin's identity: pi/4 = 4*atan(1/5) - atan(1/239).
const PI_OVER_FOUR: Bounds = Bounds {
    lo: 4 * ATAN_FIFTH.lo - ATAN_239TH.hi,
    hi: 4 * ATAN_FIFTH.hi - ATAN_239TH.lo,
};

const SCALE: u128 = 1_u128 << (u128::BITS - 1);

/// Inclusive bounds on an angle, in units of 2^-127 radians.
#[derive(Clone, Copy)]
struct Bounds {
    lo: u128,
    hi: u128,
}

/// An angle represented as `hi + lo` to retain precision lost in `hi`.
pub(crate) struct SplitAngle {
    pub hi: f64,
    pub lo: f64,
}

impl Bounds {
    // atan(1/q) = 1/q - 1/(3q^3) + 1/(5q^5) - ...
    const fn atan_reciprocal(q: u128) -> Self {
        assert!(q >= 2 && q <= 239);
        let mut power = SCALE / q;
        let mut lo = 0;
        let mut hi = 0;
        let mut denominator = 1;
        let mut subtract = false;
        loop {
            let term = power / denominator;
            if term == 0 {
                return if subtract {
                    Bounds { lo: lo - 1, hi }
                } else {
                    Bounds { lo, hi: hi + 1 }
                };
            }
            // Round each term outward; the omitted tail is below one unit.
            if subtract {
                lo -= term + 1;
                hi -= term;
            } else {
                lo += term;
                hi += term + 1;
            }
            power /= q * q;
            denominator += 2;
            subtract = !subtract;
        }
    }

    const fn split(self) -> SplitAngle {
        let scale = SCALE as f64;
        let hi = self.lo as f64 / scale;
        // The interval must identify one rounded high part.
        assert!(hi.to_bits() == (self.hi as f64 / scale).to_bits());
        let hi_scaled = (hi * scale) as u128;
        // This split only supports a nonnegative remainder.
        assert!(hi_scaled <= self.lo);
        let lo = (self.lo - hi_scaled) as f64 / scale;
        // The remaining interval must also identify one rounded low part.
        assert!(lo.to_bits() == ((self.hi - hi_scaled) as f64 / scale).to_bits());
        SplitAngle { hi, lo }
    }
}

#[cfg(test)]
mod tests {
    use super::{Bounds, ATAN_HALF, ATAN_INFINITY, ATAN_ONE, ATAN_THREE_HALVES, SCALE};

    #[test]
    fn atan_reciprocal_two_returns_bounds_for_atan_half() {
        // atan(1/2) = 1/2 - 1/24 + 1/160 - ...
        // The first two partial sums bracket the answer:
        // 1/2 - 1/24 < atan(1/2) < 1/2.
        let bounds = Bounds::atan_reciprocal(2);

        // Bounds use units of 1/SCALE; the returned interval is tighter.
        assert!(bounds.lo >= SCALE / 2 - SCALE / 24);
        assert!(bounds.lo < bounds.hi);
        assert!(bounds.hi <= SCALE / 2);
    }

    #[test]
    fn reciprocal_bounds_enclose_high_precision_reference() {
        // floor(atan(1/q) * 2^127), checked with exact rational series bounds.
        for (q, reference) in [
            (2, 78_885_552_904_014_058_877_163_976_303_324_341_056_u128),
            (5, 33_585_114_162_700_596_429_338_312_019_514_155_576),
            (239, 711_883_642_681_542_234_893_201_432_822_774_392),
        ] {
            let bounds = Bounds::atan_reciprocal(q);
            // The true scaled value lies between reference and reference + 1.
            assert!(bounds.lo <= reference);
            assert!(bounds.hi > reference);
        }
    }

    #[test]
    fn split_preserves_sub_ulp_remainder() {
        // 0.5 + 2^-55 rounds to 0.5 as f64; lo retains the difference.
        const VALUE: u128 = SCALE / 2 + (1_u128 << 72);
        const SPLIT: super::SplitAngle = Bounds {
            lo: VALUE,
            hi: VALUE,
        }
        .split();
        assert_eq!(SPLIT.hi, 0.5);
        assert_eq!(SPLIT.lo, 1.0 / (1_u64 << 55) as f64);
    }

    #[test]
    #[should_panic]
    fn split_rejects_bounds_that_round_differently() {
        Bounds {
            lo: SCALE / 2,
            hi: SCALE,
        }
        .split();
    }

    #[test]
    fn reference_angles_preserve_original_bits() {
        for (value, hi, lo) in [
            (ATAN_HALF, 0x3fdd_ac67_0561_bb4f, 0x3c7a_2b7f_222f_65e2),
            (ATAN_ONE, 0x3fe9_21fb_5444_2d18, 0x3c81_a626_3314_5c07),
            (
                ATAN_THREE_HALVES,
                0x3fef_730b_d281_f69b,
                0x3c70_0788_7af0_cbbd,
            ),
            (ATAN_INFINITY, 0x3ff9_21fb_5444_2d18, 0x3c91_a626_3314_5c07),
        ] {
            assert_eq!(value.hi.to_bits(), hi);
            assert_eq!(value.lo.to_bits(), lo);
        }
    }
}
