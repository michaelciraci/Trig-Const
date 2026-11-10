//! Arctangent reference angles, computed with integer bounds at compile time.

pub(crate) const ATAN_HALF: SplitAngle = split(atan_reciprocal(2));
pub(crate) const ATAN_ONE: SplitAngle = split(PI_OVER_FOUR);
// atan(3/2) = pi/4 + atan(1/5).
pub(crate) const ATAN_THREE_HALVES: SplitAngle = split(Bounds {
    lo: 5 * ATAN_FIFTH.lo - ATAN_239TH.hi,
    hi: 5 * ATAN_FIFTH.hi - ATAN_239TH.lo,
});
pub(crate) const ATAN_INFINITY: SplitAngle = split(Bounds {
    lo: 2 * PI_OVER_FOUR.lo,
    hi: 2 * PI_OVER_FOUR.hi,
});

const ATAN_FIFTH: Bounds = atan_reciprocal(5);
const ATAN_239TH: Bounds = atan_reciprocal(239);
// Machin's identity: pi/4 = 4*atan(1/5) - atan(1/239).
const PI_OVER_FOUR: Bounds = Bounds {
    lo: 4 * ATAN_FIFTH.lo - ATAN_239TH.hi,
    hi: 4 * ATAN_FIFTH.hi - ATAN_239TH.lo,
};

const SCALE: u128 = 1_u128 << (u128::BITS - 1);

#[derive(Clone, Copy)]
struct Bounds {
    lo: u128,
    hi: u128,
}

pub(crate) struct SplitAngle {
    pub hi: f64,
    pub lo: f64,
}

// atan(1/q) = 1/q - 1/(3q^3) + 1/(5q^5) - ...
const fn atan_reciprocal(q: u128) -> Bounds {
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

const fn split(value: Bounds) -> SplitAngle {
    let scale = SCALE as f64;
    let hi = value.lo as f64 / scale;
    assert!(hi.to_bits() == (value.hi as f64 / scale).to_bits());
    let hi_scaled = (hi * scale) as u128;
    assert!(hi_scaled <= value.lo);
    let lo = (value.lo - hi_scaled) as f64 / scale;
    assert!(lo.to_bits() == ((value.hi - hi_scaled) as f64 / scale).to_bits());
    SplitAngle { hi, lo }
}

#[cfg(test)]
mod tests {
    use super::{ATAN_HALF, ATAN_INFINITY, ATAN_ONE, ATAN_THREE_HALVES};

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
