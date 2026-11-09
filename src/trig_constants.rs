//! Pi constants computed at compile time using exact integer bounds.
//!
//! Machin's identity: pi = 16*atan(1/5) - 4*atan(1/239).
//! We retain 2432 fractional bits while computing pi, then verify that both
//! bounds give the same requested 512-bit pi/2 and 2304-bit 2/pi constants.

use crate::uint::Uint;

// One integer limb and 38 fractional limbs, with 128 guard bits beyond 2/pi.
const LIMBS: usize = 39;
const FRACTION_BITS: usize = 64 * (LIMBS - 1);

type Integer = Uint<LIMBS>;
const SCALE: Integer = Integer::power_of_two(FRACTION_BITS);

// Bounds for atan(1/q) * 2^FRACTION_BITS. Integer division rounds each
// positive term down; adding one bounds it above. Successive divisions of
// the power remain exact floors: floor(floor(a/b)/c) = floor(a/(b*c)).
const fn atan_inverse(q: u32) -> (Integer, Integer) {
    let mut power = SCALE.div_small(q);
    let mut lo = Integer::ZERO;
    let mut hi = Integer::ZERO;
    let mut denominator = 1;
    let mut subtract = false;
    loop {
        let term = power.div_small(denominator);
        let upper = term.add(Integer::from_u64(1));
        if term.is_zero() {
            // The remaining alternating series has magnitude < one unit.
            return if subtract {
                (lo.sub(Integer::from_u64(1)), hi)
            } else {
                (lo, hi.add(Integer::from_u64(1)))
            };
        }
        if subtract {
            lo = lo.sub(upper);
            hi = hi.sub(term);
        } else {
            lo = lo.add(term);
            hi = hi.add(upper);
        }
        power = power.div_small(q * q);
        denominator += 2;
        subtract = !subtract;
    }
}

const PI_BOUNDS: (Integer, Integer) = {
    let (a_lo, a_hi) = atan_inverse(5);
    let (b_lo, b_hi) = atan_inverse(239);
    (
        a_lo.mul_small(16).sub(b_hi.mul_small(4)),
        a_hi.mul_small(16).sub(b_lo.mul_small(4)),
    )
};

const fn pi_over_two() -> [u64; 9] {
    let (lo, hi) = PI_BOUNDS;
    let mut result = [0; 9];
    let mut i = 0;
    while i < result.len() {
        let bit = FRACTION_BITS + 1 - 512 + i * 64;
        result[i] = lo.window(bit);
        assert!(result[i] == hi.window(bit));
        i += 1;
    }
    result
}

const fn two_over_pi() -> [u64; 36] {
    let (lo, hi) = PI_BOUNDS;
    // Binary long division of 2/pi. After each step, remainder < hi;
    // append the next quotient bit by doubling and subtracting if needed.
    let mut remainder = SCALE.mul_small(2);
    assert!(remainder.less(lo));
    let mut quotient = Integer::ZERO;
    let mut bit = 2304;
    while bit > 0 {
        bit -= 1;
        remainder = remainder.mul_small(2);
        if !remainder.less(hi) {
            remainder = remainder.sub(hi);
            quotient.limbs[bit / 64] |= 1 << (bit % 64);
        }
    }
    // Division by pi's upper bound gives a lower quotient Q. Check that the
    // lower bound yields the same quotient without doing a second division:
    // N - Q*lo = remainder + Q*(hi-lo) must still be less than lo.
    let width = hi.sub(lo);
    let mut i = 1;
    while i < LIMBS {
        assert!(width.limbs[i] == 0);
        i += 1;
    }
    assert!(remainder.add(quotient.mul_small(width.limbs[0])).less(lo));
    let mut result = [0; 36];
    i = 0;
    while i < result.len() {
        result[i] = quotient.limbs[i];
        i += 1;
    }
    result
}

// Little-endian limbs; each integer is the floor of the scaled constant.
pub(super) const PIO2: [u64; 9] = pi_over_two(); // floor((pi/2) * 2^512)
pub(super) const INV_PIO2: [u64; 36] = two_over_pi(); // floor((2/pi) * 2^2304)
