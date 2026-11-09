//! Small, fixed-capacity unsigned integers shared by the trig implementations.
//!
//! Limbs are little-endian. There is no implicit scale or rounding mode here;
//! fixed-point interpretation belongs to `trig::Fixed`.

#[derive(Clone, Copy)]
pub(super) struct Uint<const N: usize> {
    pub(super) limbs: [u64; N],
}

impl<const N: usize> Uint<N> {
    pub(super) const ZERO: Self = Self { limbs: [0; N] };

    pub(super) const fn from_limbs(limbs: [u64; N]) -> Self {
        Self { limbs }
    }

    pub(super) const fn from_u64(value: u64) -> Self {
        let mut result = Self::ZERO;
        result.limbs[0] = value;
        result
    }

    pub(super) const fn power_of_two(bit: usize) -> Self {
        let mut result = Self::ZERO;
        result.limbs[bit / 64] = 1 << (bit % 64);
        result
    }

    pub(super) const fn less(self, rhs: Self) -> bool {
        let mut i = N;
        while i > 0 {
            i -= 1;
            if self.limbs[i] != rhs.limbs[i] {
                return self.limbs[i] < rhs.limbs[i];
            }
        }
        false
    }

    pub(super) const fn add(self, rhs: Self) -> Self {
        let mut result = Self::ZERO;
        let mut carry = 0_u128;
        let mut i = 0;
        while i < N {
            carry += self.limbs[i] as u128 + rhs.limbs[i] as u128;
            result.limbs[i] = carry as u64;
            carry >>= 64;
            i += 1;
        }
        assert!(carry == 0);
        result
    }

    pub(super) const fn sub(self, rhs: Self) -> Self {
        let mut result = Self::ZERO;
        let mut borrow = false;
        let mut i = 0;
        while i < N {
            let (value, a) = self.limbs[i].overflowing_sub(rhs.limbs[i]);
            let (value, b) = value.overflowing_sub(borrow as u64);
            result.limbs[i] = value;
            borrow = a || b;
            i += 1;
        }
        assert!(!borrow);
        result
    }

    pub(super) const fn mul_small(self, factor: u64) -> Self {
        let mut result = Self::ZERO;
        let mut carry = 0_u128;
        let mut i = 0;
        while i < N {
            carry += self.limbs[i] as u128 * factor as u128;
            result.limbs[i] = carry as u64;
            carry >>= 64;
            i += 1;
        }
        assert!(carry == 0);
        result
    }

    /// Integer quotient, rounded down. The divisor must be nonzero.
    pub(super) const fn div_small(self, divisor: u32) -> Self {
        // Small runtime values use 32-bit digits so division fits in u64.
        // The wider compile-time pi calculation uses 64-bit digits to halve
        // the interpreter's loop count. Both use the same long division.
        let digits_per_limb = if N <= 9 { 2 } else { 1 };
        let digit_bits = 64 / digits_per_limb;
        let mask = u64::MAX >> (64 - digit_bits);
        let mut result = Self::ZERO;
        let mut remainder = 0_u128;
        let mut i = digits_per_limb * N;
        while i > 0 {
            i -= 1;
            let shift = (i % digits_per_limb) * digit_bits;
            let digit = (self.limbs[i / digits_per_limb] >> shift) & mask;
            let numerator = (remainder << digit_bits) | digit as u128;
            let quotient = if digits_per_limb == 2 {
                (numerator as u64 / divisor as u64) as u128
            } else {
                numerator / divisor as u128
            };
            result.limbs[i / digits_per_limb] |= (quotient as u64) << shift;
            remainder = numerator - quotient * divisor as u128;
        }
        result
    }

    pub(super) const fn is_zero(self) -> bool {
        let mut i = 0;
        while i < N {
            if self.limbs[i] != 0 {
                return false;
            }
            i += 1;
        }
        true
    }

    pub(super) const fn bit(&self, bit: usize) -> bool {
        let word = bit / 64;
        word < N && (self.limbs[word] >> (bit % 64)) & 1 != 0
    }

    /// Read 64 bits at an arbitrary offset, zero-extending past the top limb.
    pub(super) const fn window(&self, bit: usize) -> u64 {
        let word = bit / 64;
        if word >= N {
            return 0;
        }
        let shift = bit % 64;
        let lo = self.limbs[word];
        if shift == 0 {
            lo
        } else {
            let hi = if word + 1 < N {
                self.limbs[word + 1]
            } else {
                0
            };
            (lo >> shift) | (hi << (64 - shift))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Uint;

    #[test]
    fn arithmetic_across_limbs() {
        let max = Uint::from_limbs([u64::MAX, u64::MAX, 0]);
        let one = Uint::from_u64(1);
        let next = Uint::from_limbs([0, 0, 1]);
        assert_eq!(max.add(one).limbs, next.limbs);
        assert_eq!(next.sub(one).limbs, max.limbs);
        assert_eq!(max.mul_small(2).limbs, [u64::MAX - 1, u64::MAX, 1]);
        assert!(max.less(next));
        assert!(!next.less(max));
        assert_eq!(max.window(64), u64::MAX);
        assert_eq!(max.window(65), u64::MAX >> 1);
        assert_eq!(max.window(192), 0);
        assert!(max.bit(127));
        assert!(!max.bit(128));
    }

    #[test]
    fn division_matches_native_integers() {
        // Exercise both the small and large integer division implementations.
        for value in [0, 1, u64::MAX as u128, 1_u128 << 64, u128::MAX] {
            let small = Uint::from_limbs([value as u64, (value >> 64) as u64, 0]);
            let mut large = Uint::<39>::ZERO;
            large.limbs[..3].copy_from_slice(&small.limbs);
            for divisor in [1, 3, 57121, u32::MAX] {
                let expected = value / divisor as u128;
                let expected_limbs = [expected as u64, (expected >> 64) as u64, 0];
                assert_eq!(small.div_small(divisor).limbs, expected_limbs);
                let quotient = large.div_small(divisor);
                assert_eq!(quotient.limbs[..3], expected_limbs);
                assert_eq!(quotient.limbs[3..], [0; 36]);
            }
        }
    }

    #[test]
    fn carries_and_division_across_all_limbs() {
        const MAX: Uint<39> = Uint::from_limbs([u64::MAX; 39]);
        assert_eq!(
            (const { MAX.div_small(3) }).limbs,
            [0x5555_5555_5555_5555; 39]
        );
        assert_eq!(
            (const { MAX.div_small(u32::MAX) }).limbs,
            [0x0000_0001_0000_0001; 39]
        );
        let power = Uint::<39>::power_of_two(64 * 38);
        let below = power.sub(Uint::from_u64(1));
        assert_eq!(below.limbs[..38], [u64::MAX; 38]);
        assert_eq!(below.add(Uint::from_u64(1)).limbs, power.limbs);
        assert!(power.sub(power).is_zero());
        assert_eq!(
            power.mul_small(2).limbs,
            Uint::power_of_two(64 * 38 + 1).limbs
        );
    }

    #[test]
    #[should_panic]
    fn addition_rejects_overflow() {
        Uint::from_limbs([u64::MAX; 3]).add(Uint::from_u64(1));
    }

    #[test]
    #[should_panic]
    fn subtraction_rejects_underflow() {
        Uint::<3>::ZERO.sub(Uint::from_u64(1));
    }
}
