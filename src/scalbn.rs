const SIG_BITS: u32 = f64::MANTISSA_DIGITS - 1;
const BITS: u32 = u64::BITS;
const EXP_BITS: u32 = BITS - SIG_BITS - 1;
const EXP_SAT: u32 = (1 << EXP_BITS) - 1;
const EXP_BIAS: u32 = EXP_SAT >> 1;
const EXP_MAX: i32 = EXP_BIAS as i32;
const EXP_MIN: i32 = -(EXP_MAX - 1);

/// Scale the exponent.
///
/// From N3220:
///
/// > The scalbn and scalbln functions compute `x * b^n`, where `b = FLT_RADIX` if the return type
/// > of the function is a standard floating type, or `b = 10` if the return type of the function
/// > is a decimal floating type. A range error occurs for some finite x, depending on n.
/// >
/// > [...]
/// >
/// > * `scalbn(±0, n)` returns `±0`.
/// > * `scalbn(x, 0)` returns `x`.
/// > * `scalbn(±∞, n)` returns `±∞`.
/// >
/// > If the calculation does not overflow or underflow, the returned value is exact and
/// > independent of the current rounding direction mode.
pub(crate) const fn scalbn(mut x: f64, mut n: i32) -> f64 {
    const SIG_TOTAL_BITS: u32 = SIG_BITS + 1;
    const MAX_SCALE: f64 = from_biased_exponent(EXP_BIAS << 1);
    const MIN_SCALE: f64 = f64::MIN_POSITIVE;
    const SUBNORMAL_SCALE: f64 = from_biased_exponent(SIG_TOTAL_BITS + EXP_BIAS);
    const PRESCALE: f64 = MIN_SCALE * SUBNORMAL_SCALE;
    const PRESCALE_SHIFT: i32 = -EXP_MIN - SIG_TOTAL_BITS as i32;

    // Scale in stages when 2^n cannot be represented as an f64.
    if n > EXP_MAX {
        // Worst case positive `n`: scale the minimum subnormal to near f64::MAX.
        // This can be reached by three scaling multiplications (two here and one final).
        debug_assert!(-EXP_MIN + SIG_BITS as i32 + EXP_MAX <= EXP_MAX * 3);

        x *= MAX_SCALE;
        n -= EXP_MAX;
        if n > EXP_MAX {
            x *= MAX_SCALE;
            n -= EXP_MAX;
            if n > EXP_MAX {
                n = EXP_MAX;
            }
        }
    } else if n < EXP_MIN {
        // Prescale to avoid rounding subnormal results twice.
        // Two prescales and the final factor cover f64's range.
        debug_assert!(-EXP_MIN + SIG_BITS as i32 + EXP_MAX <= PRESCALE_SHIFT * 2 + -EXP_MIN);

        x *= PRESCALE;
        n += PRESCALE_SHIFT;

        if n < EXP_MIN {
            x *= PRESCALE;
            n += PRESCALE_SHIFT;

            if n < EXP_MIN {
                n = EXP_MIN;
            }
        }
    }

    let scale = from_biased_exponent((EXP_BIAS as i32 + n) as u32);
    x * scale
}

const fn from_biased_exponent(exponent: u32) -> f64 {
    f64::from_bits((exponent as u64 & EXP_SAT as u64) << SIG_BITS)
}

#[cfg(test)]
mod tests {
    use super::scalbn;

    const MIN_SUBNORMAL: f64 = f64::from_bits(1);

    #[test]
    fn scales_across_normal_and_subnormal_range() {
        assert_eq!(
            (const { scalbn(f64::MIN_POSITIVE, -1) }),
            f64::MIN_POSITIVE / 2.0
        );
        assert_eq!((const { scalbn(MIN_SUBNORMAL, 52) }), f64::MIN_POSITIVE);

        // Halfway results round to the even subnormal.
        assert_eq!((const { scalbn(f64::from_bits(3), -1) }).to_bits(), 2);
        assert_eq!((const { scalbn(f64::from_bits(5), -1) }).to_bits(), 2);
        assert_eq!(
            (const { scalbn(-MIN_SUBNORMAL, -1) }).to_bits(),
            (-0.0_f64).to_bits()
        );
    }

    #[test]
    fn scales_extreme_exponents() {
        const TWO_POW_1023: f64 = f64::from_bits(0x7fe0_0000_0000_0000);
        assert_eq!((const { scalbn(MIN_SUBNORMAL, 2097) }), TWO_POW_1023);
        assert_eq!((const { scalbn(f64::MAX, -2097) }).to_bits(), 2);
        assert_eq!((const { scalbn(f64::MAX, -2098) }), MIN_SUBNORMAL);
        assert_eq!((const { scalbn(f64::MAX, -2099) }).to_bits(), 0);
        assert_eq!((const { scalbn(MIN_SUBNORMAL, i32::MAX) }), f64::INFINITY);
        assert_eq!((const { scalbn(f64::MAX, i32::MIN) }).to_bits(), 0);
    }

    #[test]
    fn preserves_special_values() {
        for exponent in [i32::MIN, -1, 0, 1, i32::MAX] {
            assert_eq!(scalbn(0.0, exponent).to_bits(), 0.0_f64.to_bits());
            assert_eq!(scalbn(-0.0, exponent).to_bits(), (-0.0_f64).to_bits());
            assert_eq!(scalbn(f64::INFINITY, exponent), f64::INFINITY);
            assert_eq!(scalbn(f64::NEG_INFINITY, exponent), f64::NEG_INFINITY);
            assert!(scalbn(f64::NAN, exponent).is_nan());
        }
    }
}
