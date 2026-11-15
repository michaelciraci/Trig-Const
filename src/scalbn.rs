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
    // Bits including the implicit bit
    let sig_total_bits = SIG_BITS + 1;

    // Maximum and minimum values when biased
    let exp_max = EXP_MAX;
    let exp_min = EXP_MIN;

    // 2 ^ Emax, maximum positive with null significand (0x1p1023 for f64)
    let f_exp_max = from_biased_exponent(EXP_BIAS << 1);

    // 2 ^ Emin, minimum positive normal with null significand (0x1p-1022 for f64)
    let f_exp_min = from_biased_exponent(1);

    // 2 ^ sig_total_bits, multiplier to normalize subnormals (0x1p53 for f64)
    let f_pow_subnorm = from_biased_exponent(sig_total_bits + EXP_BIAS);

    /*
     * The goal is to multiply `x` by a scale factor that applies `n`. However, there are cases
     * where `2^n` is not representable by `F` but the result should be, e.g. `x = 2^Emin` with
     * `n = -EMin + 2` (one out of range of 2^Emax). To get around this, reduce the magnitude of
     * the final scale operation by prescaling by the max/min power representable by `F`.
     */

    if n > exp_max {
        // Worst case positive `n`: scale the minimum subnormal to near f64::MAX.
        // This can be reached by three scaling multiplications (two here and one final).
        debug_assert!(-exp_min + SIG_BITS as i32 + exp_max <= exp_max * 3);

        x *= f_exp_max;
        n -= exp_max;
        if n > exp_max {
            x *= f_exp_max;
            n -= exp_max;
            if n > exp_max {
                n = exp_max;
            }
        }
    } else if n < exp_min {
        // Prescale to avoid rounding subnormal results twice.
        let mul = f_exp_min * f_pow_subnorm;
        let add = -exp_min - sig_total_bits as i32;

        // Two prescales and the final factor cover f64's range.
        debug_assert!(-exp_min + SIG_BITS as i32 + exp_max <= add * 2 + -exp_min);

        x *= mul;
        n += add;

        if n < exp_min {
            x *= mul;
            n += add;

            if n < exp_min {
                n = exp_min;
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
