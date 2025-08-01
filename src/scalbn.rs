const SIG_BITS: u32 = 52;
const BITS: u32 = 64;
const EXP_BITS: u32 = BITS - SIG_BITS - 1;
const EXP_SAT: u32 = (1 << EXP_BITS) - 1;
const EXP_BIAS: u32 = EXP_SAT >> 1;
const SIG_MASK: u64 = 4503599627370495;
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
    let zero = 0;

    // Bits including the implicit bit
    let sig_total_bits = SIG_BITS + 1;

    // Maximum and minimum values when biased
    let exp_max = EXP_MAX;
    let exp_min = EXP_MIN;

    // 2 ^ Emax, maximum positive with null significand (0x1p1023 for f64)
    let f_exp_max = from_parts(false, EXP_BIAS << 1, zero);

    // 2 ^ Emin, minimum positive normal with null significand (0x1p-1022 for f64)
    let f_exp_min = from_parts(false, 1, zero);

    // 2 ^ sig_total_bits, moltiplier to normalize subnormals (0x1p53 for f64)
    let f_pow_subnorm = from_parts(false, sig_total_bits + EXP_BIAS, zero);

    /*
     * The goal is to multiply `x` by a scale factor that applies `n`. However, there are cases
     * where `2^n` is not representable by `F` but the result should be, e.g. `x = 2^Emin` with
     * `n = -EMin + 2` (one out of range of 2^Emax). To get around this, reduce the magnitude of
     * the final scale operation by prescaling by the max/min power representable by `F`.
     */

    if n > exp_max {
        // Worse case positive `n`: `x`  is the minimum subnormal value, the result is `F::MAX`.
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
        // When scaling toward 0, the prescaling is limited to a value that does not allow `x` to
        // go subnormal. This avoids double rounding.
        if BITS > 16 {
            // `mul` s.t. `!(x * mul).is_subnormal() ∀ x`
            let mul = f_exp_min * f_pow_subnorm;
            let add = -exp_min - sig_total_bits as i32;

            // Worse case negative `n`: `x`  is the maximum positive value, the result is `F::MIN`.
            // This must be reachable by three scaling multiplications (two here and one final).
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
        } else {
            // `f16` is unique compared to other float types in that the difference between the
            // minimum exponent and the significand bits (`add = -exp_min - sig_total_bits`) is
            // small, only three. The above method depend on decrementing `n` by `add` two times;
            // for other float types this works out because `add` is a substantial fraction of
            // the exponent range. For `f16`, however, 3 is relatively small compared to the
            // exponent range (which is 39), so that requires ~10 prescale rounds rather than two.
            //
            // Work aroudn this by using a different algorithm that calculates the prescale
            // dynamically based on the maximum possible value. This adds more operations per round
            // since it needs to construct the scale, but works better in the general case.
            let add = -clamp(n + sig_total_bits as i32, exp_min, sig_total_bits as i32);
            let mul = from_parts(false, (EXP_BIAS as i32 - add) as u32, zero);

            x *= mul;
            n += add;

            if n < exp_min {
                let add = -clamp(n + sig_total_bits as i32, exp_min, sig_total_bits as i32);
                let mul = from_parts(false, (EXP_BIAS as i32 - add) as u32, zero);

                x *= mul;
                n += add;

                if n < exp_min {
                    n = exp_min;
                }
            }
        }
    }

    let scale = from_parts(false, (EXP_BIAS as i32 + n) as u32, zero);
    x * scale
}

const fn from_parts(negative: bool, exponent: u32, significand: i32) -> f64 {
    let sign = if negative { 1_u64 } else { 0 };
    f64::from_bits(
        (sign << (BITS - 1))
            | ((exponent as u64 & EXP_SAT as u64) << SIG_BITS)
            | (significand as u64 & SIG_MASK),
    )
}

const fn clamp(x: i32, min: i32, max: i32) -> i32 {
    assert!(min <= max);
    if x < min {
        min
    } else if x > max {
        max
    } else {
        x
    }
}
