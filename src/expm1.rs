/* origin: FreeBSD /usr/src/lib/msun/src/s_expm1.c */
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

const O_THRESHOLD: f64 = crate::ln(f64::MAX);
// Large negative inputs round to -1; this also triggers overflow and NaN checks.
const NEGATIVE_SATURATION_HIGH_WORD: u32 = crate::get_high_word(56.0 * core::f64::consts::LN_2);
// Below this cutoff, expm1(x) rounds to x.
const LINEAR_APPROXIMATION_HIGH_WORD: u32 = crate::get_high_word(crate::expi(2.0, -54));
const LN2_HI: f64 = crate::log_constants::LN_TWO_TRUNCATED.hi;
const LN2_LO: f64 = crate::log_constants::LN_TWO_TRUNCATED.lo;
const INVLN2: f64 = core::f64::consts::LOG2_E;
// Scaled polynomial coefficients for the reduced-argument approximation.
const Q1: f64 = -3.33333333333331316428e-02; /* BFA11111 111110F4 */
const Q2: f64 = 1.58730158725481460165e-03; /* 3F5A01A0 19FE5585 */
const Q3: f64 = -7.93650757867487942473e-05; /* BF14CE19 9EAADBB7 */
const Q4: f64 = 4.00821782732936239552e-06; /* 3ED0CFCA 86E65239 */
const Q5: f64 = -2.01099218183624371326e-07; /* BE8AFDB7 6E09C32D */

/// Computes exp(x) - 1 without cancellation near zero.
pub(crate) const fn expm1(x: f64) -> f64 {
    let high_word = crate::get_high_word(x.abs());
    let negative = x.is_sign_negative();

    if high_word >= NEGATIVE_SATURATION_HIGH_WORD {
        if x.is_nan() {
            return x;
        }
        if negative {
            return -1.0;
        }
        if x > O_THRESHOLD {
            return x * crate::constants::EXP_SCALE;
        }
    }

    // Reduce x to approximately [-ln(2)/2, ln(2)/2].
    let mut reduced = x;
    let mut correction = 0.0;
    let mut exponent = 0;
    if high_word > crate::constants::HALF_LN_2_HIGH_WORD {
        let hi;
        let lo;
        if high_word < crate::constants::THREE_HALVES_LN_2_HIGH_WORD {
            if negative {
                hi = x + LN2_HI;
                lo = -LN2_LO;
                exponent = -1;
            } else {
                hi = x - LN2_HI;
                lo = LN2_LO;
                exponent = 1;
            }
        } else {
            exponent = (INVLN2 * x + if negative { -0.5 } else { 0.5 }) as i32;
            let exponent_f64 = exponent as f64;
            // exponent * LN2_HI is exact here.
            hi = x - exponent_f64 * LN2_HI;
            lo = exponent_f64 * LN2_LO;
        }
        reduced = hi - lo;
        // Recover the rounding error in hi - lo.
        correction = (hi - reduced) - lo;
    } else if high_word < LINEAR_APPROXIMATION_HIGH_WORD {
        return x;
    }

    let half_x = 0.5 * reduced;
    let half_square = reduced * half_x;
    let polynomial = 1.0
        + half_square
            * (Q1
                + half_square * (Q2 + half_square * (Q3 + half_square * (Q4 + half_square * Q5))));
    let rational_term = 3.0 - polynomial * half_x;
    let rational_correction =
        half_square * ((polynomial - rational_term) / (6.0 - reduced * rational_term));

    // exp(reduced) - 1 ≈ reduced + half_square - reduced * rational_correction.
    if exponent == 0 {
        return reduced - (reduced * rational_correction - half_square);
    }
    let corrected_error = reduced * (rational_correction - correction) - correction;
    let corrected_error = corrected_error - half_square;

    reconstruct_expm1(reduced, corrected_error, exponent)
}

/// Reconstructs `exp(x) - 1` after reducing `x ≈ exponent * ln(2) + reduced`.
/// `exponent` is nonzero; `corrected_error` adjusts for higher-order terms and reduction error:
/// `exp(x) ≈ 2^exponent * (1 + reduced - corrected_error)`.
const fn reconstruct_expm1(reduced: f64, corrected_error: f64, exponent: i32) -> f64 {
    // Rearrange the result to retain precision for each exponent.
    match exponent {
        -1 => 0.5 * (reduced - corrected_error) - 0.5,
        0 => unreachable!(), // We check before we call this function that it's not 0
        1 => {
            if reduced < -0.25 {
                -2.0 * (corrected_error - (reduced + 0.5))
            } else {
                1.0 + 2.0 * (reduced - corrected_error)
            }
        }
        2..=56 => {
            let scale = normal_power_of_two(exponent);
            let inverse_scale = normal_power_of_two(-exponent);
            if exponent < 20 {
                (reduced - corrected_error + (1.0 - inverse_scale)) * scale
            } else {
                (reduced - (corrected_error + inverse_scale) + 1.0) * scale
            }
        }
        // 2^1024 is not representable, so scale in two steps.
        f64::MAX_EXP => (reduced - corrected_error + 1.0) * 2.0 * crate::constants::EXP_SCALE - 1.0,
        _ => (reduced - corrected_error + 1.0) * normal_power_of_two(exponent) - 1.0,
    }
}

// Construct a power of two in the normal f64 range.
const fn normal_power_of_two(exponent: i32) -> f64 {
    const EXPONENT_BIAS: i32 = f64::MAX_EXP - 1;
    const FRACTION_BITS: u32 = f64::MANTISSA_DIGITS - 1;
    f64::from_bits(((EXPONENT_BIAS + exponent) as u64) << FRACTION_BITS)
}

#[cfg(test)]
mod tests {
    use super::{expm1, normal_power_of_two, reconstruct_expm1};

    #[test]
    fn linear_approximation_cutoff() {
        let cutoff = f64::EPSILON / 4.0;
        let below = f64::from_bits(cutoff.to_bits() - 1);

        assert_eq!(expm1(below), below);
        assert_eq!(expm1(-below), -below);
        assert_eq!(expm1(cutoff), cutoff);
        assert_eq!(expm1(-cutoff), -cutoff);
    }

    #[test]
    fn argument_reduction_cutoffs() {
        // Adjacent inputs around the high-word cutoff near ln(2)/2.
        assert_eq!(expm1(0.34657359123229975), 0.41421356371988893);
        assert_eq!(expm1(0.3465735912322998), 0.41421356371988904);

        // Adjacent inputs around the high-word cutoff near 3*ln(2)/2.
        assert_eq!(expm1(-1.03972053527832), -0.6464465261231148);
        assert_eq!(expm1(-1.0397205352783203), -0.6464465261231149);
    }

    #[test]
    fn normal_power_of_two_values() {
        assert_eq!(normal_power_of_two(-2), 0.25);
        assert_eq!(normal_power_of_two(-1), 0.5);
        assert_eq!(normal_power_of_two(0), 1.0);
        assert_eq!(normal_power_of_two(1), 2.0);
        assert_eq!(normal_power_of_two(2), 4.0);
        assert_eq!(normal_power_of_two(-1022), f64::MIN_POSITIVE);
        assert_eq!(normal_power_of_two(1023), crate::expi(2.0, 1023));
    }

    #[test]
    fn reconstructs_powers_of_two_minus_one() {
        // A zero remainder and correction make exp(x) = 2^exponent.
        assert_eq!(reconstruct_expm1(0.0, 0.0, -2), -0.75);
        assert_eq!(reconstruct_expm1(0.0, 0.0, -1), -0.5);
        assert_eq!(reconstruct_expm1(0.0, 0.0, 1), 1.0);
        assert_eq!(reconstruct_expm1(0.0, 0.0, 2), 3.0);
        assert_eq!(reconstruct_expm1(0.0, 0.0, 19), 524_287.0);
        assert_eq!(reconstruct_expm1(0.0, 0.0, 20), 1_048_575.0);

        // At these magnitudes, subtracting 1 rounds away.
        assert_eq!(reconstruct_expm1(0.0, 0.0, 56), crate::expi(2.0, 56));
        assert_eq!(reconstruct_expm1(0.0, 0.0, 57), crate::expi(2.0, 57));
    }

    #[test]
    fn reconstructs_with_a_nonzero_remainder() {
        // For reduced = 0.25, corrected_error ≈ 1.25 - exp(0.25).
        assert_eq!(
            reconstruct_expm1(0.25, -0.034025416687741485, -1),
            -0.35798729165612925
        );
        assert_eq!(
            reconstruct_expm1(0.25, -0.034025416687741485, 1),
            1.568050833375483
        );
        assert_eq!(
            reconstruct_expm1(0.25, -0.034025416687741485, 2),
            4.136101666750966
        );
    }

    #[test]
    fn reconstructs_when_exponent_is_one_and_reduced_is_negative() {
        // corrected_error ≈ 1 - 0.3125 - exp(-0.3125).
        assert_eq!(
            reconstruct_expm1(-0.3125, -0.04411562894664179, 1),
            0.46323125789328357
        );
    }

    #[test]
    fn reconstructs_at_maximum_exponent() {
        // exp(-0.25) * 2^1024 is finite even though 2^1024 overflows.
        assert_eq!(
            reconstruct_expm1(-0.25, -0.028800783071404868, f64::MAX_EXP),
            1.4000448211528603e308
        );
        assert!(reconstruct_expm1(0.25, -0.034025416687741485, f64::MAX_EXP).is_infinite());
    }

    #[test]
    fn const_values() {
        // Rounded references computed with 100-digit decimal exponentials.
        assert_eq!((const { expm1(1e-20) }), 1e-20);
        assert_eq!((const { expm1(-1e-20) }), -1e-20);
        assert_eq!((const { expm1(0.5) }), 0.6487212707001282);
        assert_eq!((const { expm1(-0.5) }), -0.3934693402873666);
        assert_eq!((const { expm1(1.1) }), 2.0041660239464334);
        assert_eq!((const { expm1(-1.1) }), -0.6671289163019205);
        assert_eq!((const { expm1(20.0) }), 485165194.4097903);
        assert_eq!((const { expm1(100.0) }), 2.6881171418161356e43);
        assert_eq!((const { expm1(709.78) }), 1.7928227943945155e308);
    }

    #[test]
    fn special_values() {
        assert_eq!(expm1(-0.0).to_bits(), (-0.0_f64).to_bits());
        assert_eq!(expm1(f64::INFINITY), f64::INFINITY);
        assert_eq!(expm1(f64::NEG_INFINITY), -1.0);
        assert_eq!(expm1(710.0), f64::INFINITY);
        assert_eq!(expm1(-40.0), -1.0);
        assert!(expm1(f64::NAN).is_nan());
    }
}
