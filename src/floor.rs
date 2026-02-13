/* SPDX-License-Identifier: MIT
 * origin: musl src/math/floor.c */

//! Generic `floor` algorithm.
//!
//! Note that this uses the algorithm from musl's `floorf` rather than `floor` or `floorl` because
//! performance seems to be better (based on icount) and it does not seem to experience rounding
//! errors on i386.

const SIG_BITS: u32 = 52;
const BITS: u32 = 64;
const EXP_BITS: u32 = BITS - SIG_BITS - 1;
const EXP_SAT: u32 = (1 << EXP_BITS) - 1;
const EXP_BIAS: u32 = EXP_SAT >> 1;
const SIG_MASK: u64 = 4503599627370495;

pub const fn floor(x: f64) -> f64 {
    nightly_exp!(floor, floor_inner, x)
}

const fn floor_inner(x: f64) -> f64 {
    let zero = 0;

    let mut ix = x.to_bits();
    let e = exp_unbiased(x);

    // If the represented value has no fractional part, no truncation is needed.
    if e >= SIG_BITS as i32 {
        return x;
    }

    if e >= 0 {
        // |x| >= 1.0
        let m = SIG_MASK >> unsigned(e);
        if ix & m == zero {
            // Portion to be masked is already zero; no adjustment needed.
            return x;
        }

        if x.is_sign_negative() {
            ix += m;
        }

        ix &= !m;
        f64::from_bits(ix)
    } else if x.is_sign_positive() {
        // 0.0 <= x < 1.0; rounding down goes toward +0.0.
        0.0
    } else if ix << 1 != zero {
        // -1.0 < x < 0.0; rounding down goes toward -1.0.
        -1.0
    } else {
        // -0.0 remains unchanged
        x
    }
}

const fn ex(x: f64) -> u32 {
    (x.to_bits() >> SIG_BITS) as u32 & EXP_SAT
}

const fn signed(x: u32) -> i32 {
    x as i32
}

const fn unsigned(x: i32) -> u32 {
    x as u32
}

const fn exp_unbiased(x: f64) -> i32 {
    signed(ex(x)) - EXP_BIAS as i32
}
