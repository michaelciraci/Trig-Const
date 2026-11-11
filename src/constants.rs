//! Constants shared by related math functions.

use core::f64::consts::{FRAC_PI_4, LN_2};

pub(crate) const INFINITY_HIGH_WORD: u32 = crate::get_high_word(f64::INFINITY);
pub(crate) const PI_OVER_FOUR_HIGH_WORD: u32 = crate::get_high_word(FRAC_PI_4);
pub(crate) const HALF_LN_2_HIGH_WORD: u32 = crate::get_high_word(0.5 * LN_2);
pub(crate) const THREE_HALVES_LN_2_HIGH_WORD: u32 = crate::get_high_word(1.5 * LN_2);

// Scale used to split arguments into 24-bit chunks during pi/2 reduction.
pub(crate) const PIO2_CHUNK_SCALE: f64 = crate::expi(2.0, 24);
// Largest power of two representable as f64, used by exp and expm1.
pub(crate) const EXP_SCALE: f64 = crate::expi(2.0, (f64::MAX_EXP - 1) as isize);
