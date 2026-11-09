//! Reciprocal factorials scaled by 2^128, evaluated exactly at compile time.

// Factorials through 34! fit in u128.
const fn factorial(n: u32) -> u128 {
    assert!(n <= 34);
    let mut factorial = 1_u128;
    let mut factor = 2;
    while factor <= n {
        factorial *= factor as u128;
        factor += 1;
    }
    factorial
}

/// floor(2^128 / n!).
const fn reciprocal_factorial_128(n: u32) -> u128 {
    assert!(n >= 2);
    // 2^128 itself does not fit in u128. Since n! is even for n >= 2,
    // divide numerator and denominator by two before integer division.
    (1_u128 << 127) / (factorial(n) / 2)
}

// The first omitted terms are r^34/34! and r^35/35!. For r < 0.79 < 4/5,
// r^34 < r^4 < 1/2, since 2*4^4 < 5^4. Together with 34! > 2^127,
// this bounds both omitted terms below 2^-128, using only u128 arithmetic.
const _: () = {
    assert!(2 * 4_u128.pow(4) < 5_u128.pow(4));
    assert!(factorial(34) > (1_u128 << 127));
};

// Coefficients of r^3, r^5, ..., r^33; signs come from Horner evaluation.
pub(super) const SIN_COEFFS: [u128; 16] = [
    reciprocal_factorial_128(3),
    reciprocal_factorial_128(5),
    reciprocal_factorial_128(7),
    reciprocal_factorial_128(9),
    reciprocal_factorial_128(11),
    reciprocal_factorial_128(13),
    reciprocal_factorial_128(15),
    reciprocal_factorial_128(17),
    reciprocal_factorial_128(19),
    reciprocal_factorial_128(21),
    reciprocal_factorial_128(23),
    reciprocal_factorial_128(25),
    reciprocal_factorial_128(27),
    reciprocal_factorial_128(29),
    reciprocal_factorial_128(31),
    reciprocal_factorial_128(33),
];

// Coefficients of r^2, r^4, ..., r^32; signs come from Horner evaluation.
pub(super) const COS_COEFFS: [u128; 16] = [
    reciprocal_factorial_128(2),
    reciprocal_factorial_128(4),
    reciprocal_factorial_128(6),
    reciprocal_factorial_128(8),
    reciprocal_factorial_128(10),
    reciprocal_factorial_128(12),
    reciprocal_factorial_128(14),
    reciprocal_factorial_128(16),
    reciprocal_factorial_128(18),
    reciprocal_factorial_128(20),
    reciprocal_factorial_128(22),
    reciprocal_factorial_128(24),
    reciprocal_factorial_128(26),
    reciprocal_factorial_128(28),
    reciprocal_factorial_128(30),
    reciprocal_factorial_128(32),
];
