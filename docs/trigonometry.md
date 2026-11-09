# How sine and cosine work

The implementation uses Taylor series with **integer intervals**: every
calculation keeps lower and upper bounds on the exact answer.

## 1. Reduce the angle

Rewrite the input as:

```text
x = n * pi/2 + r
```

The remainder `r` is near `[-pi/4, pi/4]`, where Taylor series converge quickly.
The quadrant (`n mod 4`) determines whether to evaluate sine or cosine and
which sign to use. For example, `sin(1)` becomes `cos(pi/2 - 1)`.

`INV_PIO2` (2/pi) finds the quadrant and remaining fraction; `PIO2` (pi/2)
converts that fraction back to radians. Both are computed at compile time
with integer bounds. Their extra precision preserves the remainder even for
very large inputs.

## 2. Evaluate the small angle

```text
sin(r) = r - r^3/3! + r^5/5! - ...
cos(r) = 1 - r^2/2! + r^4/4! - ...
```

The first pass uses 128 fractional bits. Coefficients such as `1/3!` are
computed at compile time, and nested multiplication (Horner's method) avoids
runtime division.

Intervals account for constant approximation, arithmetic rounding, and omitted
series terms. Difficult cases restart from the original input at 256 bits,
then 512 bits, computing Taylor terms directly.

## 3. Round the answer

Round both interval endpoints to the nearest `f64`, with ties to even.
If they agree, the answer is settled. Otherwise, retry at higher precision.

If the 512-bit attempt remains unresolved, round the interval's midpoint.
If angle reduction itself is unresolved, use the nearby quadrant boundary's
value (0 or +/-1). This avoids a panic, but correct rounding is not guaranteed
in these fallback cases. No tested input has needed this final fallback.

## Code and tests

- [Uint<N>](../src/uint.rs) provides unscaled integer arithmetic.
- [Fixed<N> and intervals](../src/trig.rs) supply scaling and rounding.
  A stored integer represents `integer / 2^(64*(N-1))`.
- [Pi constants](../src/trig_constants.rs) and
  [coefficients](../src/trig_coefficients.rs) use Rust `const fn` calculations.

Run `cargo test`. Inline references include MPFR-verified values and
[published hard-rounding cases](https://www.arith2026.org/papers/Computing%20hard-to-round%20cases%20of%20sin_2C%20cos_2C%20tan%20in%20double%20precisio.pdf).
Boundary tests also check agreement between precisions.

With `nightly`, public runtime calls still delegate to std/libm; const calls
use this implementation.
