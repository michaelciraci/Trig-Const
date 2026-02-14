## trig-const

[![Rust](https://github.com/michaelciraci/Trig-Const/actions/workflows/rust.yml/badge.svg?branch=master)](https://github.com/michaelciraci/Trig-Const/actions/workflows/rust.yml)
[![](https://img.shields.io/crates/v/trig-const.svg)](https://crates.io/crates/trig-const)
[![](https://img.shields.io/crates/l/trig-const.svg)](https://crates.io/crates/trig-const)
[![](https://docs.rs/trig-const/badge.svg)](https://docs.rs/trig-const/)
![minimum rustc 1.85](https://img.shields.io/badge/rustc-1.85+-red.svg)

Rust implementation of const trig functions.

The majority of functions have been implemented using a modified version of [libm](https://crates.io/crates/libm) for const Rust.

This implementation carries forward the original MIT license.

## Project Goals
 - Correctness while const (same result as std within a rounding error)
 - no-std
 - No unsafe

## Requirements

This crate supports any compiler version back to rustc 1.85.0

```toml
[dependencies]
trig-const = "0"
```

## Example

```rust
use std::f64::consts::PI;
use trig_const::cos;

const COS_PI: f64 = cos(PI);
assert_eq!(COS_PI, -1.0);
```

```rust
use std::f64::consts::PI;
use trig_const::{atan2, cos, sin};

/// 45° in radians
const DEG_45: f64 = 45.0 * PI / 180.0;

/// Pre-computed matrix to rotate object 45°
const ROTATIONAL_MATRIX: [[f64; 3]; 3] = [
    [cos(DEG_45), 0.0, sin(DEG_45)],
    [0.0, 1.0, 0.0],
    [-sin(DEG_45), 0.0, cos(DEG_45)],
];

/// atan2 calculation
const ATAN2_0_0: f64 = atan2(0.0, 0.0);

fn main() {
    println!("{:?}", ROTATIONAL_MATRIX);
    println!("{}", ATAN2_0_0);
}
```

## Features

- `nightly`: Running in [nightly](https://rust-lang.github.io/rustup/concepts/channels.html) exposes the function [`const_eval_select`](https://doc.rust-lang.org/std/intrinsics/fn.const_eval_select.html). This allows the compiler to call a different function if the function is to be evaluated at compile-time or run-time
  - Any `const` function calls will use this library
  - Any runtime function calls will forward to `libm`.
- `std`: This feature can be used in conjunction with `nightly`. This will forward runtime function calls to `std` instead of `libm`, to take advantage of hardware intrinsics (SIMD/FPU)

## Precision

Precision will be different platform to platform. There is a precision comparison within examples, under `examples/std_cmp.rs` (to run: `cargo run --release --example std_cmp`).

On aarch64, I get:

Func   | Total Tests| Diff Count|       Max Diff
-------|------------|-----------|---------------
acos   |     2000000|     349419|    4.44089e-16
acosh  |    99000001|    8505318|    8.88178e-16
asin   |     2000000|     173790|    2.22045e-16
asinh  |    99000001|    8714913|    8.88178e-16
atan   |    50265483|    3290826|    2.22045e-16
atanh  |     1999998|     771064|    3.34115e-11
cos    |    50265483|    2173339|    1.11022e-16
cosh   |    25132742|    6498663|    2.91038e-11
ln     |    99999001|    3413955|    8.88178e-16
exp    |    20000001|    1944323|    3.63798e-12
fabs   |    20000001|          0|      0.00000e0
floor  |    20000001|          0|      0.00000e0
sin    |    50265483|    2255609|    1.11022e-16
sinh   |    25132742|    7200641|    2.91038e-11
sqrt   |    10000001|    2500953|    4.44089e-16
tan    |    50265483|   20777207|     3.72529e-9

## History

This crate was originally implemented using trigonometric Taylor series approximations, inspired by the work of Dr. Austin Henley and Dr. Stephen Marz:
  - GitHub Repo: <https://github.com/AZHenley/cosine>
  - Article: <https://austinhenley.com/blog/cosine.html>

However, several functions have since been implemented using a modified version of [libm](https://crates.io/crates/libm) for const Rust which improved precision.
