## trig-const

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

/// Pre-computed matrix to rotate object 30°
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

## History

This crate was originally implemented using trigonometric Taylor series approximations, inspired by the work of Dr. Austin Henley and Dr. Stephen Marz:
  - GitHub Repo: <https://github.com/AZHenley/cosine>
  - Article: <https://austinhenley.com/blog/cosine.html>

However, several functions have since been implemented using a modified version of [lib](https://crates.io/crates/libm) for const Rust which improved precision.
