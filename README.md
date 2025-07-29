## trig-const

Rust implementation of const trig functions.

This is implemented using Taylor series approximations. Correctness is prioritized over speed, especially considering the main use case for this crate is to expose trigonometric functions for compile-time evaluation.

The implemntation was largely inspired by the work of Dr. Austin Henley and Dr. Stephen Marz:
  - GitHub Repo: <https://github.com/AZHenley/cosine>
  - Article: <https://austinhenley.com/blog/cosine.html>

The implementation carries forward the original MIT license contained in the GitHub repo above.

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
use trig_const::{atan2, sin};

const SIN_PI_4: f64 = sin(PI / 2.0);
const ATAN2_0_0: f64 = atan2(0.0, 0.0);

fn main() {
    println!("{}\n{}", SIN_PI_4, ATAN2_0_0);
}
```
