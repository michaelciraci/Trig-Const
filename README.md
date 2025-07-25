## trig-const

Rust implementation of const trig functions.

This is implemented using a 16-term Taylor series approximation of cosine. Correctness is favored over speed, especially considering the main use case for this crate is to expose trigonometric functions for compile time.

The implemntation was largely inspired by the work of Dr. Austin Henley and Dr. Stephen Marz:
  - GitHub Repo: https://github.com/AZHenley/cosine
  - Article: https://austinhenley.com/blog/cosine.html

The implementation carries forward the original MIT license contained in the GitHub repo above.

## Requirements

This crate supports any compiler version back to rustc 1.85.1

```ignore
[dependencies]
trig-const = "0"
```

## Example

```
# use trig_const::cos;
# use std::f64::const::PI;
const COS_PI: f64 = cos(PI);
assert_eq!(COS_PI, -1.0);
```