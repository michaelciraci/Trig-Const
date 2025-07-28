use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;
use trig_const::{acos, asin, asinh, atan, cos, ln, sin, tan};

/// Benchmarks for the core trigonometric functions (sin, cos, tan).
fn bench_core_trig(c: &mut Criterion) {
    let mut group = c.benchmark_group("Core Trig");

    // Test a simple value.
    group.bench_function("cos(1.5)", |b| b.iter(|| cos(black_box(1.5))));
    group.bench_function("sin(1.5)", |b| b.iter(|| sin(black_box(1.5))));

    // Test a value that requires range reduction.
    group.bench_function("cos(10.0)", |b| b.iter(|| cos(black_box(10.0))));

    // Test tan, which involves two function calls and a division.
    group.bench_function("tan(1.0)", |b| b.iter(|| tan(black_box(1.0))));

    group.finish();
}

/// Benchmarks for the inverse trigonometric functions.
fn bench_inverse_trig(c: &mut Criterion) {
    let mut group = c.benchmark_group("Inverse Trig");

    // Test asin in its fast region (no range reduction needed).
    group.bench_function("asin(0.4)", |b| b.iter(|| asin(black_box(0.4))));
    // Test asin in its slow region (requires range reduction).
    group.bench_function("asin(0.9)", |b| b.iter(|| asin(black_box(0.9))));

    // Test acos, which derives from asin.
    group.bench_function("acos(0.5)", |b| b.iter(|| acos(black_box(0.5))));

    // Test atan in its fast region.
    group.bench_function("atan(0.5)", |b| b.iter(|| atan(black_box(0.5))));
    // Test atan at its slowest-converging input.
    group.bench_function("atan(0.99)", |b| b.iter(|| atan(black_box(0.99))));

    group.finish();
}

/// Benchmarks for logarithmic and hyperbolic functions.
fn bench_log_hyperbolic(c: &mut Criterion) {
    let mut group = c.benchmark_group("Log & Hyperbolic");

    // Test ln where its series converges relatively quickly.
    group.bench_function("ln(1.1)", |b| b.iter(|| ln(black_box(1.1))));
    // Test ln at its slowest-converging input.
    group.bench_function("ln(1.99)", |b| b.iter(|| ln(black_box(1.99))));

    // Test asinh, which depends on ln and sqrt.
    group.bench_function("asinh(2.0)", |b| b.iter(|| asinh(black_box(2.0))));

    group.finish();
}

// Register all benchmark groups with criterion's main harness.
criterion_group!(
    benches,
    bench_core_trig,
    bench_inverse_trig,
    bench_log_hyperbolic
);
criterion_main!(benches);
