use core::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_6, LN_2, PI};
use trig_const::{
    acos, acosh, asin, asinh, atan, atan2, atanh, cos, cosh, exp, ln, pow, sin, sinh, sqrt, tan,
};

const MIN_SUBNORMAL: f64 = f64::from_bits(1);
// Neighbors of 1 and the rounded pi/2.
const BELOW_ONE: f64 = f64::from_bits(1.0_f64.to_bits() - 1);
const ABOVE_ONE: f64 = f64::from_bits(1.0_f64.to_bits() + 1);
const BELOW_PI_2: f64 = f64::from_bits(FRAC_PI_2.to_bits() - 1);
const ABOVE_PI_2: f64 = f64::from_bits(FRAC_PI_2.to_bits() + 1);

const ONE_ULP: u64 = 1;
const EXACT: u64 = 0;

// Map f64 bits to numeric order so their difference counts ULPs.
fn ordered_bits(value: f64) -> u64 {
    let bits = value.to_bits();
    if value.is_sign_negative() {
        !bits
    } else {
        bits | (1 << 63)
    }
}

fn assert_within_ulps(label: &str, actual: f64, expected: f64, max_ulps: u64) {
    assert!(
        actual.is_finite() && expected.is_finite(),
        "{label}: expected finite {expected:?}, got {actual:?}"
    );
    let ulps = ordered_bits(actual).abs_diff(ordered_bits(expected));
    assert!(
        ulps <= max_ulps,
        "{label}: got {actual:?}, expected {expected:?} (difference: {ulps} ULP; limit: {max_ulps})"
    );
}

// (Description, actual, reference.)
type AccuracyCase = (&'static str, f64, f64);

// Const cases bypass nightly runtime forwarding.
fn assert_cases(cases: &[AccuracyCase], max_ulps: u64) {
    for &(label, actual, expected) in cases {
        assert_within_ulps(label, actual, expected, max_ulps);
    }
}

#[test]
fn adjacent_values_are_one_ulp_apart() {
    let one = 1.0_f64;
    let next = 1.0 + f64::EPSILON;
    assert_eq!(ordered_bits(one).abs_diff(ordered_bits(next)), ONE_ULP);
    assert_eq!(ordered_bits(-one).abs_diff(ordered_bits(-next)), ONE_ULP);
}

#[test]
fn sine_accuracy() {
    // Large inputs exercise argument reduction.
    const CASES: &[AccuracyCase] = &[
        ("min subnormal", sin(MIN_SUBNORMAL), MIN_SUBNORMAL),
        ("-0.1", sin(-0.1), -0.09983341664682815),
        ("below pi/2", sin(BELOW_PI_2), 1.0),
        ("1e20", sin(1e20), -0.6452512852657808),
        ("max finite", sin(f64::MAX), 0.004961954789184062),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn cosine_accuracy() {
    const CASES: &[AccuracyCase] = &[
        ("min subnormal", cos(MIN_SUBNORMAL), 1.0),
        ("below pi/2", cos(BELOW_PI_2), 2.83276944882399e-16),
        ("above pi/2", cos(ABOVE_PI_2), -1.6081226496766366e-16),
        ("1e20", cos(1e20), 0.7639704044417283),
        ("max finite", cos(f64::MAX), -0.9999876894265599),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn tangent_accuracy() {
    // Probe both sides of the pi/2 pole.
    const CASES: &[AccuracyCase] = &[
        ("min subnormal", tan(MIN_SUBNORMAL), MIN_SUBNORMAL),
        ("below pi/2", tan(BELOW_PI_2), 3_530_114_321_217_157.5),
        ("above pi/2", tan(ABOVE_PI_2), -6_218_431_163_823_738.0),
        (
            "-10.9955742230756",
            tan(-10.9955742230756),
            -15_506_598.359_940_352,
        ),
        ("max finite", tan(f64::MAX), -0.004962015874444895),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn arcsine_accuracy() {
    const CASES: &[AccuracyCase] = &[
        ("min subnormal", asin(MIN_SUBNORMAL), MIN_SUBNORMAL),
        ("-0.5", asin(-0.5), -FRAC_PI_6),
        ("just below 1", asin(BELOW_ONE), 1.5707963118937354),
        ("just above -1", asin(-BELOW_ONE), -1.5707963118937354),
        ("0.975", asin(0.975), 1.3467210414930773),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn arccosine_accuracy() {
    const CASES: &[AccuracyCase] = &[
        ("just above -1", acos(-BELOW_ONE), 3.141592638688632),
        ("-0.5", acos(-0.5), 2.0943951023931957),
        ("0.5", acos(0.5), FRAC_PI_3),
        ("just below 1", acos(BELOW_ONE), 1.4901161193847656e-8),
        ("0.975", acos(0.975), 0.22407528530181936),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn arctangent_accuracy() {
    const CASES: &[AccuracyCase] = &[
        ("min subnormal", atan(MIN_SUBNORMAL), MIN_SUBNORMAL),
        ("0.4375", atan(0.4375), 0.4124104415973873),
        ("0.5", atan(0.5), 0.4636476090008061),
        ("1.5", atan(1.5), 0.982793723247329),
        ("max finite", atan(f64::MAX), FRAC_PI_2),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn atan2_accuracy() {
    const CASES: &[AccuracyCase] = &[
        ("1, -1", atan2(1.0, -1.0), 2.356194490192345),
        ("-1, -1", atan2(-1.0, -1.0), -2.356194490192345),
        (
            "min subnormal, min normal",
            atan2(MIN_SUBNORMAL, f64::MIN_POSITIVE),
            f64::EPSILON,
        ),
        ("max, 1", atan2(f64::MAX, 1.0), FRAC_PI_2),
        ("1, -max", atan2(1.0, -f64::MAX), PI),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn inverse_hyperbolic_sine_accuracy() {
    const CASES: &[AccuracyCase] = &[
        ("min subnormal", asinh(MIN_SUBNORMAL), MIN_SUBNORMAL),
        ("-0.5", asinh(-0.5), -0.48121182505960347),
        ("2.0", asinh(2.0), 1.4436354751788103),
        ("max finite", asinh(f64::MAX), 710.475860073944),
        ("negative max finite", asinh(-f64::MAX), -710.475860073944),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn inverse_hyperbolic_cosine_accuracy() {
    const CASES: &[AccuracyCase] = &[
        ("just above 1", acosh(ABOVE_ONE), 2.1073424255447014e-8),
        ("1.5", acosh(1.5), 0.9624236501192069),
        ("2.0", acosh(2.0), 1.3169578969248168),
        ("1e20", acosh(1e20), 46.74484904044086),
        ("max finite", acosh(f64::MAX), 710.475860073944),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn inverse_hyperbolic_tangent_accuracy() {
    const CASES: &[AccuracyCase] = &[
        ("min subnormal", atanh(MIN_SUBNORMAL), MIN_SUBNORMAL),
        ("-0.5", atanh(-0.5), -0.5493061443340549),
        ("0.5", atanh(0.5), 0.5493061443340549),
        ("just below 1", atanh(BELOW_ONE), 18.714973875118524),
        ("just above -1", atanh(-BELOW_ONE), -18.714973875118524),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn logarithm_accuracy() {
    const CASES: &[AccuracyCase] = &[
        ("min subnormal", ln(MIN_SUBNORMAL), -744.4400719213812),
        ("min normal", ln(f64::MIN_POSITIVE), -708.3964185322641),
        ("0.5", ln(0.5), -LN_2),
        ("just above 1", ln(ABOVE_ONE), 2.2204460492503128e-16),
        ("max finite", ln(f64::MAX), 709.782712893384),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn exponential_accuracy() {
    const CASES: &[AccuracyCase] = &[
        ("underflow boundary", exp(-745.1332191019411), MIN_SUBNORMAL),
        ("subnormal result", exp(-744.0), f64::from_bits(2)),
        (
            "normal/subnormal boundary",
            exp(-708.3964185322641),
            2.2250738585072626e-308,
        ),
        ("0.5", exp(0.5), 1.6487212707001282),
        ("709.78", exp(709.78), 1.7928227943945155e308),
        (
            "last finite result",
            exp(709.782712893384),
            1.7976931348622732e308,
        ),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn hyperbolic_sine_accuracy() {
    const CASES: &[AccuracyCase] = &[
        ("min subnormal", sinh(MIN_SUBNORMAL), MIN_SUBNORMAL),
        ("-0.5", sinh(-0.5), -0.5210953054937474),
        ("1.0", sinh(1.0), 1.1752011936438014),
        ("709.78", sinh(709.78), 8.964113971972578e307),
        ("710.47", sinh(710.47), 1.7871893267684048e308),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn hyperbolic_cosine_accuracy() {
    const CASES: &[AccuracyCase] = &[
        ("min subnormal", cosh(MIN_SUBNORMAL), 1.0),
        ("0.5", cosh(0.5), 1.1276259652063807),
        ("1.0", cosh(1.0), 1.5430806348152437),
        ("709.78", cosh(709.78), 8.964113971972578e307),
        ("710.47", cosh(710.47), 1.7871893267684048e308),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn square_root_accuracy() {
    const CASES: &[AccuracyCase] = &[
        (
            "min subnormal",
            sqrt(MIN_SUBNORMAL),
            2.2227587494850775e-162,
        ),
        (
            "min normal",
            sqrt(f64::MIN_POSITIVE),
            1.4916681462400413e-154,
        ),
        ("2.0", sqrt(2.0), core::f64::consts::SQRT_2),
        ("max finite", sqrt(f64::MAX), 1.3407807929942596e154),
    ];
    assert_cases(CASES, EXACT);
}

#[test]
fn power_accuracy() {
    const CASES: &[AccuracyCase] = &[
        ("2^-1074", pow(2.0, -1074.0), MIN_SUBNORMAL),
        ("(1e155)^-2", pow(1e155, -2.0), 1e-310),
        (
            "(1 + ulp)^2^31",
            pow(ABOVE_ONE, 2_147_483_648.0),
            1.000000476837272,
        ),
        ("0.5^1022", pow(0.5, 1022.0), f64::MIN_POSITIVE),
        ("(-1.5)^3", pow(-1.5, 3.0), -3.375),
        (
            "1.000001^100000",
            pow(1.000001, 100_000.0),
            1.105170862808048,
        ),
    ];
    assert_cases(CASES, ONE_ULP);
}

#[test]
fn underflow_boundaries() {
    // One ULP below the smallest subnormal is zero.
    for (label, actual) in [
        ("sin", const { sin(MIN_SUBNORMAL) }),
        ("tan", const { tan(MIN_SUBNORMAL) }),
        ("asin", const { asin(MIN_SUBNORMAL) }),
        ("atan", const { atan(MIN_SUBNORMAL) }),
        ("asinh", const { asinh(MIN_SUBNORMAL) }),
        ("atanh", const { atanh(MIN_SUBNORMAL) }),
        ("sinh", const { sinh(MIN_SUBNORMAL) }),
        ("exp", const { exp(-745.1332191019411) }),
        ("pow", const { pow(2.0, -1074.0) }),
    ] {
        assert_eq!(actual.to_bits(), MIN_SUBNORMAL.to_bits(), "{label}");
    }
    assert_eq!(
        (const { exp(-745.1332191019412) }).to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!((const { pow(2.0, -1075.0) }).to_bits(), 0.0_f64.to_bits());
}

#[test]
fn exponential_overflow_boundary() {
    assert_eq!((const { exp(709.7827128933841) }), f64::INFINITY);
}

#[test]
fn trig_special_values() {
    // Float equality doesn't distinguish signed zero.
    assert_eq!((const { sin(0.0) }).to_bits(), 0.0_f64.to_bits());
    assert_eq!((const { sin(-0.0) }).to_bits(), (-0.0_f64).to_bits());
    assert_eq!((const { tan(0.0) }).to_bits(), 0.0_f64.to_bits());
    assert_eq!((const { tan(-0.0) }).to_bits(), (-0.0_f64).to_bits());
    assert_eq!((const { cos(-0.0) }), 1.0);
    assert!(const { sin(f64::INFINITY) }.is_nan());
    assert!(const { cos(f64::NEG_INFINITY) }.is_nan());
    assert!(const { tan(f64::NAN) }.is_nan());
}

#[test]
fn atan2_special_values() {
    assert_eq!((const { atan2(0.0, 0.0) }).to_bits(), 0.0_f64.to_bits());
    assert_eq!((const { atan2(-0.0, 0.0) }).to_bits(), (-0.0_f64).to_bits());
    assert_eq!((const { atan2(0.0, -0.0) }), PI);
    assert_eq!((const { atan2(-0.0, -0.0) }), -PI);
    assert_eq!((const { atan2(f64::INFINITY, 1.0) }), FRAC_PI_2);
    assert_eq!((const { atan2(f64::NEG_INFINITY, 1.0) }), -FRAC_PI_2);
    assert!(const { atan2(f64::NAN, 1.0) }.is_nan());
}
