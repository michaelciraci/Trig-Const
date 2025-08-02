use std::f64::consts::PI;

const STEP: f64 = 0.000001;

/// Run a sweep of precision comparisons against std
/// Precision will be different platform to platform
/// On my computer, I get:
/// ```ignore
/// Func   | Total Tests| Diff Count|       Max Diff
/// acos   |     2000000|     349419|    4.44089e-16
/// acosh  |    99000001|    8505318|    8.88178e-16
/// asin   |     2000000|     173790|    2.22045e-16
/// asinh  |    99000001|    8714913|    8.88178e-16
/// atan   |    50265483|    3290826|    2.22045e-16
/// atanh  |     1999998|     771064|    3.34115e-11
/// cos    |    50265483|    2173339|    1.11022e-16
/// cosh   |    25132742|    6498663|    2.91038e-11
/// ln     |    99999001|    3413955|    8.88178e-16
/// exp    |    20000001|    1944323|    3.63798e-12
/// fabs   |    20000001|          0|      0.00000e0
/// floor  |    20000001|          0|      0.00000e0
/// sin    |    50265483|    2255609|    1.11022e-16
/// sinh   |    25132742|    7200641|    2.91038e-11
/// sqrt   |    10000001|    2500953|    4.44089e-16
/// tan    |    50265483|   20777207|     3.72529e-9
/// ```
fn main() {
    println!(
        "{:<7}|{:>12}|{:>11}|{:>15}",
        "Func", "Total Tests", "Diff Count", "Max Diff"
    );
    let tests = [
        CompareArgs {
            name: "acos".to_string(),
            start: -1.0,
            stop: 1.0,
            step: STEP,
            std_fn: |x| x.acos(),
            const_fn: |x: f64| trig_const::acos(x),
        },
        CompareArgs {
            name: "acosh".to_string(),
            start: 1.0,
            stop: 100.0,
            step: STEP,
            std_fn: |x| x.acosh(),
            const_fn: |x: f64| trig_const::acosh(x),
        },
        CompareArgs {
            name: "asin".to_string(),
            start: -1.0,
            stop: 1.0,
            step: STEP,
            std_fn: |x| x.asin(),
            const_fn: |x: f64| trig_const::asin(x),
        },
        CompareArgs {
            name: "asinh".to_string(),
            start: 1.0,
            stop: 100.0,
            step: STEP,
            std_fn: |x| x.asinh(),
            const_fn: |x: f64| trig_const::asinh(x),
        },
        CompareArgs {
            name: "atan".to_string(),
            start: -8.0 * PI,
            stop: 8.0 * PI,
            step: STEP,
            std_fn: |x| x.atan(),
            const_fn: |x: f64| trig_const::atan(x),
        },
        CompareArgs {
            name: "atanh".to_string(),
            start: -1.0 + STEP,
            stop: 1.0 - STEP,
            step: STEP,
            std_fn: |x| x.atanh(),
            const_fn: |x| trig_const::atanh(x),
        },
        CompareArgs {
            name: "cos".to_string(),
            start: -8.0 * PI,
            stop: 8.0 * PI,
            step: STEP,
            std_fn: |x| x.cos(),
            const_fn: |x: f64| trig_const::cos(x),
        },
        CompareArgs {
            name: "cosh".to_string(),
            start: -4.0 * PI,
            stop: 4.0 * PI,
            step: STEP,
            std_fn: |x| x.cosh(),
            const_fn: |x| trig_const::cosh(x),
        },
        CompareArgs {
            name: "ln".to_string(),
            start: 0.001,
            stop: 100.0,
            step: STEP,
            std_fn: |x| x.ln(),
            const_fn: |x: f64| trig_const::ln(x),
        },
        CompareArgs {
            name: "exp".to_string(),
            start: -10.0,
            stop: 10.0,
            step: STEP,
            std_fn: |x| x.exp(),
            const_fn: |x: f64| trig_const::exp(x),
        },
        CompareArgs {
            name: "fabs".to_string(),
            start: -10.0,
            stop: 10.0,
            step: STEP,
            std_fn: |x| x.abs(),
            const_fn: |x: f64| trig_const::fabs(x),
        },
        CompareArgs {
            name: "floor".to_string(),
            start: -10.0,
            stop: 10.0,
            step: STEP,
            std_fn: |x| x.floor(),
            const_fn: |x: f64| trig_const::floor(x),
        },
        CompareArgs {
            name: "sin".to_string(),
            start: -8.0 * PI,
            stop: 8.0 * PI,
            step: STEP,
            std_fn: |x| x.sin(),
            const_fn: |x: f64| trig_const::sin(x),
        },
        CompareArgs {
            name: "sinh".to_string(),
            start: -4.0 * PI,
            stop: 4.0 * PI,
            step: STEP,
            std_fn: |x| x.sinh(),
            const_fn: |x| trig_const::sinh(x),
        },
        CompareArgs {
            name: "sqrt".to_string(),
            start: 0.0,
            stop: 10.0,
            step: STEP,
            std_fn: |x| x.sqrt(),
            const_fn: |x: f64| trig_const::sqrt(x),
        },
        CompareArgs {
            name: "tan".to_string(),
            start: -8.0 * PI,
            stop: 8.0 * PI,
            step: STEP,
            std_fn: |x| x.tan(),
            const_fn: |x: f64| trig_const::tan(x),
        },
    ];

    for test in tests {
        let diff = compare_functions(&test);
        println!(
            "{:<7}|{:>12}|{:>11}|{:>15.5e}",
            test.name, diff.total_tests, diff.diff_tests, diff.max_diff
        );
    }
}

fn compare_functions(c: &CompareArgs) -> DiffCounter {
    let mut const_metric = DiffCounter::default();

    for x in float_loop(c.start, c.stop, c.step) {
        let real = (c.std_fn)(x);
        let const_result = (c.const_fn)(x);

        const_metric.add_metric(real, const_result);
    }

    const_metric
}

#[derive(Debug, Default)]
struct DiffCounter {
    total_tests: usize,
    diff_tests: usize,
    max_diff: f64,
}

impl DiffCounter {
    fn add_metric(&mut self, real: f64, actual: f64) {
        self.total_tests += 1;
        let diff = (real - actual).abs();
        if diff != 0.0 {
            self.diff_tests += 1;
            self.max_diff = f64::max(self.max_diff, diff);
        }
    }
}

fn float_loop(start: f64, stop: f64, step: f64) -> impl Iterator<Item = f64> {
    core::iter::successors(Some(start), move |prev| {
        let next = prev + step;
        (next < stop).then_some(next)
    })
}

struct CompareArgs {
    name: String,
    start: f64,
    stop: f64,
    step: f64,
    std_fn: fn(f64) -> f64,
    const_fn: fn(f64) -> f64,
}
