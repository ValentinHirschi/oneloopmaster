//! Real-mass four-mass D0 that depended on the order of its propagators.
//!
//! OneLOop's boxf4 chart (Denner-Nierste-Scharf), which four real masses take
//! outside the contour box's timelike configurations, was wrong in about 2% of
//! the relabellings of random four-mass boxes (OneLOop 3.6 itself in 118 of
//! the 6040 relabellings that took it in a sweep of 400 boxes; the expressions
//! in 156, in 54 boxes, as some lips were also decided differently at 128
//! bits). It is Denner and Dittmaier's eq. (3.61) (arXiv:1005.2076), which holds
//! only with the r-function in the p23 slot real, positive and exactly real;
//! and OneLOop's phase-tracked products `x r13 r_ij` differ from that formula's
//! principal `x r13` and explicit `eta(x, r13) Log R` term by multiples of
//! `(2 pi i)^2`. `box_integral::box_four_masses_chart` puts an invariant below
//! its pseudo-threshold in the p23 slot and evaluates eq. (3.61) with exact
//! infinitesimal decisions; boxes without such an invariant (all invariants
//! above their pseudo-thresholds, no real momenta) keep OneLOop's chart, which
//! needed the continued dilogarithm's unit-circle test on a real part in
//! `Float`. The rows of tests/data/four_mass_orderings.txt cover each mechanism.
#![cfg(feature = "generated-evaluators")]
use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend,
};
use symbolica::domains::float::{Complex, DoubleFloat, Float, RealLike};

type C = Complex<f64>;

/// Relative tolerance against the references (whose own error is below
/// 1e-15 here): row 9, four unit masses at invariants of order 1e5, keeps 10.7
/// digits in binary64. At 128 bits every route and image matches to 1e-15, and
/// all 24 images agree with each other to 1e-30.
const TOLERANCE: f64 = 5e-11;
const PRECISE: f64 = 1e-15;
const IMAGE_AGREEMENT: f64 = 1e-30;

fn rows() -> Vec<(Vec<C>, C)> {
    include_str!("data/four_mass_orderings.txt")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| {
            let v: Vec<f64> = line
                .split_whitespace()
                .map(|x| x.parse().unwrap())
                .collect();
            assert_eq!(v[0], 4.0);
            let mut args: Vec<C> = v[2..8].iter().map(|&x| C::new(x, 0.)).collect();
            args.extend(v[8..16].chunks_exact(2).map(|z| C::new(z[0], z[1])));
            args.push(C::new(v[1], 0.));
            (args, C::new(v[16], v[17]))
        })
        .collect()
}

/// The box under all 24 propagator relabellings.
fn images(args: &[C]) -> Vec<Vec<C>> {
    let invariant = |i: usize, j: usize| match (i.min(j), i.max(j)) {
        (0, 1) => args[0],
        (1, 2) => args[1],
        (2, 3) => args[2],
        (0, 3) => args[3],
        (0, 2) => args[4],
        (1, 3) => args[5],
        _ => unreachable!(),
    };
    let mut images = Vec::new();
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                let d = 6 - a - b - c;
                if a == b || b == c || a == c || d > 3 || d == a || d == b || d == c {
                    continue;
                }
                let mut image = vec![
                    invariant(a, b),
                    invariant(b, c),
                    invariant(c, d),
                    invariant(d, a),
                    invariant(a, c),
                    invariant(b, d),
                ];
                image.extend([a, b, c, d].map(|k| args[6 + k]));
                image.push(args[10]);
                images.push(image);
            }
        }
    }
    assert_eq!(images.len(), 24);
    images
}

fn relative(value: C, expected: C) -> f64 {
    (value.re - expected.re).hypot(value.im - expected.im) / expected.re.hypot(expected.im)
}

fn on_stack(action: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(256 * 1024 * 1024)
        .spawn(action)
        .unwrap()
        .join()
        .unwrap();
}

/// Binary64 routes (native with and without the hand-written kernels, the
/// interpreted expressions and SymJIT) and DoubleFloat, in every relabelling.
#[test]
fn four_mass_boxes_are_independent_of_the_labelling_in_binary64() {
    on_stack(|| {
        let rows = rows();
        assert_eq!(rows.len(), 17);
        let mut failures = Vec::new();
        let mut native: Vec<_> = [true, false]
            .map(|handwritten| {
                let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
                evaluator.set_handwritten(handwritten);
                (handwritten, evaluator)
            })
            .into();
        let mut double = NativeEvaluator::<DoubleFloat>::new(ScalarIntegral::D0).unwrap();
        for (line, (args, expected)) in rows.iter().enumerate() {
            for (number, image) in images(args).into_iter().enumerate() {
                let mut check = |route: String, value: C| {
                    let error = relative(value, *expected);
                    if !(error <= TOLERANCE) {
                        failures.push(format!(
                            "row {line} image {number} {route}: {value:?} vs {expected:?} ({error:e})"
                        ));
                    }
                };
                for (handwritten, evaluator) in &mut native {
                    let mut output = [C::new(f64::NAN, f64::NAN); 3];
                    evaluator.evaluate(&image, &mut output).unwrap();
                    check(format!("native (hand-written: {handwritten})"), output[0]);
                }
                for backend in [EvaluationBackend::Expression, EvaluationBackend::SymJit] {
                    let mut output = [C::new(f64::NAN, f64::NAN); 3];
                    evaluate_with_backend(ScalarIntegral::D0, &image, &mut output, backend)
                        .unwrap();
                    check(backend.name().to_string(), output[0]);
                }
                let input: Vec<_> = image
                    .iter()
                    .map(|z| Complex::new(DoubleFloat::from(z.re), DoubleFloat::from(z.im)))
                    .collect();
                let zero = Complex::new(DoubleFloat::from(0.), DoubleFloat::from(0.));
                let mut output = [zero; 3];
                double.evaluate(&input, &mut output).unwrap();
                check(
                    "DoubleFloat".into(),
                    C::new(output[0].re.to_f64(), output[0].im.to_f64()),
                );
            }
        }
        assert!(
            failures.is_empty(),
            "{} failures:\n{}",
            failures.len(),
            failures.join("\n")
        );
    });
}

/// Native code and interpreted expressions at 128 bits: every relabelling
/// matches the reference, and all 24 agree with each other to 1e-30.
#[test]
fn four_mass_boxes_are_independent_of_the_labelling_at_128_bits() {
    on_stack(|| {
        let bits = 128;
        let mut failures = Vec::new();
        let mut native =
            NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::D0, bits).unwrap();
        let mut expression = PrecisionEvaluator::with_binary_precision_and_backend(
            ScalarIntegral::D0,
            bits,
            EvaluationBackend::Expression,
        )
        .unwrap();
        for (line, (args, expected)) in rows().iter().enumerate() {
            for route in ["native", "expression"] {
                let mut values: Vec<Complex<Float>> = Vec::new();
                for image in images(args) {
                    let input: Vec<_> = image
                        .iter()
                        .map(|z| {
                            Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im))
                        })
                        .collect();
                    let mut output = core::array::from_fn::<_, 3, _>(|_| {
                        Complex::new(Float::new(bits), Float::new(bits))
                    });
                    if route == "native" {
                        native.evaluate(&input, &mut output).unwrap();
                    } else {
                        expression.evaluate(&input, &mut output).unwrap();
                    }
                    values.push(output[0].clone());
                }
                let first = &values[0];
                let scale = first.re.to_f64().hypot(first.im.to_f64());
                for (number, value) in values.iter().enumerate() {
                    let error = relative(C::new(value.re.to_f64(), value.im.to_f64()), *expected);
                    let spread = (value.re.clone() - &first.re)
                        .to_f64()
                        .hypot((value.im.clone() - &first.im).to_f64())
                        / scale;
                    if !(error <= PRECISE) || !(spread <= IMAGE_AGREEMENT) {
                        failures.push(format!(
                            "row {line} image {number} {route} 128 bits: {} {} vs {expected:?} \
                             (error {error:e}, spread from image 0 {spread:e})",
                            value.re, value.im
                        ));
                    }
                }
            }
        }
        assert!(
            failures.is_empty(),
            "{} failures:\n{}",
            failures.len(),
            failures.join("\n")
        );
    });
}
