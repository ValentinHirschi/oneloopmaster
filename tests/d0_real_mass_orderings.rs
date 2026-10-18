//! Real-mass D0 that depended on the order of its propagators (NinjaRS):
//! p² = (20602, 6775, 113296, -106553), s12 = 80235, s23 = -4769,
//! m² = (8315.178393760001, 100, 8315.178393760001, 0), mu² = 8315.178393760001.
//! In the 8 relabellings with these channels the finite three-mass box took
//! OneLOop's real-mass boxf33 chart and had the real part +1.714058e-8 instead
//! of -9.718748e-9, on every route and at every precision except DoubleFloat in
//! the 4 cyclic images (OneLOop 3.6 gives the same wrong value in double and
//! quadruple precision; COLLIER -9.718748e-9 - 8.341413e-9 i). OneLOop takes its
//! contour representation (boxc) for real masses only in timelike
//! configurations; elsewhere boxf33 depends on the labelling. In a sweep of
//! 400 random real-mass boxes, 39 of 154 three-mass boxes were wrong in some of
//! their 24 relabellings, always where the contour box supports the
//! kinematics, which three-mass boxes now use wherever it does. (Four-mass
//! boxes keep OneLOop's selection: its boxf4 chart was wrong in 12 of the
//! 1432 relabellings that took it, in 8 of 88 boxes, but the contour box costs
//! 1.5 to 4 times as much on the recorded NinjaRS streams.)
//! References (tests/data/real_mass_orderings.txt): OneLOop 3.6 in quadruple
//! precision with vanishing widths, which is independent of the real-mass
//! charts.
#![cfg(feature = "generated-evaluators")]
use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend,
};
use symbolica::domains::float::{Complex, DoubleFloat, Float, RealLike};

type C = Complex<f64>;

/// Relative tolerance against the extrapolated references (whose own error
/// is about 1e-13), and between the 128-bit images of one box.
const TOLERANCE: f64 = 1e-11;
const IMAGE_AGREEMENT: f64 = 1e-30;

fn rows() -> Vec<(Vec<C>, C)> {
    include_str!("data/real_mass_orderings.txt")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| {
            let v: Vec<f64> = line.split_whitespace().map(|x| x.parse().unwrap()).collect();
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
fn real_mass_boxes_are_independent_of_the_labelling_in_binary64() {
    on_stack(|| {
        let rows = rows();
        assert!(rows.len() >= 13);
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
        assert!(failures.is_empty(), "{} failures:\n{}", failures.len(), failures.join("\n"));
    });
}

/// Native code and interpreted expressions at 128 bits: every relabelling
/// matches the reference, and all 24 agree with each other to 1e-30.
#[test]
fn real_mass_boxes_are_independent_of_the_labelling_at_128_bits() {
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
                        .map(|z| Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im)))
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
                    let error =
                        relative(C::new(value.re.to_f64(), value.im.to_f64()), *expected);
                    let spread = (value.re.clone() - &first.re)
                        .to_f64()
                        .hypot((value.im.clone() - &first.im).to_f64())
                        / scale;
                    if !(error <= TOLERANCE) || !(spread <= IMAGE_AGREEMENT) {
                        failures.push(format!(
                            "row {line} image {number} {route} 128 bits: {} {} vs {expected:?} \
                             (error {error:e}, spread from image 0 {spread:e})",
                            value.re, value.im
                        ));
                    }
                }
            }
        }
        assert!(failures.is_empty(), "{} failures:\n{}", failures.len(), failures.join("\n"));
    });
}
