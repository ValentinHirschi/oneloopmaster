//! One-mass D0 with a complex mass (OneLOop's boxf1), in every relabelling.
//! OneLOop takes the logarithm of a dilogarithm argument `X r34/r13` as
//! `log X + log r34 - log r13` with roots `X` normalized by the complex
//! `sqrt(m²)`; with a width the normalization turns `X` across the negative
//! real axis where Denner and Dittmaier's root `u = X sqrt(m²)/|sqrt(m²)|`
//! does not, and their one-mass formula takes `log u + log(Y23/Y02)`. In
//! `4 1.0 P 14.78 12.26 -18.32 15.39 -17.48 0 0 0 0 0 0 3.018 -0.302` OneLOop
//! (both precisions) and this crate jumped at P = -0.055 from -0.0250 to
//! 0.0254 + 1.87i, and at P = -1e-6 gave 6.11 + 383i and 3.03 + 192i for
//! -0.0599 - 8.2e-4i; 131 of 3000 random one-mass boxes with |Im m²| from 1e-3
//! to 0.5 of Re m² were wrong by 0.3 to 7 of the value.
//! References (tests/data/one_mass_width.txt): Denner and Dittmaier's formula
//! at 40 digits, which COLLIER's two branches in quadruple precision and a
//! direct Feynman-parameter integral reproduce, independently of OneLOop.
#![cfg(feature = "generated-evaluators")]
use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend,
};
use symbolica::domains::float::{Complex, DoubleFloat, Float, RealLike};

type C = Complex<f64>;

/// Relative tolerance in binary64 (references are within 1e-16 of the
/// 40-digit values).
const TOLERANCE: f64 = 1e-12;
/// DoubleFloat and 128 bits against the rounded references.
const PRECISE: f64 = 1e-15;
/// Agreement of the 128-bit relabellings of one box.
const IMAGE_AGREEMENT: f64 = 1e-30;
/// Rows 2 and 3 (P = -1e-6 and -2.3e-10) are next to a collinear limit and
/// to `r13 r24 = r12 r34` at once. In the labelling of the data they keep
/// 4e-16 in binary64; relabelled, the small invariant becomes `p2²` or
/// `s12` of the one-mass formula, which kept 2e-10. The loss predates these
/// sheets: the same box with a real mass loses 4e-6 there, and with
/// Im m² = -0.003 these expressions lose 1.3e-7 (native and interpreted)
/// and 2.9e-7 (SymJIT), where 09b2c74, still on the right sheet at that
/// Im m², lost 2.3e-7 (native), 3.8e-6 (interpreted) and 3.0e-7 (SymJIT).
/// The bound 1e-9 is empirical, 5 times the 2e-10 observed. The 128-bit
/// relabellings spread by 2e-29 (the invariant's 1e-10 of the scale times
/// 2^-128).
const COLLINEAR_ROWS: [usize; 2] = [2, 3];
const COLLINEAR_RELABELLED: f64 = 1e-9;
const COLLINEAR_IMAGE_AGREEMENT: f64 = 1e-27;

/// A 0C00 box next to a collinear limit (p2² = 2.4e-10) with the mass next
/// to its shells in p4² and s23 (m² - p4² = -4.4e-5 - 7.3e-5i, m² - s23 =
/// 4e-12 - 7.3e-5i): a root X of the one-mass quadratic is within 7e-8 of
/// 1/r14, and the mass pair's argument X r14 next to one on a wound sheet.
/// Its complement formed from the rounded root lost 1.2e-10 of the value in
/// 8 of the 24 relabellings, on every binary64 route (2.8e-10 SymJIT), where
/// OneLOop's sheets (09b2c74) kept 8e-16; a one-ulp change of any input
/// moves the value by at most 2.7e-15. Reference: native Float at 256 bits,
/// on which 09b2c74 and f2adce4 agree to 1e-73, rounded to binary64; COLLIER
/// 1.2.9 in quadruple precision (COLI and DD agree, fed the exact decimal
/// expansions of the inputs) agrees within 1.2e-18 in every relabelling.
const MASS_PAIR_NEXT_TO_ONE: &str = "4 1.0 18.592973581222886 2.3921721013489905e-10 \
    221.42383370666434 0.04722941995105156 1.1932387948634814e-07 0.04718503480698597 \
    0.0 0.0 0.0 0.0 0.0 0.0 0.0471850348110201 -7.343284330107032e-05 \
    0.28392146259104049021502716246973941336 -0.026104059475840481041527149036706690315 \
    0.0 0.0 0.0 0.0";

fn parse(line: &str) -> (Vec<C>, C) {
    let v: Vec<f64> = line.split_whitespace().map(|x| x.parse().unwrap()).collect();
    assert_eq!(v[0], 4.0);
    assert_eq!(v.len(), 22);
    let mut args: Vec<C> = v[2..8].iter().map(|&x| C::new(x, 0.)).collect();
    args.extend(v[8..16].chunks_exact(2).map(|z| C::new(z[0], z[1])));
    args.push(C::new(v[1], 0.));
    (args, C::new(v[16], v[17]))
}

fn rows() -> Vec<(Vec<C>, C)> {
    include_str!("data/one_mass_width.txt")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(parse)
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
fn complex_one_mass_boxes_match_the_references_in_binary64() {
    on_stack(|| {
        let rows = rows();
        assert_eq!(rows.len(), 15);
        let failures = binary64_failures(&rows, |line, number| {
            if number > 0 && COLLINEAR_ROWS.contains(&line) {
                COLLINEAR_RELABELLED
            } else {
                TOLERANCE
            }
        });
        assert!(failures.is_empty(), "{} failures:\n{}", failures.len(), failures.join("\n"));
    });
}

/// The same routes for [MASS_PAIR_NEXT_TO_ONE], in every relabelling.
#[test]
fn complex_one_mass_box_with_its_mass_pair_next_to_one_matches_the_reference_in_binary64() {
    on_stack(|| {
        let failures = binary64_failures(&[parse(MASS_PAIR_NEXT_TO_ONE)], |_, _| TOLERANCE);
        assert!(failures.is_empty(), "{} failures:\n{}", failures.len(), failures.join("\n"));
    });
}

/// Every relabelling of `rows` on the binary64 routes within `tolerance(row,
/// relabelling)`, and on DoubleFloat within [PRECISE].
fn binary64_failures(
    rows: &[(Vec<C>, C)],
    tolerance: impl Fn(usize, usize) -> f64,
) -> Vec<String> {
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
            let binary64 = tolerance(line, number);
            let mut check = |route: String, value: C, bound: f64| {
                let error = relative(value, *expected);
                if !(error <= bound) {
                    failures.push(format!(
                        "row {line} image {number} {route}: {value:?} vs {expected:?} ({error:e})"
                    ));
                }
            };
            for (handwritten, evaluator) in &mut native {
                let mut output = [C::new(f64::NAN, f64::NAN); 3];
                evaluator.evaluate(&image, &mut output).unwrap();
                check(format!("native (hand-written: {handwritten})"), output[0], binary64);
            }
            for backend in [EvaluationBackend::Expression, EvaluationBackend::SymJit] {
                let mut output = [C::new(f64::NAN, f64::NAN); 3];
                evaluate_with_backend(ScalarIntegral::D0, &image, &mut output, backend)
                    .unwrap();
                check(backend.name().to_string(), output[0], binary64);
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
                PRECISE,
            );
        }
    }
    failures
}

/// Native code and interpreted expressions at 128 bits: every relabelling
/// matches the reference, and all 24 agree with each other to 1e-30.
#[test]
fn complex_one_mass_boxes_match_the_references_at_128_bits() {
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
                    let agreement = if COLLINEAR_ROWS.contains(&line) {
                        COLLINEAR_IMAGE_AGREEMENT
                    } else {
                        IMAGE_AGREEMENT
                    };
                    if !(error <= PRECISE) || !(spread <= agreement) {
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
