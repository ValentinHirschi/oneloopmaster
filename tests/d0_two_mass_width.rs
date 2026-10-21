//! Two-mass D0 with complex masses (OneLOop's boxf2, and boxf5 for opposite
//! masses), in every relabelling. OneLOop takes the logarithm of a
//! dilogarithm argument `X r34/r13` as `log X + log r34 - log r13` with roots
//! `X` normalized by the complex square roots of the masses; with widths the
//! normalization turns `X` across the negative real axis where Denner and
//! Dittmaier's root `u = X sqrt(m4²)/sqrt(m3²)` does not, and on that sheet
//! the dilogarithm of an argument next to one picks the wrong side of its
//! branch point (by `4 pi²` in the terms). In `4 1.0 0.0462 -4.297 12.27
//! 10.74 -31.21 -13.06 0 0 0 0 0.2161 -0.00173 1.556 -0.00161` OneLOop (both
//! precisions) and this crate gave 0.0782 + 0.0831i for -7.16e-4 + 0.0831i;
//! 0.4% to 2.8% of random two-mass boxes with |Im m²| from 1e-3 to 0.5 of
//! Re m², depending on how the masses and invariants are sampled (STATUS.md),
//! were wrong by 0.6 to 17 of the value, in every mass pattern, with
//! one real mass, with the widths divided by 1e6, and next to collinear limits
//! (by up to 9e4).
//! References (tests/data/two_mass_width.txt): Denner and Dittmaier's
//! two-mass formula at 40 digits, which COLLIER's two branches in quadruple
//! precision reproduce, independently of OneLOop.
//!
//! On those sheets the mass pairs' arguments `X r14` and `X/r14` reach one
//! where both massive lines sit at or next to their shells (a real mass and
//! a complex one whose real part is), and `1 - X r14` formed from the rounded
//! argument lost up to 0.16 of the value in binary64, where OneLOop's sheets
//! kept 6e-16 (tests/data/two_mass_shell.txt, references from the same
//! formula, which COLLIER's COLI branch reproduces).
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

const WIDTH_ROWS: &str = include_str!("data/two_mass_width.txt");
const SHELL_ROWS: &str = include_str!("data/two_mass_shell.txt");
const CONDITIONING_ROWS: &str = include_str!("data/two_mass_conditioning.txt");

/// A data row: the inputs, the finite part rounded to binary64 and, where
/// given, its real and imaginary parts to 90 digits.
struct Row {
    args: Vec<C>,
    expected: C,
    digits: Option<[String; 2]>,
}

fn rows(data: &str) -> Vec<Row> {
    data.lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split_whitespace().collect();
            assert!(fields.len() == 22 || fields.len() == 24);
            let v: Vec<f64> = fields[..22].iter().map(|x| x.parse().unwrap()).collect();
            assert_eq!(v[0], 4.0);
            let mut args: Vec<C> = v[2..8].iter().map(|&x| C::new(x, 0.)).collect();
            args.extend(v[8..16].chunks_exact(2).map(|z| C::new(z[0], z[1])));
            args.push(C::new(v[1], 0.));
            Row {
                args,
                expected: C::new(v[16], v[17]),
                digits: (fields.len() == 24).then(|| [fields[22].into(), fields[23].into()]),
            }
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
fn complex_two_mass_boxes_match_the_references_in_binary64() {
    on_stack(|| {
        let rows = rows(WIDTH_ROWS);
        assert_eq!(rows.len(), 17);
        check_binary64(&rows);
    });
}

/// The same routes for boxes with both massive lines at or next to their
/// shells (tests/data/two_mass_shell.txt).
#[test]
fn complex_two_mass_boxes_next_to_their_shells_match_the_references_in_binary64() {
    on_stack(|| {
        let rows = rows(SHELL_ROWS);
        assert_eq!(rows.len(), 9);
        check_binary64(&rows);
    });
}

/// The same routes for well-conditioned boxes whose Denner-Dittmaier terms
/// cancelled (tests/data/two_mass_conditioning.txt): a mass pair's argument
/// next to one with the other pair's complement between 1/16 and 1/4, and
/// roots on opposite lips of the negative real axis taken across it from the
/// farther one.
#[test]
fn complex_two_mass_boxes_with_cancelling_sheet_terms_match_the_references_in_binary64() {
    on_stack(|| {
        let rows = rows(CONDITIONING_ROWS);
        assert_eq!(rows.len(), 4);
        check_binary64(&rows);
    });
}

fn check_binary64(rows: &[Row]) {
    let mut failures = Vec::new();
    let mut native: Vec<_> = [true, false]
        .map(|handwritten| {
            let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
            evaluator.set_handwritten(handwritten);
            (handwritten, evaluator)
        })
        .into();
    let mut double = NativeEvaluator::<DoubleFloat>::new(ScalarIntegral::D0).unwrap();
    for (line, Row { args, expected, .. }) in rows.iter().enumerate() {
        for (number, image) in images(args).into_iter().enumerate() {
            let mut check = |route: String, value: C, tolerance: f64| {
                let error = relative(value, *expected);
                if !(error <= tolerance) {
                    failures.push(format!(
                        "row {line} image {number} {route}: {value:?} vs {expected:?} ({error:e})"
                    ));
                }
            };
            for (handwritten, evaluator) in &mut native {
                let mut output = [C::new(f64::NAN, f64::NAN); 3];
                evaluator.evaluate(&image, &mut output).unwrap();
                check(
                    format!("native (hand-written: {handwritten})"),
                    output[0],
                    TOLERANCE,
                );
            }
            for backend in [EvaluationBackend::Expression, EvaluationBackend::SymJit] {
                let mut output = [C::new(f64::NAN, f64::NAN); 3];
                evaluate_with_backend(ScalarIntegral::D0, &image, &mut output, backend).unwrap();
                check(backend.name().to_string(), output[0], TOLERANCE);
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
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// Native code and interpreted expressions at 128 bits: every relabelling
/// matches the reference, and all 24 agree with each other to 1e-30.
#[test]
fn complex_two_mass_boxes_match_the_references_at_128_bits() {
    on_stack(|| check_precise(&rows(WIDTH_ROWS), 128, Some(IMAGE_AGREEMENT)));
}

/// The boxes next to their shells, native and interpreted at 128 and 256
/// bits, against the 90-digit references. Symbolica's `Float` tracks
/// precision: the difference of an invariant and a mass next to it keeps only
/// its significant bits, and the result reports the precision that is left
/// (down to 80 of 128 and 208 of 256 bits here; see
/// tests/permutation_regressions.rs). The input-cancellation floor: a
/// difference of two binary64 inputs cancels at most 53 bits, so every
/// relabelling must track at least `bits - 56` bits and match the reference
/// to `2^(56 - bits)` (at least `bits - 56` accurate bits), and the binary64
/// reference to 1e-15.
///
/// The value is not checked against the precision it tracks: numerica's
/// tracked precision is an estimate, not a bound. On row 0 the tracked bits
/// exceed the accurate ones by up to 3.5, 8.4, 8.3 and 7.1 bits at 128, 192,
/// 256 and 320 bits, against a 1024-bit reference. At 09b2c74 (before the
/// two-mass sheets and the single-rounded products) row 0 falls at most 5.85
/// bits short, but row 6 (native) 8.23, 9.09, 9.63 and 7.75 bits (at 256
/// bits 253 tracked, 243.37 accurate); the rows' 90-digit references cannot
/// measure more than about 299 accurate bits. The minimum accuracy over the
/// relabellings is unchanged: at least 202.5 bits at 256 (202.3 at 09b2c74)
/// and 267.0 at 320 (266.9). Single relabellings lose up to 40.3 bits against
/// 09b2c74 (row 0), all above the floor and the same since f2adce4.
#[test]
fn complex_two_mass_boxes_next_to_their_shells_match_the_references_at_128_and_256_bits() {
    on_stack(|| {
        let rows = rows(SHELL_ROWS);
        check_precise(&rows, 128, None);
        check_precise(&rows, 256, None);
    });
}

/// `agreement`: the relabellings' bound on their spread from the first; none
/// for the input-cancellation floor against the rows' 90-digit references.
fn check_precise(rows: &[Row], bits: u32, agreement: Option<f64>) {
    let mut failures = Vec::new();
    let mut native =
        NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::D0, bits).unwrap();
    let mut expression = PrecisionEvaluator::with_binary_precision_and_backend(
        ScalarIntegral::D0,
        bits,
        EvaluationBackend::Expression,
    )
    .unwrap();
    for (line, row) in rows.iter().enumerate() {
        let expected = &row.expected;
        for route in ["native", "expression"] {
            let mut values: Vec<Complex<Float>> = Vec::new();
            for image in images(&row.args) {
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
                let error = relative(C::new(value.re.to_f64(), value.im.to_f64()), *expected);
                let (bound, deviation) = match agreement {
                    Some(agreement) => {
                        let spread = (value.re.clone() - &first.re)
                            .to_f64()
                            .hypot((value.im.clone() - &first.im).to_f64())
                            / scale;
                        (agreement, spread)
                    }
                    None => {
                        let [re, im] = row
                            .digits
                            .as_ref()
                            .unwrap()
                            .each_ref()
                            .map(|digits| Float::parse(digits, Some(bits + 64)).unwrap());
                        let tracked = value.re.prec().min(value.im.prec());
                        let precise = (value.re.clone() - &re)
                            .to_f64()
                            .hypot((value.im.clone() - &im).to_f64())
                            / re.to_f64().hypot(im.to_f64());
                        // The input-cancellation floor: at least `bits - 56`
                        // tracked and accurate bits (not `tracked` accurate
                        // bits; see the test's comment).
                        let bound = if tracked + 56 >= bits {
                            2f64.powi(56 - bits as i32)
                        } else {
                            f64::NAN
                        };
                        (bound, precise)
                    }
                };
                if !(error <= PRECISE) || !(deviation <= bound) {
                    failures.push(format!(
                        "row {line} image {number} {route} {bits} bits: {} {} vs \
                             {expected:?} (error {error:e}, deviation {deviation:e}, bound \
                             {bound:e}; tracked bits {}, {})",
                        value.re,
                        value.im,
                        value.re.prec(),
                        value.im.prec()
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
}
