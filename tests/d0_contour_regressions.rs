//! Contour-box D0 (OneLOop's boxc) where the line offsets of the vanishing-Gram
//! change (`box_complex::line_product`) and its second contour lost accuracy
//! against 09b2c74, and degenerate boxes next to a Feynman-parameter corner.
//!
//! - The line products replace a pole's offset `y - z` from a line edge's
//!   root by Vieta's product `a1 (Q(u) - ratio B(u))/K`, exact in `a1`. A pole
//!   also lies next to the root when the conic nearly vanishes where the line
//!   meets the edge; there the polynomial cancels at least as much as `y - z`,
//!   and the line's term cancels against the conic's term on the same edge,
//!   which keeps `y - z`. The products are now taken only where `a1 = g + beta
//!   h` itself cancels (a line nearly parallel to the contour direction).
//! - With `p4 = 0` and a massless `m2` next to the line, a pole lies next to
//!   the corner (1, 0), at an `x` next to 1 that binary64 resolves only to eps.
//!   `zero_external` takes the reflected labelling, which puts the corner at
//!   the origin, and `r1` takes the complements `1 - x` of its dilogarithms'
//!   arguments from the poles.
//!
//! References (tests/data/contour_regressions.txt): COLLIER 1.2.9 in
//! quadruple precision, independent of OneLOop and this crate.
#![cfg(feature = "generated-evaluators")]
use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend,
};
use symbolica::domains::float::{Complex, DoubleFloat, Float, RealLike};

type C = Complex<f64>;

/// Relative tolerance in binary64, relative to |D0|.
const BINARY64: f64 = 1e-11;
/// DoubleFloat, 128 and 256 bits against the references: COLLIER's own
/// relabellings of row 2 spread by 7e-16.
const PRECISE: f64 = 1e-15;
/// Agreement of the 256-bit relabellings of one box (and of the 128-bit ones
/// of rows 1 and 2).
const IMAGE_AGREEMENT: f64 = 1e-30;
/// Rows 3 to 7 put a pole within `t²` of a Feynman-parameter corner (`t` the
/// relative distance of their invariant from its on-shell value, one ulp to
/// 1e-4): at 128 bits its offset keeps `2^-128/t²` of itself, so their
/// relabellings agree only to 1e-16 to 1e-29 (row 3: 4e-31 from the corner,
/// 7e-9 of the offset). Row 3 at 128 bits is within 1.2e-11 (native) and
/// 3.8e-11 (interpreted) of the reference, beyond 1e-15 in 16 relabellings
/// (in 20 at 09b2c74, with the same bounds): Float tracks the precision of
/// the poles' small components through their cancellations in the labelling
/// that puts the corner at the origin (2.4e-7 at 09b2c74 in that labelling
/// with `p3` 3e-16, relative, below `s12`).
const ROW_3_AT_128_BITS: f64 = 1e-10;

fn rows(data: &str) -> Vec<(Vec<C>, C)> {
    data.lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| {
            let v: Vec<f64> = line
                .split_whitespace()
                .map(|x| x.parse().unwrap())
                .collect();
            assert_eq!(v[0], 4.0);
            assert_eq!(v.len(), 22);
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
    let error =
        (value.re - expected.re).hypot(value.im - expected.im) / expected.re.hypot(expected.im);
    if error.is_nan() { f64::INFINITY } else { error }
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
fn check_binary64(data: &str, count: usize, tolerance: f64) {
    let rows = rows(data);
    assert_eq!(rows.len(), count);
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
            let mut check = |route: String, value: C, tolerance: f64| {
                let error = relative(value, *expected);
                if !(error <= tolerance) {
                    failures.push(format!(
                        "row {} image {number} {route}: {value:?} vs {expected:?} ({error:e})",
                        line + 1
                    ));
                }
            };
            for (handwritten, evaluator) in &mut native {
                let mut output = [C::new(f64::NAN, f64::NAN); 3];
                evaluator.evaluate(&image, &mut output).unwrap();
                check(
                    format!("native (hand-written: {handwritten})"),
                    output[0],
                    tolerance,
                );
            }
            for backend in [EvaluationBackend::Expression, EvaluationBackend::SymJit] {
                let mut output = [C::new(f64::NAN, f64::NAN); 3];
                evaluate_with_backend(ScalarIntegral::D0, &image, &mut output, backend).unwrap();
                check(backend.name().to_string(), output[0], tolerance);
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

/// Every binary64 route within 1e-11 of the references in all 24
/// relabellings (now at most 7.7e-12), and DoubleFloat within 1e-15. With the
/// line products everywhere (9a898ee): rows 1 and 2 up to 3.9e-10 and 1.7e-11
/// (2.9e-12 and 5.8e-13 at 09b2c74), row 3 NaN in 12 relabellings (wrong by
/// up to 130 times at 09b2c74), rows 4 to 7 wrong by up to 170 times (71 at
/// 09b2c74).
#[test]
fn contour_box_regressions_match_collier_in_binary64() {
    on_stack(|| check_binary64(include_str!("data/contour_regressions.txt"), 7, BINARY64));
}

/// Native code and interpreted expressions at 128 and 256 bits: every
/// relabelling matches the reference (row 3 at 128 bits within
/// [ROW_3_AT_128_BITS]), and all 24 agree with each other to 1e-30 at 256 bits
/// (and at 128 bits for rows 1 and 2, see [ROW_3_AT_128_BITS]).
#[test]
fn contour_box_regressions_match_collier_at_high_precision() {
    on_stack(|| {
        let mut failures = Vec::new();
        for bits in [128, 256] {
            let mut native =
                NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::D0, bits).unwrap();
            let mut expression = PrecisionEvaluator::with_binary_precision_and_backend(
                ScalarIntegral::D0,
                bits,
                EvaluationBackend::Expression,
            )
            .unwrap();
            for (line, (args, expected)) in rows(include_str!("data/contour_regressions.txt"))
                .iter()
                .enumerate()
            {
                for route in ["native", "expression"] {
                    let mut values: Vec<Complex<Float>> = Vec::new();
                    for image in images(args) {
                        let input: Vec<_> = image
                            .iter()
                            .map(|z| {
                                Complex::new(
                                    Float::with_val(bits, z.re),
                                    Float::with_val(bits, z.im),
                                )
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
                        let error =
                            relative(C::new(value.re.to_f64(), value.im.to_f64()), *expected);
                        let spread = (value.re.clone() - &first.re)
                            .to_f64()
                            .hypot((value.im.clone() - &first.im).to_f64())
                            / scale;
                        let (tolerance, agreement) = match (bits, line + 1) {
                            (128, 3) => (ROW_3_AT_128_BITS, f64::INFINITY),
                            (128, 4..) => (PRECISE, f64::INFINITY),
                            _ => (PRECISE, IMAGE_AGREEMENT),
                        };
                        if !(error <= tolerance) || !(spread <= agreement) {
                            failures.push(format!(
                                "row {} image {number} {route} {bits} bits: {} {} vs {expected:?} \
                                 (error {error:e}, spread from image 0 {spread:e})",
                                line + 1,
                                value.re,
                                value.im
                            ));
                        }
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

/// Boxes binary64 still gets wrong (references: COLLIER 1.2.9 in quadruple
/// precision, COLI branch):
///
/// - A Gram determinant at the rounding level with four equal complex masses
///   (`tests/generate_near_threshold_rows.py`, seed 201): the selected
///   contour's `a1 = g + beta h` rounds to zero, where it is 1.5e-13 at 256
///   bits, and the poles, real in binary64, lose the imaginary parts of order
///   `a1 Im m²` that choose their lips. 5.1e-3 on the hand-written and
///   interpreted routes in every relabelling (0.34 at 09b2c74, 0.55 with the
///   line products everywhere); the generated code and SymJIT keep 8e-14.
/// - Three real masses, `p1²` one ulp below `m2² = m4² = p3²`, the massless
///   `m1` between the legs `p1` and `p4`, and `s23 = 0`: in the labellings
///   the contour box takes, the conic's value and gradient vanish at the
///   massless corner to within rounding. Every binary64 route was wrong in
///   every relabelling at 09b2c74 (NaN in 8, by 1.4e-2 in 4 and by 0.3 in
///   12); now NaN in 16 and 1.4e-2 in 8.
const KNOWN_FAILURES: &str = "\
4 1.0 2.7819050504399936 2.3889525049365044 2.388952504936503 105.83204530964998 51.91801797559696 9.555810019746128 2.2406833551152996 -0.0005207803841710876 2.2406833551152996 -0.0005207803841710876 2.2406833551152996 -0.0005207803841710876 2.2406833551152996 -0.0005207803841710876 -1.83619824164963746073e-2 -6.24123636690162522126e-3 0.0 0.0 0.0 0.0
4 1.0 0.1113435943976259 108.23851883958666 0.11134359439762591 -0.012515862621692797 9.0 0.0 0.0 0.0 0.11134359439762591 0.0 9.0 0.0 0.11134359439762591 0.0 -2.77080305527435389896e+1 2.11066704520210922998e+1 0.0 0.0 0.0 0.0
";

#[test]
#[ignore = "known binary64 failures: a Gram determinant at the rounding level, a conic singular at a massless corner"]
fn contour_box_known_failures_match_collier_in_binary64() {
    on_stack(|| check_binary64(KNOWN_FAILURES, 2, BINARY64));
}
