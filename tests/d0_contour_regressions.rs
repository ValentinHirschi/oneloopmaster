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
//!   the origin, unless one of the reflected chart's conics has a double
//!   contour direction, and `r1` takes the complements `1 - x` of its
//!   dilogarithms' arguments from the poles.
//! - With the Gram determinant at the rounding level one contour's `a1`
//!   vanishes and the other's `K` with it (tests/data/contour_gram_rounding.txt):
//!   a `K` that rounds to zero puts its far pole 2^52 beyond the near one, and
//!   where the selected contour's `a1` keeps fewer than 24 bits the regular
//!   chart of the other contour is evaluated too and taken unless the two
//!   agree to binary64's half precision (`box_complex::tfun_chart`). Where
//!   no chart is regular, binary64 returns NaN unless the chart used keeps
//!   10 bits under a move of its line by one rounding unit (a flat chart) or
//!   13 bits with its K at the rounding level (a K that rounds to zero),
//!   instead of a silently wrong value (tests/data/contour_unresolved.txt).
//!   Where the other chart's K is itself a rounding residue, the selected
//!   chart is kept where K at its rounding level moves the other more than
//!   2^10 times as much as the selected one moves with its line, if it keeps
//!   20 bits under that move (rows 1 and 3 of tests/data/contour_fresh_sweep.txt).
//! - At an exactly double real pole of a line edge no line product is taken
//!   (`box_complex::r1_line`).
//! - Where a line product disagrees with the rounded offset of a pole that
//!   carries the error itself, the line's and the conic's terms share that
//!   error and cancel it with `y - z`, which the product breaks; binary64
//!   probes such a pole and keeps the variant that moves less
//!   (`box_complex::tfun_selected`, tests/data/contour_pole_errors.txt).
//! - The line products are kept only where every edge whose offset cancels
//!   has one, and not on an edge whose pole lies within 8 offsets of an end
//!   of the edge: there `y - z` keeps the rounded pole's cancellations that
//!   the products break (`box_complex::products_kept`,
//!   tests/data/contour_line_products.txt). A product more than 2^10 times
//!   `a1` from its root is probed as one that disagrees with its pole (row 2
//!   of tests/data/contour_fresh_sweep.txt).
//! - The reflection of `zero_external` needs the line nearer (1, 0) than
//!   (0, 0) by more than 1/128: where both corners are equally near it moves
//!   no pole away from a corner, and its chart lost up to 1.2 times the value.
//!   It is taken only where a `tfun` of the original chart selects a pole that
//!   the rounding resolves poorly next to the corner's vertex: elsewhere the
//!   reflected chart has nothing to gain and lost up to 1e-11 to its own
//!   cancellations. Nor is it taken where the reflected chart selects a flat
//!   contour without a regular alternative (in `Float` any flat contour) and
//!   the original selects none (rows 4 and 5 of
//!   tests/data/contour_fresh_sweep.txt).
//!
//! References (tests/data/contour_regressions.txt,
//! tests/data/contour_gram_rounding.txt, tests/data/contour_unresolved.txt,
//! tests/data/contour_pole_errors.txt, tests/data/contour_line_products.txt,
//! tests/data/contour_fresh_sweep.txt): COLLIER 1.2.9 in quadruple precision,
//! independent of OneLOop and this crate.
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

/// The rows' arguments, references and, where a 23rd field gives one, their
/// own binary64 tolerance.
fn rows(data: &str) -> Vec<(Vec<C>, C, Option<f64>)> {
    data.lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| {
            let v: Vec<f64> = line
                .split_whitespace()
                .map(|x| x.parse().unwrap())
                .collect();
            assert_eq!(v[0], 4.0);
            assert!(v.len() == 22 || v.len() == 23);
            let mut args: Vec<C> = v[2..8].iter().map(|&x| C::new(x, 0.)).collect();
            args.extend(v[8..16].chunks_exact(2).map(|z| C::new(z[0], z[1])));
            args.push(C::new(v[1], 0.));
            (args, C::new(v[16], v[17]), v.get(22).copied())
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

/// The data rows `range` of `data` (comments skipped).
fn data_lines(data: &str, range: std::ops::Range<usize>) -> String {
    data.lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .skip(range.start)
        .take(range.len())
        .map(|line| format!("{line}\n"))
        .collect()
}

/// Binary64 routes (native with and without the hand-written kernels, the
/// interpreted expressions and SymJIT) and DoubleFloat, in every relabelling;
/// with `nonfinite`, a non-finite binary64 value passes (a box whose charts
/// binary64 cannot resolve must not be silently wrong).
fn check_binary64(
    data: &str,
    count: usize,
    tolerance: f64,
    double_tolerance: f64,
    nonfinite: bool,
) {
    let all: Vec<usize> = (0..24).collect();
    check_binary64_images(data, count, &all, tolerance, double_tolerance, nonfinite)
}

/// [check_binary64] in the relabellings `selected` of [images] only.
fn check_binary64_images(
    data: &str,
    count: usize,
    selected: &[usize],
    tolerance: f64,
    double_tolerance: f64,
    nonfinite: bool,
) {
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
    for (line, (args, expected, own)) in rows.iter().enumerate() {
        let tolerance = own.unwrap_or(tolerance);
        for (number, image) in images(args).into_iter().enumerate() {
            if !selected.contains(&number) {
                continue;
            }
            let mut check = |route: String, value: C, tolerance: f64| {
                let error = relative(value, *expected);
                let detected = nonfinite && !(value.re.is_finite() && value.im.is_finite());
                if !(error <= tolerance) && !detected {
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
            let value = C::new(output[0].re.to_f64(), output[0].im.to_f64());
            let error = relative(value, *expected);
            if !(error <= double_tolerance) {
                failures.push(format!(
                    "row {} image {number} DoubleFloat: {value:?} vs {expected:?} ({error:e})",
                    line + 1
                ));
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

/// Every binary64 route within 1e-11 of the references in all 24
/// relabellings (now at most 7.7e-12), and DoubleFloat within 1e-15. With the
/// line products everywhere (9a898ee): rows 1 and 2 up to 3.9e-10 and 1.7e-11
/// (2.9e-12 and 5.8e-13 at 09b2c74), row 3 NaN in 12 relabellings (wrong by
/// up to 130 times at 09b2c74), rows 4 to 7 wrong by up to 170 times (71 at
/// 09b2c74).
#[test]
fn contour_box_regressions_match_collier_in_binary64() {
    on_stack(|| {
        check_binary64(
            include_str!("data/contour_regressions.txt"),
            7,
            BINARY64,
            PRECISE,
            false,
        )
    });
}

/// Native code and interpreted expressions at 128 and 256 bits: every
/// relabelling matches the reference (row 3 at 128 bits within
/// [ROW_3_AT_128_BITS]), and all 24 agree with each other to 1e-30 at 256 bits
/// (and at 128 bits for rows 1 and 2, see [ROW_3_AT_128_BITS]).
#[test]
fn contour_box_regressions_match_collier_at_high_precision() {
    on_stack(|| {
        check_precise(
            include_str!("data/contour_regressions.txt"),
            |bits, row| match (bits, row) {
                (128, 3) => (ROW_3_AT_128_BITS, f64::INFINITY),
                (128, 4..) => (PRECISE, f64::INFINITY),
                _ => (PRECISE, IMAGE_AGREEMENT),
            },
        )
    });
}

/// Every binary64 route within 1e-11 of the references in all 24
/// relabellings (row 2 within 1e-8 and row 9 within 1e-10, see the data
/// file), where 09b2c74 was non-finite in 12 to 20 relabellings of rows 1, 2,
/// 4 and 5 and wrong by 1.4e-2 to 5.7 times the value on the others, and
/// f2adce4 by 5.1e-3 to 140 times on all but row 2 (4.7e-3); DoubleFloat
/// within 1e-14 (2.3e-15 on row 6, as at 09b2c74).
#[test]
fn contour_boxes_with_a_rounding_level_gram_determinant_match_collier_in_binary64() {
    on_stack(|| {
        check_binary64(
            include_str!("data/contour_gram_rounding.txt"),
            9,
            BINARY64,
            1e-14,
            false,
        )
    });
}

/// Boxes whose binary64 chart the rounding cannot resolve (rows 1 and 2 of
/// tests/data/contour_unresolved.txt): every binary64 route non-finite (the
/// caller detects it and rescues it in higher precision) or within 1e-11 in
/// all 24 relabellings. Now NaN in all of them; 09b2c74 too, and f2adce4 was
/// finite and wrong by up to 6.5e-4 and 1.9 on every route. DoubleFloat within
/// 1e-14 (4.3e-15 on row 2).
#[test]
fn contour_boxes_with_unresolved_charts_are_not_silently_wrong() {
    on_stack(|| {
        let data = data_lines(include_str!("data/contour_unresolved.txt"), 0..2);
        check_binary64(&data, 2, BINARY64, 1e-14, true)
    });
}

/// Boxes with a flat chart that keeps its value under the move of its line
/// (rows 3 and 4 of tests/data/contour_unresolved.txt): every binary64 route
/// finite and within 1e-11 in all 24 relabellings (now at most 2.3e-12; 2.3e-12
/// at 09b2c74), DoubleFloat within 1e-14.
#[test]
fn contour_boxes_with_resolved_flat_charts_stay_finite() {
    on_stack(|| {
        let data = data_lines(include_str!("data/contour_unresolved.txt"), 2..4);
        check_binary64(&data, 2, BINARY64, 1e-14, false)
    });
}

/// The same boxes in native code and interpreted expressions at 128 and 256
/// bits: within 1e-15 of the references, the 24 relabellings within 1e-30 of
/// each other at 256 bits.
#[test]
fn contour_boxes_with_a_rounding_level_gram_determinant_match_collier_at_high_precision() {
    on_stack(|| {
        check_precise(include_str!("data/contour_gram_rounding.txt"), |bits, _| {
            if bits == 256 {
                (PRECISE, IMAGE_AGREEMENT)
            } else {
                (PRECISE, f64::INFINITY)
            }
        })
    });
}

/// Boxes of the second verifier's fresh sweep of 59b96e5 (rows 1 to 3 of
/// tests/data/contour_fresh_sweep.txt): every binary64 route within the row's
/// tolerance in all 24 relabellings, DoubleFloat within 1e-14. Row 1 (an exact
/// vanishing Gram determinant whose other chart's K rounds to zero) was wrong
/// by 5.3e-6 in all 24 at 59b96e5 and keeps 1.0e-7 hand-written and 2.5e-8 on
/// the other routes, as 09b2c74; row 2 (one line edge's offset 1.3e4 times
/// `a1`) was wrong by up to 1.4e-5 and keeps 6.4e-10, 2.4e-14 in the 12
/// reflected relabellings; row 3 (an exact vanishing Gram determinant whose
/// other chart's K is a rounding residue) was NaN and keeps 6.8e-9, as 09b2c74.
#[test]
fn contour_boxes_of_a_fresh_sweep_match_collier_in_binary64() {
    on_stack(|| {
        let data = data_lines(include_str!("data/contour_fresh_sweep.txt"), 0..3);
        check_binary64(&data, 3, BINARY64, 1e-14, false)
    });
}

/// The relabellings of row 4 of tests/data/contour_fresh_sweep.txt that
/// `zero_external` may reflect (see the data file).
const REFLECTION_UNRESOLVED_IMAGES: [usize; 12] = [0, 1, 2, 9, 10, 11, 15, 16, 17, 18, 19, 20];

/// A box whose reflected zero-external chart selects a flat contour without a
/// regular alternative while the original selects none (row 4 of
/// tests/data/contour_fresh_sweep.txt): in the relabellings that reach the
/// reflection test, every binary64 route within 1e-9 (now and at 09b2c74
/// 9.7e-11; reflected, 59b96e5 was wrong by 2.9e-2 or NaN), DoubleFloat within
/// 1e-14. The other 12 relabellings, which put the box into that orientation
/// directly, are wrong by 2.8e-2 or NaN, as 09b2c74 (3.4e-2 to 0.23).
#[test]
fn zero_external_keeps_its_orientation_where_the_reflected_chart_is_unresolved() {
    on_stack(|| {
        let data = data_lines(include_str!("data/contour_fresh_sweep.txt"), 3..4);
        check_binary64_images(&data, 1, &REFLECTION_UNRESOLVED_IMAGES, BINARY64, 1e-14, false)
    });
}

/// The relabellings of row 5 of tests/data/contour_fresh_sweep.txt that keep
/// their orientation at 128 bits (see the data file).
const FLOAT_REFLECTION_IMAGES: [usize; 12] = [0, 1, 4, 8, 9, 11, 13, 16, 17, 18, 20, 21];

/// A box with an exact vanishing Gram determinant whose reflected
/// zero-external chart selects a flat contour (row 5 of
/// tests/data/contour_fresh_sweep.txt): binary64 within the row's tolerance in
/// all 24 relabellings (5.5e2 at 09b2c74), and native code and interpreted
/// expressions at 128 bits within 1e-15 in the relabellings that reach the
/// reflection test, where `Float` keeps the original orientation (3.3e-16 at
/// 64 tracked bits, as 09b2c74; 59b96e5 reflected and kept 20 tracked bits,
/// wrong by up to 1.9e-5). The other 12 keep 1.2e-6 to 1.9e-5 at 20 tracked
/// bits (09b2c74 3.7e-2 to 1.4e6 at 2 to 11).
#[test]
fn zero_external_keeps_its_orientation_in_float_where_the_reflected_contour_is_flat() {
    on_stack(|| {
        let data = data_lines(include_str!("data/contour_fresh_sweep.txt"), 4..5);
        check_binary64(&data, 1, BINARY64, 1e-14, false);
        let (args, expected, _) = rows(&data).remove(0);
        let bits = 128;
        let mut native =
            NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::D0, bits).unwrap();
        let mut expression = PrecisionEvaluator::with_binary_precision_and_backend(
            ScalarIntegral::D0,
            bits,
            EvaluationBackend::Expression,
        )
        .unwrap();
        let mut failures = Vec::new();
        for (number, image) in images(&args).into_iter().enumerate() {
            if !FLOAT_REFLECTION_IMAGES.contains(&number) {
                continue;
            }
            let input: Vec<_> = image
                .iter()
                .map(|z| Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im)))
                .collect();
            for route in ["native", "expression"] {
                let mut output = core::array::from_fn::<_, 3, _>(|_| {
                    Complex::new(Float::new(bits), Float::new(bits))
                });
                if route == "native" {
                    native.evaluate(&input, &mut output).unwrap();
                } else {
                    expression.evaluate(&input, &mut output).unwrap();
                }
                let value = C::new(output[0].re.to_f64(), output[0].im.to_f64());
                let error = relative(value, expected);
                if !(error <= PRECISE) {
                    failures.push(format!("image {number} {route}: {value:?} vs {expected:?} ({error:e})"));
                }
            }
        }
        assert!(failures.is_empty(), "{} failures:\n{}", failures.len(), failures.join("\n"));
    });
}

/// Native code and interpreted expressions at 128 and 256 bits, in every
/// relabelling: `bounds(bits, row)` gives the tolerance against the
/// reference and the bound on the spread from the first relabelling.
fn check_precise(data: &str, bounds: impl Fn(u32, usize) -> (f64, f64)) {
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
        for (line, (args, expected, _)) in rows(data).iter().enumerate() {
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
                let (tolerance, agreement) = bounds(bits, line + 1);
                for (number, value) in values.iter().enumerate() {
                    let error = relative(C::new(value.re.to_f64(), value.im.to_f64()), *expected);
                    let spread = (value.re.clone() - &first.re)
                        .to_f64()
                        .hypot((value.im.clone() - &first.im).to_f64())
                        / scale;
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
}

/// Boxes whose line products disagree with the rounded offsets of a pole that
/// carries the error itself (tests/data/contour_pole_errors.txt): every
/// binary64 route within the row's tolerance in all 24 relabellings, where
/// f2adce4 lost up to 1.2e-10 (see the data file); DoubleFloat within 1e-14.
#[test]
fn contour_boxes_with_inaccurate_poles_match_collier_in_binary64() {
    on_stack(|| {
        check_binary64(
            include_str!("data/contour_pole_errors.txt"),
            5,
            BINARY64,
            1e-14,
            false,
        )
    });
}

/// Boxes whose line products lost to the rounded `y - z` where the poles'
/// offsets from the line edges' roots do not scale with `a1`
/// (tests/data/contour_line_products.txt, `box_complex::products_kept`):
/// every binary64 route within the row's tolerance in all 24 relabellings,
/// where cb5370d lost up to 1.3e-5 (about 10 times 09b2c74 on rows 2 to 7,
/// see the data file); DoubleFloat within 1e-14. Row 8 keeps every product of
/// its flat charts (SymJIT was wrong by 1.2 there without them).
#[test]
fn contour_boxes_with_unscaled_line_offsets_match_collier_in_binary64() {
    on_stack(|| {
        check_binary64(
            include_str!("data/contour_line_products.txt"),
            8,
            BINARY64,
            1e-14,
            false,
        )
    });
}

/// A box of the NinjaRS streams' oracle sweep (reference: COLLIER 1.2.9 in
/// quadruple precision, COLI branch, which this crate's native 256-bit value
/// matches) whose contour chart has an exactly double real pole on a line
/// edge in the relabellings below (`y1 = y2 = -6.131598e-3`, `a1` cancelling by
/// 2.6e-3, `K` nonzero). There the line product, the square of one offset, over
/// the other offset gave 7.6e18 times the value on the interpreted route at
/// f2adce4 and on the generated route at 82c87d1, where 09b2c74 and the
/// other routes were non-finite; they must stay non-finite or right on the
/// hand-written, generated and interpreted routes. Not tested, and as at
/// 09b2c74: SymJIT is wrong by 3.4e2 to 3.7e2 in these four relabellings,
/// and every binary64 route by 1 to 20 times the value in 12 of the other 20.
const DOUBLE_REAL_POLE: &str = "\
4 1.0 0.0 -0.8946129797968647 -24.09653762111702 0.08938220565009944 0.1941275888116677 0.08938220565009944 0.6383507753193867 0.0 0.6383507753193867 0.0 0.0 0.0 0.25 0.0 0.186223741107587548155091032027640752 0.0 0.0 0.0 0.0 0.0
";

#[test]
fn contour_box_with_a_double_real_pole_is_not_silently_wrong() {
    on_stack(|| {
        let (args, expected, _) = rows(DOUBLE_REAL_POLE).remove(0);
        let images = images(&args);
        let mut native: Vec<_> = [true, false]
            .map(|handwritten| {
                let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
                evaluator.set_handwritten(handwritten);
                (handwritten, evaluator)
            })
            .into();
        let mut failures = Vec::new();
        for number in [0, 9, 16, 18] {
            let image = &images[number];
            let mut values = Vec::new();
            for (handwritten, evaluator) in &mut native {
                let mut output = [C::new(f64::NAN, f64::NAN); 3];
                evaluator.evaluate(image, &mut output).unwrap();
                values.push((format!("native (hand-written: {handwritten})"), output[0]));
            }
            let mut output = [C::new(f64::NAN, f64::NAN); 3];
            evaluate_with_backend(
                ScalarIntegral::D0,
                image,
                &mut output,
                EvaluationBackend::Expression,
            )
            .unwrap();
            values.push((EvaluationBackend::Expression.name().to_string(), output[0]));
            for (route, value) in values {
                let finite = value.re.is_finite() && value.im.is_finite();
                let error = relative(value, expected);
                if finite && !(error <= BINARY64) {
                    failures.push(format!(
                        "image {number} {route}: {value:?} vs {expected:?} ({error:e})"
                    ));
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

/// Boxes of the verifier's sweep where `zero_external` may reflect the chart
/// (a massless `m2` next to the line), in their own labelling (references:
/// COLLIER 1.2.9 in quadruple precision, COLI and DD agreeing on the first,
/// COLI and this crate's native 256-bit value on the others). In the first
/// (v1 15934) the line passes the corners (0, 0) and (1, 0) equally near
/// within rounding, and the reflected chart was wrong by 0.26 to 1.2 times the
/// value on every binary64 route (82c87d1 to 5c996fb), where the original
/// orientation keeps 5.4e-14 (09b2c74: 3.4e-6 to 0.37). In the others the line
/// is 2 % and 3 % nearer (1, 0) (v1 12400, 12507), and the reflection keeps
/// 1.3e-13 and 1.2e-15 where the original orientation loses 1.4e-6 and 2.3e-7.
/// Every binary64 route within 1e-12, DoubleFloat within 1e-14.
const REFLECTION_MARGIN: &str = "\
4 1.0 -9.260750350034394 0.01332834072541728 0.01804580527802355 0.0 22.25741310392177 -31.4907142762668 2.652287073296478 -0.00027624146126022564 0.0 0.0 0.643732097409299 -0.0019642405295094835 2.652287073296478 -0.00027624146126022564 -0.00198082435694084673454364950341565646 0.0287927817444659747843837737425150650 0.0 0.0 0.0 0.0
4 1.0 6.712461550561756 6.71246155001818 6.712461550562316 0.0 -0.019320419936932016 0.0909621142299991 6.712461550562316 0.0 0.0 0.0 6.712461550562316 0.0 6.712461550562316 0.0 0.526261206935409259760819801758524851 0.0 0.0 0.0 0.0 0.0
4 1.0 0.0 0.5021464029637546 4.701945181048143 4.701945181048143 -177.2091619818234 187.00978766741036 1.322201366963372 0.0 0.5021464029637546 0.0 0.0 0.0 4.701945181048143 -3.130180320024094e-10 -0.00879197763498300466943027664188170128 0.00582470293409468688038602773748266245 0.0 0.0 0.0 0.0
";

#[test]
fn zero_external_reflects_only_a_line_nearer_its_corner() {
    on_stack(|| check_reflection_rows(REFLECTION_MARGIN));
}

/// Boxes of the verifier's sweeps where the line passes nearer the corner
/// (1, 0) by more than 1/128 and `zero_external` reflects the chart only where
/// one of its `tfun` selects a pole next to the corner's vertex that the
/// rounding resolves poorly (`box_complex::tfun_vertex_suspicious`), in their
/// own labelling (references: COLLIER 1.2.9 in quadruple precision, COLI and
/// DD agreeing, or COLI and this crate's native 256-bit value). The first
/// four have no such pole and keep the original orientation: reflected (at
/// cb5370d) they lost 4.1e-12 to 1.0e-11 (v1 1804: a far pole of the nearly
/// vanishing `b = p4²`), 3.3e-13 to 7.2e-12 (v1 7832: a line product), 7.9e-13
/// to 3.0e-12 (v1 16613) and 8.1e-13 to 1.3e-12 (g1 10856: a `K` cancelling
/// by 2e3 and 400) on the binary64 routes. The last three are reflected:
/// unreflected (09b2c74) they lose 1.8e-11 to 2.2e-11 (v1 17278: a pole whose
/// rounding the offset from the vertex magnifies 1e4 times), 9.4e-12 to
/// 9.8e-12 (v1 18471: a pole 38 times nearer the vertex than its size) and
/// 3.2e-12 to 4.1e-12 (SW 29392 of the NinjaRS streams' sweep: 20 times).
/// Every binary64 route within its row's tolerance (the last field),
/// DoubleFloat within 1e-14.
const REFLECTION_POLE: &str = "\
4 1.0 0.28656384818680486 42.02584986987939 -9.077630159488835e-13 0.0 1.097330042297487e-13 -0.26115558711995274 13.830744098936528 -0.0009531705655522072 0.0 0.0 42.02584986987939 -3.062623934004004e-07 12.466532305719312 -7.942213616841566e-13 0.00306554602027375734616617387437222126 1.52644889366157907794829027880865826E-7 0.0 0.0 0.0 0.0 2e-12
4 1.0 83.3855670021432 0.0 11.232603227340993 83.38548078042452 -0.9024314250092198 -1.59358277821813 83.3855670020681 0.0 0.03957245932010741 0.0 0.03957245932010741 0.0 0.0 0.0 0.00927510205593330776947712802887193334 0.0845509034072344636972574282900693038 0.0 0.0 0.0 0.0 1e-12
4 1.0 2.948377372752437e-11 28.31952850318703 0.0 0.015949891160530368 -0.17922683605460168 -229.35946438621744 0.0 0.0 28.31952850318703 -3.1837729083046333e-11 72.0747200576466 0.0 72.0747200576466 0.0 0.000132413763514520412429270062876527252 7.19828033227632298322485743265276321E-17 0.0 0.0 0.0 0.0 1e-12
4 1.0 188.6865735725467 0.0 0.7175303256548632 -0.12566377922328004 -2.596072948553744 -0.011750231886380956 15.811414894439205 0.0 7.195596015821104 -0.05897417412953987 9.001577285728034 0.0 0.0 0.0 0.000912544938334547156413868319884028920 0.00484161524092646691201311421080879498 0.0 0.0 0.0 0.0 2e-13
4 1.0 0.0 0.01767368687841736 -3.251557176123042e-13 0.01767368868492391 -0.08618317740844156 260.1220143634163 0.01767368687841736 -1.016157680684762e-06 0.01767368687841736 -1.016157680684762e-06 0.0 0.0 1.9271531991107602 -6.479744509773929e-06 -1.23985213200944342635652881892295716 1.55916141513682332420475476068907249 0.0 0.0 0.0 0.0 1e-12
4 1.0 1.3104711065925037 1.3104711065930064 0.0 0.013927560143652096 2.745048083797347 0.24764122273995193 0.0 0.0 1.3104711065925037 -0.033826877553513365 2.553002723799046 -8.659481203644983e-06 0.013927560143652096 0.0 -30.7952611290646049686034280959006578 22.1744616936444843648885335846340354 0.0 0.0 0.0 0.0 1e-13
4 1.0 1467.1852045767862 0.243012274094481 -0.13365167391458466 13.512190149234277 -49.81300097992969 0.0 12.363261727598141 0.0 0.243012274094481 0.0 0.0 0.0 0.243012274094481 0.0 -0.0116948527312188023630420582502082495 0.0396641692032404853740058663027177328 0.0 0.0 0.0 0.0 1e-12
";

#[test]
fn zero_external_reflects_only_a_chart_with_a_pole_at_its_corner() {
    on_stack(|| check_reflection_rows(REFLECTION_POLE));
}

/// Every binary64 route (hand-written and generated native code, the
/// interpreted expressions and SymJIT) within the row's tolerance (1e-12 where
/// none is given) and DoubleFloat within 1e-14, in the rows' own labelling.
fn check_reflection_rows(data: &'static str) {
    let mut native: Vec<_> = [true, false]
        .map(|handwritten| {
            let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
            evaluator.set_handwritten(handwritten);
            (handwritten, evaluator)
        })
        .into();
    let mut double = NativeEvaluator::<DoubleFloat>::new(ScalarIntegral::D0).unwrap();
    let mut failures = Vec::new();
    for (line, (args, expected, tolerance)) in rows(data).iter().enumerate() {
        let binary64 = tolerance.unwrap_or(1e-12);
        let mut values = Vec::new();
        for (handwritten, evaluator) in &mut native {
            let mut output = [C::new(f64::NAN, f64::NAN); 3];
            evaluator.evaluate(args, &mut output).unwrap();
            values.push((
                format!("native (hand-written: {handwritten})"),
                output[0],
                binary64,
            ));
        }
        for backend in [EvaluationBackend::Expression, EvaluationBackend::SymJit] {
            let mut output = [C::new(f64::NAN, f64::NAN); 3];
            evaluate_with_backend(ScalarIntegral::D0, args, &mut output, backend).unwrap();
            values.push((backend.name().to_string(), output[0], binary64));
        }
        let input: Vec<_> = args
            .iter()
            .map(|z| Complex::new(DoubleFloat::from(z.re), DoubleFloat::from(z.im)))
            .collect();
        let zero = Complex::new(DoubleFloat::from(0.), DoubleFloat::from(0.));
        let mut output = [zero; 3];
        double.evaluate(&input, &mut output).unwrap();
        let value = C::new(output[0].re.to_f64(), output[0].im.to_f64());
        values.push(("DoubleFloat".to_string(), value, 1e-14));
        for (route, value, tolerance) in values {
            let error = relative(value, *expected);
            if !(error <= tolerance) {
                failures.push(format!(
                    "row {} {route}: {value:?} vs {expected:?} ({error:e})",
                    line + 1
                ));
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

/// A box binary64 still gets wrong (reference: COLLIER 1.2.9 in quadruple
/// precision, COLI branch): three real masses, `p1²` one ulp below
/// `m2² = m4² = p3²`, the massless `m1` between the legs `p1` and `p4`, and
/// `s23 = 0`: in the labellings the contour box takes, the conic's value and
/// gradient vanish at the massless corner to within rounding. Every binary64
/// route was wrong in every relabelling at 09b2c74 (NaN in 8, by 1.4e-2 in 4
/// and by 0.3 in 12); now NaN in 16 and 1.4e-2 in 8. (The box with a Gram
/// determinant at the rounding level that was a known failure here, 5.1e-3,
/// is row 3 of tests/data/contour_gram_rounding.txt.)
const KNOWN_FAILURES: &str = "\
4 1.0 0.1113435943976259 108.23851883958666 0.11134359439762591 -0.012515862621692797 9.0 0.0 0.0 0.0 0.11134359439762591 0.0 9.0 0.0 0.11134359439762591 0.0 -2.77080305527435389896e+1 2.11066704520210922998e+1 0.0 0.0 0.0 0.0
";

#[test]
#[ignore = "known binary64 failure: a conic singular at a massless corner"]
fn contour_box_known_failures_match_collier_in_binary64() {
    on_stack(|| check_binary64(KNOWN_FAILURES, 1, BINARY64, PRECISE, false));
}
