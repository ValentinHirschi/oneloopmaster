//! Regression points for continuation, root-labelling and binary64 defects
//! found by the OneLOop sweep (tests/oracle_sweep.py), checked against
//! OneLOop 3.6 in quadruple precision and under the integrals' symmetries.
//!
//! tests/data/permutation_regressions.txt holds rows that the expressions
//! before this fix evaluated wrongly: B0/dB0 with one massless line (NaN in
//! binary64), two-mass C0 with real roots on the wrong lip, at thresholds and
//! with a zero root out of place, the 't Hooft-Veltman C0 continuation, the
//! one-mass and massless infrared boxes, and the finite one-, three- and
//! four-mass boxes (root lips, chart choice and exact double roots). References
//! are taken where OneLOop's double and quadruple precision results agree,
//! except for the real-mass three-mass boxes of lines 471-474 and 585-588,
//! where both agree on the wrong imaginary part of OneLOop's real-mass boxf33
//! chart (-26.894 and -0.15379): their references are the contour box at 256
//! bits (native code and interpreted expressions agree), which OneLOop's own
//! complex-mass route approaches as the widths vanish
//! (tests/d0_real_mass_orderings.rs).
#![cfg(feature = "generated-evaluators")]
use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend,
};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::{
        float::{Complex, Float, RealLike},
        rational::Rational,
    },
};

type C = Complex<f64>;

struct Row {
    family: ScalarIntegral,
    args: Vec<C>,
    expected: [C; 3],
    line: usize,
}

fn rows() -> Vec<Row> {
    include_str!("data/permutation_regressions.txt")
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.starts_with('#') && !line.trim().is_empty())
        .map(|(index, line)| {
            let values: Vec<f64> = line
                .split_whitespace()
                .map(|v| v.parse().unwrap())
                .collect();
            let (family, momenta, masses) = match values[0] as i64 {
                1 => (ScalarIntegral::A0, 0, 1),
                2 => (ScalarIntegral::B0, 1, 2),
                -2 => (ScalarIntegral::DB0, 1, 2),
                3 => (ScalarIntegral::C0, 3, 3),
                4 => (ScalarIntegral::D0, 6, 4),
                kind => panic!("unknown kind {kind}"),
            };
            let mut args: Vec<C> = values[2..2 + momenta]
                .iter()
                .map(|&v| C::new(v, 0.))
                .collect();
            let masses_end = 2 + momenta + 2 * masses;
            args.extend(
                values[2 + momenta..masses_end]
                    .chunks_exact(2)
                    .map(|v| C::new(v[0], v[1])),
            );
            args.push(C::new(values[1], 0.));
            assert_eq!(values.len(), masses_end + 6, "line {}", index + 1);
            let expected = core::array::from_fn(|i| {
                C::new(values[masses_end + 2 * i], values[masses_end + 2 * i + 1])
            });
            Row {
                family,
                args,
                expected,
                line: index + 1,
            }
        })
        .collect()
}

/// Largest coefficient error relative to the largest reference coefficient.
fn error(actual: &[C], expected: &[C; 3]) -> f64 {
    let scale = expected
        .iter()
        .fold(0.0_f64, |s, z| s.max(z.re.hypot(z.im)));
    let mut worst = 0.0_f64;
    for (a, e) in actual.iter().zip(expected) {
        if !a.re.is_finite() || !a.im.is_finite() {
            return f64::INFINITY;
        }
        worst = worst.max((a.re - e.re).hypot(a.im - e.im));
    }
    worst / scale
}

fn on_stack(action: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(action)
        .unwrap()
        .join()
        .unwrap();
}

/// Binary64 evaluators; `handwritten` selects the hand-written C0/D0 sector
/// kernels (the binary64 default) or the generated code everywhere.
fn family_evaluators(handwritten: bool) -> [NativeEvaluator<f64>; 5] {
    [
        ScalarIntegral::A0,
        ScalarIntegral::B0,
        ScalarIntegral::DB0,
        ScalarIntegral::C0,
        ScalarIntegral::D0,
    ]
    .map(|family| {
        let mut evaluator = NativeEvaluator::<f64>::new(family).unwrap();
        evaluator.set_handwritten(handwritten);
        evaluator
    })
}

fn index(family: ScalarIntegral) -> usize {
    family as usize
}

/// Binary64 tolerance relative to the largest reference coefficient. The
/// worst rows (1.4e-10 to 4.2e-10, depending on the evaluator's operation
/// order) are two-mass triangles with a small complex mass.
const BINARY64_TOLERANCE: f64 = 5e-10;

/// Lines 471-474, a real-mass three-mass box exactly at the pseudo-threshold
/// p3 = m3 next to the massless fourth line, keep nine digits in binary64
/// (5.8e-10 to 7.8e-10 by evaluator; exact at 256 bits) in the contour box,
/// which three-mass boxes now use wherever it supports the kinematics. The
/// boxf33 chart used before had the wrong imaginary part here (-26.894
/// instead of +27.690), as in some relabellings of about a quarter of
/// real-mass three-mass boxes (tests/d0_real_mass_orderings.rs).
const CONTOUR_BINARY64_LIMITED: [usize; 4] = [471, 472, 473, 474];
const CONTOUR_BINARY64_TOLERANCE: f64 = 2e-9;

/// The binary64 tolerance of a fixture line.
fn binary64_tolerance(line: usize, tolerance: f64) -> f64 {
    if CONTOUR_BINARY64_LIMITED.contains(&line) {
        CONTOUR_BINARY64_TOLERANCE
    } else {
        tolerance
    }
}

#[test]
fn native_binary64_matches_oneloop() {
    on_stack(|| {
        let rows = rows();
        assert!(rows.len() > 600);
        let mut failures = Vec::new();
        for handwritten in [true, false] {
            let mut evaluators = family_evaluators(handwritten);
            for row in &rows {
                let mut output = [C::new(f64::NAN, f64::NAN); 3];
                evaluators[index(row.family)]
                    .evaluate(&row.args, &mut output)
                    .unwrap();
                let error = error(&output, &row.expected);
                if !(error <= binary64_tolerance(row.line, BINARY64_TOLERANCE)) {
                    failures.push(format!(
                        "line {} (hand-written kernels: {handwritten}): {output:?} vs {:?} ({error:e})",
                        row.line, row.expected
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

/// Rows whose binary64 values are correct but whose 256-bit native values are
/// not: two-mass triangles at Gram-singular on-shell kinematics (lines 79-81,
/// 109-111, 181-183, 224-226, 247-249, 286-288) or with two vanishing
/// invariants (106-108, 244-246). Symbolica's `Float` tracks precision and
/// rounds the result of a cancellation to the bits it considers significant;
/// at these removable singularities the cancellations are exact, but operands
/// rounded to a few bits spoil the result. The original expressions fail all
/// of them too; which of them fail moves with the evaluator's operation order.
const BINARY64_ONLY: [usize; 24] = [
    79, 80, 81, 106, 107, 108, 109, 110, 111, 181, 182, 183, 224, 225, 226, 244, 245, 246, 247,
    248, 249, 286, 287, 288,
];

/// The same limitation for the interpreted expressions at 256 bits, whose
/// operation order differs from the native code's.
const EXPRESSION_BINARY64_ONLY: [usize; 18] = [
    79, 80, 81, 109, 110, 111, 181, 182, 183, 224, 225, 226, 247, 248, 249, 286, 287, 288,
];

#[test]
fn native_256_bits_matches_oneloop() {
    on_stack(|| {
        let bits = 256;
        let mut failures = Vec::new();
        let mut evaluators: Vec<Option<NativeEvaluator<Float>>> = (0..5).map(|_| None).collect();
        for row in rows()
            .into_iter()
            .filter(|row| !BINARY64_ONLY.contains(&row.line))
        {
            let evaluator = evaluators[index(row.family)].get_or_insert_with(|| {
                NativeEvaluator::<Float>::with_binary_precision(row.family, bits).unwrap()
            });
            let input = row
                .args
                .iter()
                .map(|z| Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im)))
                .collect::<Vec<_>>();
            let mut output = std::array::from_fn::<_, 3, _>(|_| {
                Complex::new(Float::new(bits), Float::new(bits))
            });
            evaluator.evaluate(&input, &mut output).unwrap();
            let output = output.map(|z| C::new(z.re.to_f64(), z.im.to_f64()));
            // The references are quadruple-precision values rounded to binary64.
            let error = error(&output, &row.expected);
            if !(error <= 1e-13) {
                failures.push(format!(
                    "line {}: {output:?} vs {:?} ({error:e})",
                    row.line, row.expected
                ));
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

/// SymJIT compiles the same expressions with its own operation order and
/// packed complex arithmetic; on the most ill-conditioned rows (near-degenerate
/// two-mass triangles and boxes) this loses up to 1.6e-8 more than the
/// interpreter. A wrong branch or sheet would be an O(1) error.
const SYMJIT_TOLERANCE: f64 = 1e-7;

#[test]
fn expression_and_symjit_backends_match_oneloop() {
    on_stack(|| {
        let mut failures = Vec::new();
        for (backend, tolerance) in [
            (EvaluationBackend::Expression, BINARY64_TOLERANCE),
            (EvaluationBackend::SymJit, SYMJIT_TOLERANCE),
        ] {
            for row in rows() {
                let mut output = [C::new(f64::NAN, f64::NAN); 3];
                evaluate_with_backend(row.family, &row.args, &mut output, backend).unwrap();
                let error = error(&output, &row.expected);
                if !(error <= binary64_tolerance(row.line, tolerance)) {
                    failures.push(format!(
                        "line {} ({}): {output:?} vs {:?} ({error:e})",
                        row.line,
                        backend.name(),
                        row.expected
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

/// The interpreted expressions at 256 bits; see [BINARY64_ONLY] for the rows
/// skipped (the interpreter rounds tracked precision like the native code).
#[test]
fn expression_256_bits_matches_oneloop() {
    on_stack(|| {
        let bits = 256;
        let mut failures = Vec::new();
        let mut evaluators: Vec<Option<PrecisionEvaluator>> = (0..5).map(|_| None).collect();
        for row in rows()
            .into_iter()
            .filter(|row| !EXPRESSION_BINARY64_ONLY.contains(&row.line))
        {
            let evaluator = evaluators[index(row.family)].get_or_insert_with(|| {
                PrecisionEvaluator::with_binary_precision_and_backend(
                    row.family,
                    bits,
                    EvaluationBackend::Expression,
                )
                .unwrap()
            });
            let input = row
                .args
                .iter()
                .map(|z| Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im)))
                .collect::<Vec<_>>();
            let mut output = std::array::from_fn::<_, 3, _>(|_| {
                Complex::new(Float::new(bits), Float::new(bits))
            });
            evaluator.evaluate(&input, &mut output).unwrap();
            let output = output.map(|z| C::new(z.re.to_f64(), z.im.to_f64()));
            let error = error(&output, &row.expected);
            if !(error <= 1e-13) {
                failures.push(format!(
                    "line {}: {output:?} vs {:?} ({error:e})",
                    row.line, row.expected
                ));
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

/// Images of a triangle under its rotations and reflections.
fn triangle_images(args: &[C]) -> Vec<Vec<C>> {
    let (p, m, mu) = (&args[0..3], &args[3..6], args[6]);
    let mut images = Vec::new();
    for shift in 0..3 {
        let pr: Vec<C> = (0..3).map(|i| p[(i + shift) % 3]).collect();
        let mr: Vec<C> = (0..3).map(|i| m[(i + shift) % 3]).collect();
        images.push([pr.clone(), mr.clone(), vec![mu]].concat());
        images.push(
            [
                vec![pr[2], pr[1], pr[0]],
                vec![mr[0], mr[2], mr[1]],
                vec![mu],
            ]
            .concat(),
        );
    }
    images
}

/// Images of a box under all 24 propagator relabelings.
fn box_images(args: &[C]) -> Vec<Vec<C>> {
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

/// `(data line, image index)` pairs at which OneLOop itself is unreliable.
fn skipped_images() -> std::collections::HashSet<(usize, usize)> {
    include_str!("data/permutation_image_skips.txt")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| {
            let mut fields = line.split_whitespace().map(|v| v.parse().unwrap());
            (fields.next().unwrap(), fields.next().unwrap())
        })
        .collect()
}

/// Relabelings change the chart (root labels, sectors, channel exchange) and
/// with it the binary64 conditioning; the worst stable image loses 1e-8.
const IMAGE_TOLERANCE: f64 = 2e-8;

#[test]
fn triangles_and_boxes_are_invariant_under_relabelling() {
    on_stack(|| {
        let skipped = skipped_images();
        assert!(skipped.len() < 1000);
        let mut failures = Vec::new();
        let mut checked = 0;
        for handwritten in [true, false] {
            let mut evaluators = family_evaluators(handwritten);
            for row in rows() {
                let images = match row.family {
                    ScalarIntegral::C0 => triangle_images(&row.args),
                    ScalarIntegral::D0 => box_images(&row.args),
                    ScalarIntegral::B0 | ScalarIntegral::DB0 => {
                        let mut swapped = row.args.clone();
                        swapped.swap(1, 2);
                        vec![row.args.clone(), swapped]
                    }
                    ScalarIntegral::A0 => continue,
                };
                for (number, image) in images.into_iter().enumerate() {
                    if skipped.contains(&(row.line, number)) {
                        continue;
                    }
                    checked += 1;
                    let mut output = [C::new(f64::NAN, f64::NAN); 3];
                    evaluators[index(row.family)]
                        .evaluate(&image, &mut output)
                        .unwrap();
                    let error = error(&output, &row.expected);
                    if !(error <= binary64_tolerance(row.line, IMAGE_TOLERANCE)) {
                        failures.push(format!(
                        "line {} image {number} (hand-written kernels: {handwritten}) {image:?}: {output:?} ({error:e})",
                        row.line
                    ));
                    }
                }
            }
        }
        assert!(checked > 15000, "{checked} images checked");
        assert!(
            failures.is_empty(),
            "{} failures:\n{}",
            failures.len(),
            failures.join("\n")
        );
    });
}

/// A nonzero second mass far below the first is as hard as a massless line
/// for the quadratic formula; the roots must keep their digits.
#[test]
fn bubbles_with_a_tiny_second_mass_stay_finite() {
    on_stack(|| {
        let mut b0 = NativeEvaluator::<f64>::new(ScalarIntegral::B0).unwrap();
        let mut db0 = NativeEvaluator::<f64>::new(ScalarIntegral::DB0).unwrap();
        for (p, m0) in [
            (0.6309573444801934, 0.1),
            (-0.6309573444801934, 1.),
            (0.15848931924611134, 1.),
        ] {
            let mut massless = [C::new(0., 0.); 3];
            let mut tiny = [C::new(0., 0.); 3];
            let args = |m1: f64| {
                [
                    C::new(p, 0.),
                    C::new(m0, 0.),
                    C::new(m1, 0.),
                    C::new(1., 0.),
                ]
            };
            for (evaluator, tolerance) in [(&mut b0, 1e-12), (&mut db0, 1e-9)] {
                evaluator.evaluate(&args(0.), &mut massless).unwrap();
                evaluator.evaluate(&args(1e-20), &mut tiny).unwrap();
                for (a, b) in tiny.iter().zip(&massless) {
                    assert!(
                        a.re.is_finite() && a.im.is_finite(),
                        "p={p}, m0={m0}: {tiny:?}"
                    );
                    assert!(
                        (a.re - b.re).hypot(a.im - b.im) <= tolerance * (1. + b.re.hypot(b.im))
                    );
                }
            }
        }
    });
}

/// Zero-momentum tensor bubbles with real masses of opposite sign: the
/// anchor logarithm sits on its cut and must take OneLOop's lip. The finite
/// parts are (-1)^n ∫ x^n [-log(2(1-x) - x - i0)] dx for m0² = 2, m1² = -1.
#[test]
fn zero_momentum_tensor_bubbles_with_opposite_sign_masses() {
    let expected = [
        (-0.429_300_626_542_234_4, -0.872_664_625_997_164_8),
        (0.301_911_389_574_326_4, 0.736_916_795_286_494_7),
        (-0.213_455_694_787_163_2, -0.630_257_785_442_396_8),
        (0.153_843_037_219_820_4, 0.545_576_995_808_597_8),
    ];
    let [p, m0, m1, mu] = ["perm_reg_p", "perm_reg_m0", "perm_reg_m1", "perm_reg_mu"]
        .map(|name| symbolica::symbol!(name).to_atom());
    let builders: [fn(&Atom, &Atom, &Atom, &Atom) -> oneloop::LaurentSeries; 4] =
        [oneloop::b1, oneloop::b11, oneloop::b111, oneloop::b1111];
    for (build, (re, im)) in builders.into_iter().zip(expected) {
        let finite = build(&p, &m0, &m1, &mu).into_coefficients()[0].clone();
        let mut evaluator = Atom::evaluator_multiple(
            &[finite.as_view()],
            &[p.clone(), m0.clone(), m1.clone(), mu.clone()],
        )
        .build()
        .unwrap()
        .map_coeff(&|v: &Complex<Rational>| C::new(v.re.to_f64(), v.im.to_f64()));
        let mut output = [C::new(f64::NAN, f64::NAN)];
        evaluator.evaluate(
            &[
                C::new(0., 0.),
                C::new(2., 0.),
                C::new(-1., 0.),
                C::new(1., 0.),
            ],
            &mut output,
        );
        assert!(
            (output[0].re - re).abs() < 1e-12 && (output[0].im - im).abs() < 1e-12,
            "{output:?} vs ({re}, {im})"
        );
    }
}

/// D0 at exactly real masses in the contour-box regions, against the limit of
/// vanishing widths (tests/data/real_mass_limits.txt). The residue logarithms
/// of the contour box are exactly real there, and their lips must follow the
/// widths; before this was made exact, roundoff chose them and moved results
/// by multiples of 2 pi i / k. OneLOop 3.6 is itself wrong at some of these
/// points (its real-mass value differs from its own vanishing-width limit).
#[test]
fn real_mass_contour_boxes_match_the_vanishing_width_limit() {
    on_stack(|| {
        let rows: Vec<(Vec<C>, [C; 3])> = include_str!("data/real_mass_limits.txt")
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
                (
                    args,
                    core::array::from_fn(|i| C::new(v[16 + 2 * i], v[17 + 2 * i])),
                )
            })
            .collect();
        assert!(rows.len() >= 40);
        let mut failures = Vec::new();
        let mut check = |route: &str, args: &[C], output: [C; 3], expected: &[C; 3]| {
            let error = error(&output, expected);
            if !(error <= 1e-9) {
                failures.push(format!(
                    "{route}: {args:?}: {output:?} vs {expected:?} ({error:e})"
                ));
            }
        };
        for handwritten in [true, false] {
            let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
            evaluator.set_handwritten(handwritten);
            for (args, expected) in &rows {
                let mut output = [C::new(f64::NAN, f64::NAN); 3];
                evaluator.evaluate(args, &mut output).unwrap();
                check(
                    &format!("native (hand-written: {handwritten})"),
                    args,
                    output,
                    expected,
                );
            }
        }
        for backend in [EvaluationBackend::Expression, EvaluationBackend::SymJit] {
            for (args, expected) in &rows {
                let mut output = [C::new(f64::NAN, f64::NAN); 3];
                evaluate_with_backend(ScalarIntegral::D0, args, &mut output, backend).unwrap();
                check(backend.name(), args, output, expected);
            }
        }
        let bits = 256;
        let mut native =
            NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::D0, bits).unwrap();
        let mut expression = PrecisionEvaluator::with_binary_precision_and_backend(
            ScalarIntegral::D0,
            bits,
            EvaluationBackend::Expression,
        )
        .unwrap();
        for (args, expected) in &rows {
            let input = args
                .iter()
                .map(|z| Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im)))
                .collect::<Vec<_>>();
            let fresh = || {
                std::array::from_fn::<_, 3, _>(|_| Complex::new(Float::new(bits), Float::new(bits)))
            };
            let mut output = fresh();
            native.evaluate(&input, &mut output).unwrap();
            check(
                "native 256 bits",
                args,
                output.map(|z| C::new(z.re.to_f64(), z.im.to_f64())),
                expected,
            );
            let mut output = fresh();
            expression.evaluate(&input, &mut output).unwrap();
            check(
                "expression 256 bits",
                args,
                output.map(|z| C::new(z.re.to_f64(), z.im.to_f64())),
                expected,
            );
        }
        assert!(
            failures.is_empty(),
            "{} failures:\n{}",
            failures.len(),
            failures.join("\n")
        );
    });
}
