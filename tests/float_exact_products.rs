//! C0 and D0 rows whose `Float` values were wrong while binary64 was right
//! (tests/data/float_exact_products.txt): on the native route and in the
//! interpreted expressions, at precisions that varied from row to row, by
//! 5e-2 to 4e6 of the value in the D0 oracle sweeps (seeds 1-4: 62 native and
//! 44 interpreted rows at 256 bits) and in the C0 sweeps (46 and 33).
//!
//! numerica rounded the more precise factor of a `Float` product to the
//! other's precision before multiplying, so `a b` and `b a` could differ in
//! the last bit and `z conj(z)` kept an imaginary residue; and a factor `-1`,
//! with which evaluators form `a - b`, rounded `b` to the working precision
//! although square roots carry one more bit, so `x - |x|` (the sign test) and
//! `z - conj(z)` (the imaginary part) kept one-bit residues. Exact sign and
//! zero tests branched on them: the order of the roots of a quadratic, the
//! lip of a logarithm, the sheet of a dilogarithm. Symbolica's numerica now
//! rounds products and quotients once and keeps factors of exactly `±1`
//! exact.
//!
//! Every route is checked in all propagator relabellings: binary64 (native
//! with and without the hand-written kernels, interpreted, SymJIT),
//! DoubleFloat, and the native and interpreted `Float` routes at 128 to 512
//! bits, which must match the references and agree with each other.
#![cfg(feature = "generated-evaluators")]
use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend,
};
use symbolica::domains::float::{Complex, DoubleFloat, Float, FloatLike, RealLike};

type C = Complex<f64>;

/// Binary64 against the references.
const TOLERANCE: f64 = 1e-11;
/// DoubleFloat against the references.
const DOUBLE: f64 = 1e-15;
/// Working precisions of the `Float` routes. The residues fell on different
/// rows at different precisions.
const BITS: [u32; 5] = [128, 192, 256, 320, 512];
/// Agreement of the relabellings with each other, unless the values track
/// fewer bits (some rows lose up to 70 bits to cancellations).
const AGREEMENT: f64 = 1e-30;
/// Tracked bits of the dominant component: the wrong values tracked 1 to 20.
const TRACKED: u32 = 48;

/// The relative uncertainty that `bits` tracked bits claim, with 8 bits of
/// slack for numerica's estimate.
fn claimed(bits: u32) -> f64 {
    2f64.powi(8 - bits as i32)
}

struct Row {
    family: ScalarIntegral,
    args: Vec<C>,
    reference: (String, String),
    tolerance: f64,
    /// Only the cyclic relabellings (see the data file).
    rotations_only: bool,
    /// The interpreted `Float` route is an open defect here (STATUS.md),
    /// checked by the ignored [open_interpreted_float_defect] instead.
    open_interpreted: bool,
}

fn rows() -> Vec<Row> {
    include_str!("data/float_exact_products.txt")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| {
            let tokens: Vec<&str> = line.split_whitespace().collect();
            let (values, reference) = tokens.split_at(tokens.len() - 4);
            let v: Vec<f64> = values.iter().map(|x| x.parse().unwrap()).collect();
            let (family, momenta) = match v[0] as i64 {
                3 => (ScalarIntegral::C0, 3),
                4 => (ScalarIntegral::D0, 6),
                kind => panic!("unexpected kind {kind}"),
            };
            let masses = if momenta == 3 { 3 } else { 4 };
            assert_eq!(v.len(), 2 + momenta + 2 * masses, "malformed row {line}");
            let mut args: Vec<C> = v[2..2 + momenta].iter().map(|&x| C::new(x, 0.)).collect();
            args.extend(v[2 + momenta..].chunks_exact(2).map(|z| C::new(z[0], z[1])));
            args.push(C::new(v[1], 0.));
            Row {
                family,
                args,
                reference: (reference[0].to_string(), reference[1].to_string()),
                tolerance: reference[2].parse().unwrap(),
                rotations_only: reference[3].split(',').any(|flag| match flag {
                    "all" | "open-interpreted" => false,
                    "rotations" => true,
                    other => panic!("unexpected flag {other}"),
                }),
                open_interpreted: reference[3]
                    .split(',')
                    .any(|flag| flag == "open-interpreted"),
            }
        })
        .collect()
}

/// A triangle in its 6 or a box in its 24 propagator relabellings.
fn images(row: &Row) -> Vec<Vec<C>> {
    let args = &row.args;
    let (n, masses) = match row.family {
        ScalarIntegral::C0 => (3, 3),
        _ => (4, 6),
    };
    // The invariant between propagators i and j.
    let invariant = |i: usize, j: usize| {
        let index = if n == 3 {
            match (i.min(j), i.max(j)) {
                (0, 1) => 0,
                (1, 2) => 1,
                _ => 2,
            }
        } else {
            match (i.min(j), i.max(j)) {
                (0, 1) => 0,
                (1, 2) => 1,
                (2, 3) => 2,
                (0, 3) => 3,
                (0, 2) => 4,
                _ => 5,
            }
        };
        args[index]
    };
    let mut images = Vec::new();
    let mut order: Vec<usize> = (0..n).collect();
    permutations(&mut order, 0, &mut |order| {
        let cyclic = (0..n).all(|k| order[(k + 1) % n] == (order[k] + 1) % n);
        if row.rotations_only && !cyclic {
            return;
        }
        let mut image: Vec<C> = (0..n)
            .map(|k| invariant(order[k], order[(k + 1) % n]))
            .collect();
        if n == 4 {
            image.push(invariant(order[0], order[2]));
            image.push(invariant(order[1], order[3]));
        }
        image.extend(order.iter().map(|&k| args[masses + k]));
        image.push(args[masses + n]);
        images.push(image);
    });
    let count = match (row.rotations_only, n) {
        (true, _) => n,
        (false, 3) => 6,
        _ => 24,
    };
    assert_eq!(images.len(), count);
    images
}

fn permutations(order: &mut Vec<usize>, start: usize, visit: &mut impl FnMut(&[usize])) {
    if start == order.len() {
        visit(order);
        return;
    }
    for i in start..order.len() {
        order.swap(start, i);
        permutations(order, start + 1, visit);
        order.swap(start, i);
    }
}

fn reference_f64(row: &Row) -> C {
    C::new(
        row.reference.0.parse().unwrap(),
        row.reference.1.parse().unwrap(),
    )
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

fn report(failures: Vec<String>) {
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// The binary64 routes and DoubleFloat, which were right before the fix and
/// must stay so.
#[test]
fn binary64_and_double_float_routes_match_the_references() {
    on_stack(|| {
        let mut failures = Vec::new();
        for (line, row) in rows().iter().enumerate() {
            let expected = reference_f64(row);
            let mut native: Vec<_> = [true, false]
                .map(|handwritten| {
                    let mut evaluator = NativeEvaluator::<f64>::new(row.family).unwrap();
                    evaluator.set_handwritten(handwritten);
                    (handwritten, evaluator)
                })
                .into();
            let mut double = NativeEvaluator::<DoubleFloat>::new(row.family).unwrap();
            for (number, image) in images(row).into_iter().enumerate() {
                let mut check = |route: String, value: C, tolerance: f64| {
                    let error = relative(value, expected);
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
                    evaluate_with_backend(row.family, &image, &mut output, backend).unwrap();
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
                    DOUBLE.max(row.tolerance),
                );
            }
        }
        report(failures);
    });
}

/// The native and interpreted `Float` routes at 128 to 512 bits, in every
/// relabelling: at least 48 bits are tracked, and within the precision they
/// track (or the row's tolerance: 1e-32 for the 40-digit references, 1e-14
/// for COLLIER's quadruple precision) the finite part matches the reference
/// and the first relabelling (to 1e-30 otherwise); the poles vanish.
#[test]
fn float_routes_match_the_references_at_every_precision() {
    on_stack(|| report(check_float_routes(&BITS, false)));
}

/// The interpreted `Float` route on the rows where it is still an open
/// defect: the two-mass C0 at vanishing invariants keeps 14 to 19 tracked
/// bits at 112, 128 and 224 bits (wrong by 5.3e-3 to 1.2e-2), where the
/// native route keeps them all (STATUS.md, open limitations).
#[test]
#[ignore = "open: interpreted Float C0 at vanishing invariants (STATUS.md)"]
fn open_interpreted_float_defect() {
    on_stack(|| report(check_float_routes(&[112, 128, 224], true)));
}

/// Checks of [float_routes_match_the_references_at_every_precision]; with
/// `open` only the interpreted route of the rows with an open defect there.
fn check_float_routes(precisions: &[u32], open: bool) -> Vec<String> {
    let mut failures = Vec::new();
    let rows = rows();
    for &bits in precisions {
        for family in [ScalarIntegral::C0, ScalarIntegral::D0] {
            let mut native = NativeEvaluator::<Float>::with_binary_precision(family, bits).unwrap();
            let mut expression = PrecisionEvaluator::with_binary_precision_and_backend(
                family,
                bits,
                EvaluationBackend::Expression,
            )
            .unwrap();
            for (line, row) in rows.iter().enumerate().filter(|(_, r)| r.family == family) {
                let reference = Complex::new(
                    Float::parse(&row.reference.0, Some(bits + 64)).unwrap(),
                    Float::parse(&row.reference.1, Some(bits + 64)).unwrap(),
                );
                let scale = reference.re.to_f64().hypot(reference.im.to_f64());
                for route in ["native", "expression"] {
                    if open != (route == "expression" && row.open_interpreted) {
                        continue;
                    }
                    let mut first: Option<Complex<Float>> = None;
                    for (number, image) in images(row).into_iter().enumerate() {
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
                        let value = output[0].clone();
                        let distance = |other: &Complex<Float>| {
                            (value.re.clone() - &other.re)
                                .to_f64()
                                .hypot((value.im.clone() - &other.im).to_f64())
                                / scale
                        };
                        // The precision of the dominant component (the other can
                        // be a rounding residue of an analytic zero).
                        let tracked_bits = |z: &Complex<Float>| {
                            if z.re.to_f64().abs() >= z.im.to_f64().abs() {
                                z.re.get_precision()
                            } else {
                                z.im.get_precision()
                            }
                        };
                        let tracked = tracked_bits(&value);
                        let error = distance(&reference);
                        let spread = first.as_ref().map_or(0., distance);
                        let first_tracked = first.as_ref().map_or(tracked, tracked_bits);
                        let poles = output[1..].iter().all(|pole| {
                            pole.re.to_f64().hypot(pole.im.to_f64()) <= AGREEMENT * scale
                        });
                        if !(error <= row.tolerance.max(claimed(tracked)))
                            || !(spread <= AGREEMENT.max(claimed(tracked.min(first_tracked))))
                            || !poles
                            || tracked < TRACKED
                        {
                            failures.push(format!(
                                "row {line} image {number} {route} {bits} bits: {} {} \
                                 (error {error:e}, spread from image 0 {spread:e}, \
                                 {tracked} tracked bits, poles vanish: {poles})",
                                value.re, value.im
                            ));
                        }
                        first.get_or_insert(value);
                    }
                }
            }
        }
    }
    failures
}
