//! The contour box next to a pseudo-threshold beside a massless line.
//!
//! Each quadratic edge polynomial of the contour box (`box_complex::s3`) is the
//! Feynman polynomial of a two-point subgraph, `p x² + (mi - mj - p) x + mj`,
//! whose discriminant is the Källén function `λ(p, mi, mj)` and whose values
//! at the ends are the masses. Formed from the rounded coefficients, the
//! discriminant kept a rounding residue at a pseudo-threshold next to a
//! massless line (`p = mi`, `mj = 0`, where `λ = (p - mi)²`), whose square
//! root moved the roots by sqrt(eps): binary64 lost 8 to 10 digits there.
//! Taking the double-root forms whenever the discriminant was within rounding
//! of zero restored the threshold but cost the roots their separation next to
//! it (0.56 r relative error at p = mi (1 + r), 2.8e-7 at r = 5e-7). The
//! callers now pass `λ` formed from the invariants as
//! `(p - (√mi + √mj)²)(p - (√mi - √mj)²)` and the mass at the far end, so the
//! threshold is decided exactly and a massless end is an exact root.
//!
//! References: native code in `Float` at 192 bits on the same binary64 inputs.
#![cfg(feature = "generated-evaluators")]
use oneloop::{EvaluationBackend, NativeEvaluator, ScalarIntegral, evaluate_with_backend};
use symbolica::domains::float::{Complex, Float, RealLike};

type C = Complex<f64>;

const BITS: u32 = 192;

fn on_stack(test: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(test)
        .unwrap()
        .join()
        .unwrap();
}

/// `kind mu² p.. (re im)..` rows as evaluator inputs (mu² last).
fn input(row: &str) -> Vec<C> {
    let v: Vec<f64> = row.split_whitespace().map(|x| x.parse().unwrap()).collect();
    assert_eq!(v.len(), 16);
    let mut args: Vec<C> = v[2..8].iter().map(|&x| C::new(x, 0.0)).collect();
    args.extend(v[8..].chunks_exact(2).map(|m| C::new(m[0], m[1])));
    args.push(C::new(v[1], 0.0));
    args
}

/// The worst relative error of the finite part over every binary64 route.
fn binary64_errors(inputs: &[Vec<C>]) -> Vec<String> {
    let mut precise = NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::D0, BITS)
        .unwrap();
    let mut native = [true, false].map(|handwritten| {
        let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
        evaluator.set_handwritten(handwritten);
        (handwritten, evaluator)
    });
    let mut errors = Vec::new();
    for args in inputs {
        let high = args
            .iter()
            .map(|z| Complex::new(Float::with_val(BITS, z.re), Float::with_val(BITS, z.im)))
            .collect::<Vec<_>>();
        let mut out =
            core::array::from_fn::<_, 3, _>(|_| Complex::new(Float::new(BITS), Float::new(BITS)));
        precise.evaluate(&high, &mut out).unwrap();
        let reference = C::new(out[0].re.to_f64(), out[0].im.to_f64());
        let mut routes = Vec::new();
        for (handwritten, evaluator) in &mut native {
            let mut output = [C::new(f64::NAN, f64::NAN); 3];
            evaluator.evaluate(args, &mut output).unwrap();
            routes.push((format!("native (hand-written: {handwritten})"), output[0]));
        }
        for backend in [EvaluationBackend::Expression, EvaluationBackend::SymJit] {
            let mut output = [C::new(f64::NAN, f64::NAN); 3];
            evaluate_with_backend(ScalarIntegral::D0, args, &mut output, backend).unwrap();
            routes.push((backend.name().to_string(), output[0]));
        }
        for (route, value) in routes {
            let difference = value - reference;
            let error = difference.re.hypot(difference.im) / reference.re.hypot(reference.im);
            errors.push(format!("{error:.1e} {route} p0={:?}: {value:?} vs {reference:?}", args[0].re));
        }
    }
    errors
}

fn check(inputs: &[Vec<C>], tolerance: f64) {
    let errors = binary64_errors(inputs);
    let failures = errors
        .iter()
        .filter(|line| !(line.split(' ').next().unwrap().parse::<f64>().unwrap() <= tolerance))
        .collect::<Vec<_>>();
    assert!(failures.is_empty(), "{failures:#?}");
}

/// p² = (M (1 + r), 81738, -181550, -222558), s12 = -43301, s23 = -155866,
/// m² = (M, 0, 100, 100), mu² = M = 8315.178393760001: p1² next to m1² beside
/// the massless second line, r from 1e-15 to 1e-5 on both sides and r = 0.
/// Every binary64 route is within 1e-12 (it was up to 2.8e-7 at r = 5e-7,
/// 1e-10 to 1e-9 for r beyond 1e-6, and 1e-14 at the threshold itself).
#[test]
fn pseudo_threshold_neighbourhood_in_binary64() {
    on_stack(|| {
        const M: f64 = 8315.178393760001;
        let rest = "81738 -181550 -222558 -43301 -155866 8315.178393760001 0 0 0 100 0 100 0";
        let mut inputs = vec![input(&format!("4 {M} {M} {rest}"))];
        for exponent in -15..=-5 {
            for sign in [1.0, -1.0] {
                let p = M * (1.0 + sign * 10f64.powi(exponent));
                inputs.push(input(&format!("4 {M} {p} {rest}")));
            }
        }
        check(&inputs, 1e-12);
    });
}

/// Sweep boxes (`tests/oracle_sweep.py`) with an invariant one ulp from the
/// pseudo-threshold beside a massless line, one with real masses (s12 next to
/// m1²) and one with a complex mass (p2² next to m3²). The second root then
/// rounds onto the end of the integration (or next to it, on either side),
/// and the rounded coefficients no longer give the sign of its slope: with
/// the roots separated (not collapsed within rounding), binary64 was O(1)
/// wrong until both roots continued from the exact end values. Every binary64
/// route within 1e-12.
#[test]
fn within_an_ulp_of_the_pseudo_threshold_in_binary64() {
    on_stack(|| {
        let rows = [
            "4 1.0 0.03588045152205349 0.0 -0.21764653930803976 4.3806008259221025 0.9053892696515414 0.09607948475970099 0.9053892696515415 0.0 7.306691271758426 0.0 0.0 0.0 0.9053892696515415 0.0",
            "4 1.0 -0.07941272858779677 1.4043727320632262 -28.59543039506344 9.0 9.0 0.0 9.0 -0.01438797842123132 0.0 0.0 1.4043727320632264 0.0 2.25 0.0",
        ];
        check(&rows.map(input), 1e-12);
    });
}

/// p² = (6.686564391522424, 4, 0, M), s12 = 0, s23 = 4, m² = (M, 0.1067...,
/// 0.5221..., M), M = 0.6479...: relabelled to the zero-external chart with
/// p4² = s23, where the contour line's coefficient g + h vanishes. The rounded
/// g + h (1e-16) put the root of one edge at 1e15, which one of the two
/// equivalent contours (complex conjugate betas, chosen by a tie in their
/// distances) evaluated O(1) wrong in binary64; the line's sums are now exact.
#[test]
fn vanishing_line_coefficient_in_binary64() {
    on_stack(|| {
        let row = "4 0.34422897492297583 6.686564391522424 4.0 0.0 0.6479305255028737 0.0 4.0 0.6479305255028737 0.0 0.10671797382773325 0.0 0.5221278889743961 0.0 0.6479305255028737 0.0";
        check(&[input(row)], 1e-12);
    });
}

/// A box with a vanishing channel (s23 = 0) in its four cyclic labellings,
/// two of which take the contour box's zero-external chart: p4² is within
/// two ulps of m1² beside the massless fourth line. Edge quadratics that end
/// on that line have the root 1 there, which binary64 computed as 1 + 9e-15
/// and which then flipped endpoint decisions (2.2e-10 to 2.4e-10 before).
#[test]
fn massless_end_roots_in_binary64() {
    on_stack(|| {
        let rows = [
            "4 1.0 -150.13830018861105 28.715479639249917 5.379501694598977 3.1562378331595013 -424.76251272146203 0.0 3.156237833159502 0.0 3.4180033385407063 0.0 3.4180033385407063 0.0 0.0 0.0",
            "4 1.0 -424.76251272146203 28.715479639249917 0.0 3.1562378331595013 -150.13830018861105 5.379501694598977 3.156237833159502 0.0 3.4180033385407063 0.0 3.4180033385407063 0.0 0.0 0.0",
            "4 1.0 3.1562378331595013 0.0 28.715479639249917 -424.76251272146203 -150.13830018861105 5.379501694598977 3.156237833159502 0.0 0.0 0.0 3.4180033385407063 0.0 3.4180033385407063 0.0",
            "4 1.0 3.1562378331595013 5.379501694598977 28.715479639249917 -150.13830018861105 -424.76251272146203 0.0 3.156237833159502 0.0 0.0 0.0 3.4180033385407063 0.0 3.4180033385407063 0.0",
        ];
        check(&rows.map(input), 1e-11);
    });
}
