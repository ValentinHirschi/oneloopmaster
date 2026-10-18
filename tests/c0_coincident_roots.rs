//! A two-mass C0 whose trif2 roots nearly coincide (NinjaRS syn-n6-r6, the
//! triangle behind lane 0's error): p² = (-237088.34352377, -236743.32762192388,
//! 0), m² = (29929 - 259.5 i, 0, 6464.16 - 168.84 i), mu² = 10^4, with
//! sqrt(lambda)/|p²| = 1.5e-3. The divided dilogarithm differences of the two
//! roots cancelled by that much, and those of the reflected labelling, whose
//! second root is zero, cancelled against Li2(1) = pi²/6. Binary64 was off by
//! 2.4e-13 (native, canonical labelling), 2.3e-14 (native, reflected) and
//! 3.8e-13 (SymJIT); OneLOop 3.6 by 3.3e-15. The term-wise divided difference
//! (`sheet_exact::near_difference`) and logarithms next to one formed from
//! their offsets leave 2e-15 to 9e-15 on every binary64 route. Reference:
//! 256 and 400 bits, which agree with OneLOop 3.6 in quadruple precision.
#![cfg(feature = "generated-evaluators")]
use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend,
};
use symbolica::domains::float::{Complex, Float, RealLike};

const MOMENTA: [f64; 3] = [-237088.34352377, -236743.32762192388, 0.0];
const MASSES: [(f64, f64); 3] = [
    (29929.0, -259.5),
    (0.0, 0.0),
    (6464.160000000001, -168.84000000000003),
];
const EXPECTED: [f64; 2] = [-1.148_249_753_268_316_364_585_722e-5, -5.243_707_138_558_474_789_137e-8];
const EXPECTED_DIGITS: [&str; 2] = ["-1.148249753268316364585722e-5", "-5.243707138558474789137e-8"];

/// Relative tolerances in binary64 and at 256 bits.
const BINARY64: f64 = 1.5e-14;
const PRECISE: f64 = 1e-24;

/// The six labellings: cyclic rotations and their reflections.
fn images() -> Vec<(String, [Complex<f64>; 7])> {
    let mut images = Vec::new();
    for rotation in 0..3 {
        for reflected in [false, true] {
            let mut p: [f64; 3] = core::array::from_fn(|i| MOMENTA[(i + rotation) % 3]);
            let mut m: [(f64, f64); 3] = core::array::from_fn(|i| MASSES[(i + rotation) % 3]);
            if reflected {
                p = [p[2], p[1], p[0]];
                m = [m[0], m[2], m[1]];
            }
            let label = format!("rotation {rotation}{}", if reflected { ", reflected" } else { "" });
            let input = [
                Complex::new(p[0], 0.0),
                Complex::new(p[1], 0.0),
                Complex::new(p[2], 0.0),
                Complex::new(m[0].0, m[0].1),
                Complex::new(m[1].0, m[1].1),
                Complex::new(m[2].0, m[2].1),
                Complex::new(1e4, 0.0),
            ];
            images.push((label, input));
        }
    }
    images
}

fn error(value: Complex<f64>) -> f64 {
    (value.re - EXPECTED[0]).hypot(value.im - EXPECTED[1]) / EXPECTED[0].hypot(EXPECTED[1])
}

fn precise_error(value: &Complex<Float>) -> f64 {
    let bits = value.re.prec();
    let [re, im] = EXPECTED_DIGITS.map(|digits| Float::parse(digits, Some(bits)).unwrap());
    let difference = (value.re.clone() - &re).to_f64().hypot((value.im.clone() - &im).to_f64());
    difference / EXPECTED[0].hypot(EXPECTED[1])
}

/// Every route: native binary64 (with and without the hand-written kernels,
/// which leave two-mass triangles to the generated code), the interpreted
/// expressions and SymJIT in binary64, and native code and interpreted
/// expressions at 256 bits.
#[test]
fn nearly_coincident_triangle_roots_match_oneloop_on_every_route() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(check_every_route)
        .unwrap()
        .join()
        .unwrap();
}

fn check_every_route() {
    let mut failures = Vec::new();
    let mut check = |route: String, error: f64, tolerance: f64| {
        if !(error <= tolerance) {
            failures.push(format!("{route}: relative error {error:e}"));
        }
    };
    let bits = 256;
    let fresh =
        || core::array::from_fn::<_, 3, _>(|_| Complex::new(Float::new(bits), Float::new(bits)));
    let mut native_256 =
        NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::C0, bits).unwrap();
    let mut expression_256 = PrecisionEvaluator::with_binary_precision_and_backend(
        ScalarIntegral::C0,
        bits,
        EvaluationBackend::Expression,
    )
    .unwrap();
    for (label, input) in images() {
        for handwritten in [true, false] {
            let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::C0).unwrap();
            evaluator.set_handwritten(handwritten);
            let mut output = [Complex::new(0.0, 0.0); 3];
            evaluator.evaluate(&input, &mut output).unwrap();
            check(
                format!("{label}: native (hand-written: {handwritten})"),
                error(output[0]),
                BINARY64,
            );
        }
        for backend in [EvaluationBackend::Expression, EvaluationBackend::SymJit] {
            let mut output = [Complex::new(0.0, 0.0); 3];
            evaluate_with_backend(ScalarIntegral::C0, &input, &mut output, backend).unwrap();
            check(format!("{label}: {}", backend.name()), error(output[0]), BINARY64);
        }
        let input = input.map(|z| Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im)));
        let mut output = fresh();
        native_256.evaluate(&input, &mut output).unwrap();
        check(format!("{label}: native 256 bits"), precise_error(&output[0]), PRECISE);
        let mut output = fresh();
        expression_256.evaluate(&input, &mut output).unwrap();
        check(format!("{label}: expression 256 bits"), precise_error(&output[0]), PRECISE);
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
