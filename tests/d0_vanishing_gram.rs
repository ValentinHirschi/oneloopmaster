//! A finite four-mass D0 with equal real masses near a vanishing Gram
//! determinant (NinjaRS gg -> HHH, p² = (202500, 15625, 15625, 15625)), which
//! the contour representation (boxc) evaluates. Its residue term `r0` takes
//! `log(q)/(1 - q)` of quotients of nearly coincident poles lying on opposite
//! lips of a cut: the sheet winding contributes `i pi n/(1 - q)`, large and
//! cancelling against the other terms, and `1 - q` formed from the rounded
//! quotient cost binary64 five digits (6.8e-10 relative; OneLOop 7.9e-15). The
//! winding is now added as `i pi n/(y1 - y2)` with the difference of the poles,
//! leaving 4e-12 to 1.4e-11. That remainder belongs to the orientation the
//! canonical (sorted-lambda) relabelling selects, where poles of the linear
//! `s3` terms lie within 1e-5 of their zeros; most other orientations give
//! 1e-14 (OneLOop's own, unpermuted one 3e-14). Reference: OneLOop 3.6 in
//! quadruple precision.
#![cfg(feature = "generated-evaluators")]
use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend,
};
use symbolica::domains::float::{Complex, Float, RealLike};

/// `[p1², p2², p3², p4², s12, s23, m0², m1², m2², m3², mu²]`.
const ARGUMENTS: [f64; 11] = [
    202500.0,
    15625.0,
    15625.0,
    15625.0,
    93305.84790585532,
    62501.31001345166,
    29929.0,
    29929.0,
    29929.0,
    29929.0,
    1e4,
];
const EXPECTED: [f64; 2] = [-1.907_030_378_722_972_01e-10, 8.834_687_875_568_863_24e-10];

/// Relative tolerances in binary64 and at 256 bits.
const BINARY64: f64 = 5e-11;
const PRECISE: f64 = 1e-15;

/// Every evaluation route: native binary64 (hand-written kernels and generated
/// code), the interpreted expressions and SymJIT in binary64, and 256-bit
/// native code and interpreted expressions.
#[test]
fn box_near_a_vanishing_gram_determinant_matches_oneloop_on_every_route() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(check_every_route)
        .unwrap()
        .join()
        .unwrap();
}

fn check_every_route() {
    let mut failures = Vec::new();
    let mut check = |route: &str, value: Complex<f64>, tolerance: f64| {
        let error =
            (value.re - EXPECTED[0]).hypot(value.im - EXPECTED[1]) / EXPECTED[0].hypot(EXPECTED[1]);
        if !(error <= tolerance) {
            failures.push(format!("{route}: {value:?} vs {EXPECTED:?} ({error:e})"));
        }
    };
    let input = ARGUMENTS.map(|x| Complex::new(x, 0.0));
    for handwritten in [true, false] {
        let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
        evaluator.set_handwritten(handwritten);
        let mut output = [Complex::new(0.0, 0.0); 3];
        evaluator.evaluate(&input, &mut output).unwrap();
        check(
            &format!("native (hand-written: {handwritten})"),
            output[0],
            BINARY64,
        );
    }
    for backend in [EvaluationBackend::Expression, EvaluationBackend::SymJit] {
        let mut output = [Complex::new(0.0, 0.0); 3];
        evaluate_with_backend(ScalarIntegral::D0, &input, &mut output, backend).unwrap();
        check(backend.name(), output[0], BINARY64);
    }
    let bits = 256;
    let input = ARGUMENTS.map(|x| Complex::new(Float::with_val(bits, x), Float::new(bits)));
    let fresh =
        || core::array::from_fn::<_, 3, _>(|_| Complex::new(Float::new(bits), Float::new(bits)));
    let mut native =
        NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::D0, bits).unwrap();
    let mut output = fresh();
    native.evaluate(&input, &mut output).unwrap();
    let value = Complex::new(output[0].re.to_f64(), output[0].im.to_f64());
    check("native 256 bits", value, PRECISE);
    let mut expression = PrecisionEvaluator::with_binary_precision_and_backend(
        ScalarIntegral::D0,
        bits,
        EvaluationBackend::Expression,
    )
    .unwrap();
    let mut output = fresh();
    expression.evaluate(&input, &mut output).unwrap();
    let value = Complex::new(output[0].re.to_f64(), output[0].im.to_f64());
    check("expression 256 bits", value, PRECISE);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
