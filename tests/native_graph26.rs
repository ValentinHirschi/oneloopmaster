//! A gg -> HHH top-loop box at which the 2024 OneLOop development snapshot
//! (OneLOopBridge `avh_olo.f90`, dated 2024-11-23) returns a different real
//! part (-2.0606e-9, also in its quad precision). OneLOop 3.4 (bundled with
//! Ninja), OneLOop 3.6 and both COLLIER branches (COLI and DD) agree with the
//! value below, as do this crate's 128-, 256- and 512-bit evaluations.
#![cfg(feature = "generated-evaluators")]
use oneloop::{NativeEvaluator, ScalarIntegral};
use symbolica::domains::float::Complex;

#[test]
fn graph26_box_matches_independent_references() {
    // D0(p1²..p4², s12, s23; m², m², m², m²; mu²) at sqrt(s) = 600 GeV,
    // propagators {1,2,3,4} of the pentagon (gluon p1² = 0, Higgs p² = 125²).
    let input = [
        0.0, 15625.0, 15625.0, -104375.0, -104375.0, 135625.0, 29929.0, 29929.0, 29929.0,
        29929.0, 10000.0,
    ]
    .map(|x| Complex::new(x, 0.0));
    let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
    let mut output = [Complex::new(0.0, 0.0); 3];
    evaluator.evaluate(&input, &mut output).unwrap();
    // COLLIER COLI: (2.86378976453069544e-10, 1.85745093171721722e-10).
    let expected = Complex::new(2.863_789_764_530_695_6e-10, 1.857_450_931_717_218_9e-10);
    let error = (output[0].re - expected.re).hypot(output[0].im - expected.im);
    assert!(error < 1e-12 * expected.re.hypot(expected.im), "{:?}", output[0]);
    assert_eq!(output[1], Complex::new(0.0, 0.0));
    assert_eq!(output[2], Complex::new(0.0, 0.0));
}
