use oneloop::{
    NativeTensorEvaluator, PoleCoefficient, PoleSelection, ScalarIntegral, TensorIntegral,
};
use symbolica::domains::float::{Complex, DoubleFloat, Float, FloatLike, RealLike};

type C = Complex<f64>;

fn c(re: f64) -> C {
    C::new(re, 0.)
}

fn bits(values: &[C]) -> Vec<(u64, u64)> {
    values
        .iter()
        .map(|z| (z.re.to_bits(), z.im.to_bits()))
        .collect()
}

fn close(a: C, b: C) {
    assert!(
        (a.re - b.re).hypot(a.im - b.im) <= 2e-13 * (1. + b.re.hypot(b.im)),
        "{a:?} != {b:?}"
    );
}

#[test]
fn native_tensor_scaleless_bubbles_keep_both_origins_through_rank_four() {
    for family in [
        TensorIntegral::B11,
        TensorIntegral::B111,
        TensorIntegral::B1111,
    ] {
        for binary64 in [false, true] {
            let mut evaluator = NativeTensorEvaluator::<f64>::new(family).unwrap();
            evaluator.set_binary64(binary64).unwrap();
            let point = [c(0.), c(0.), c(0.), c(13.)];
            let mut output = vec![c(99.); family.outputs()];
            for selection in [
                PoleSelection::BOTH,
                PoleSelection::UV,
                PoleSelection::IR,
                PoleSelection::NONE,
                PoleSelection::BOTH,
            ] {
                evaluator.set_pole_selection(selection);
                assert_eq!(evaluator.pole_selection(), selection);
                evaluator.evaluate(&point, &mut output).unwrap();
                for (rank, triple) in output.chunks_exact(3).enumerate() {
                    let residue = (-1_f64).powi(rank as i32) / (rank + 1) as f64;
                    let expected = match selection {
                        PoleSelection::UV => residue,
                        PoleSelection::IR => -residue,
                        _ => 0.,
                    };
                    assert_eq!(triple[0], c(0.));
                    assert_eq!(triple[2], c(0.));
                    close(triple[1], c(expected));
                }
            }
        }
    }
}

#[test]
fn tensor_selection_preserves_all_finite_terms_and_combined_bits() {
    for (family, point) in [
        (TensorIntegral::An, vec![C::new(7., -0.2), c(3.)]),
        (
            TensorIntegral::B1111,
            vec![c(-11.), C::new(4., -0.3), c(9.), c(3.)],
        ),
        (TensorIntegral::B1111, vec![c(0.), c(4.), c(9.), c(3.)]),
        (TensorIntegral::B1111, vec![c(17.), c(0.), c(0.), c(3.)]),
    ] {
        let mut evaluator = NativeTensorEvaluator::<f64>::new(family).unwrap();
        let mut total = vec![c(0.); family.outputs()];
        evaluator.evaluate(&point, &mut total).unwrap();
        for selection in [
            PoleSelection::UV,
            PoleSelection::IR,
            PoleSelection::NONE,
            PoleSelection::BOTH,
        ] {
            evaluator.set_pole_selection(selection);
            let mut selected = vec![c(0.); family.outputs()];
            evaluator.evaluate(&point, &mut selected).unwrap();
            for (a, b) in selected.chunks_exact(3).zip(total.chunks_exact(3)) {
                assert_eq!(bits(&a[..1]), bits(&b[..1]));
                close(a[1], if selection.uv { b[1] } else { c(0.) });
                close(a[2], c(0.));
            }
            if selection == PoleSelection::BOTH {
                assert_eq!(bits(&selected), bits(&total));
            }
        }
    }
}

#[test]
fn multiprecision_bubble_residues_are_not_rounded_through_binary64() {
    let z = Complex::new(Float::new(256), Float::new(256));
    let one = Complex::new(Float::with_val(256, 1), Float::new(256));
    let mut evaluator =
        NativeTensorEvaluator::<Float>::with_binary_precision(TensorIntegral::B1111, 256)
            .unwrap()
            .with_pole_selection(PoleSelection::UV);
    let input = [z.clone(), z.clone(), z.clone(), one];
    let mut uv = vec![z.clone(); 15];
    let mut ir = uv.clone();
    evaluator.evaluate(&input, &mut uv).unwrap();
    evaluator.set_pole_selection(PoleSelection::IR);
    evaluator.evaluate(&input, &mut ir).unwrap();
    for rank in 0..=4 {
        let pole = &uv[3 * rank + 1];
        assert!((pole.clone() + &ir[3 * rank + 1]).is_fully_zero());
        let sign = if rank % 2 == 0 { 1 } else { -1 };
        let expected = Float::with_val(256, sign) / Float::with_val(256, rank + 1);
        assert_eq!(pole.re.to_rational(), expected.to_rational());
        assert!(pole.im.is_fully_zero());
        if rank == 2 || rank == 4 {
            assert_ne!(
                pole.re.to_rational(),
                Float::with_val(256, pole.re.to_f64()).to_rational()
            );
        }
    }
    let zero = Complex::new(DoubleFloat::from(0.), DoubleFloat::from(0.));
    let one = Complex::new(DoubleFloat::from(1.), DoubleFloat::from(0.));
    let mut evaluator = NativeTensorEvaluator::<DoubleFloat>::new(TensorIntegral::B11)
        .unwrap()
        .with_pole_selection(PoleSelection::UV);
    let input = [zero, zero, zero, one];
    let mut uv = [zero; 9];
    let mut ir = [zero; 9];
    evaluator.evaluate(&input, &mut uv).unwrap();
    evaluator.set_pole_selection(PoleSelection::IR);
    evaluator.evaluate(&input, &mut ir).unwrap();
    assert!((uv[7] + ir[7]).is_fully_zero());
    assert_eq!(uv[7].re, DoubleFloat::from(1.) / DoubleFloat::from(3.));
    assert_ne!(uv[7].re, DoubleFloat::from(1_f64 / 3.));
}

#[test]
fn selection_keeps_transactional_shape_and_domain_errors() {
    let mut evaluator = NativeTensorEvaluator::<f64>::new(TensorIntegral::B1111)
        .unwrap()
        .with_pole_selection(PoleSelection::UV);
    let mut output = vec![C::new(91., -37.); 30];
    let before = bits(&output);
    let mut input = vec![c(-11.), c(4.), c(9.), c(3.)];
    input.extend([c(-11.), C::new(4., 0.2), c(9.), c(3.)]);
    assert!(evaluator.evaluate_batch(&input, &mut output, 2).is_err());
    assert_eq!(bits(&output), before);
    assert!(
        evaluator
            .evaluate_batch(&input[..7], &mut output, 2)
            .is_err()
    );
    assert_eq!(bits(&output), before);
    evaluator.evaluate_batch(&[], &mut [], 0).unwrap();
    let mut triple = [c(9.); 3];
    assert!(
        PoleSelection::UV
            .apply(PoleCoefficient::BubbleMoment { rank: 5 }, &mut triple)
            .is_err()
    );
    assert_eq!(triple, [c(9.); 3]);
}

#[test]
fn public_selector_supports_interpreted_results_without_changing_finite_part() {
    let mut triangle = [C::new(2., -3.), c(4.), c(5.)];
    PoleSelection::UV
        .apply(PoleCoefficient::Scalar(ScalarIntegral::C0), &mut triangle)
        .unwrap();
    assert_eq!(triangle, [C::new(2., -3.), c(0.), c(0.)]);
    let mut derivative = [c(7.), c(-0.125), c(0.)];
    PoleSelection::IR
        .apply(
            PoleCoefficient::Scalar(ScalarIntegral::DB0),
            &mut derivative,
        )
        .unwrap();
    assert_eq!(derivative, [c(7.), c(-0.125), c(0.)]);
}

#[cfg(feature = "generated-evaluators")]
#[test]
fn scalar_selection_covers_scaleless_bubble_and_infrared_derivative() {
    use oneloop::NativeEvaluator;
    for (family, input, uv_expected, ir_expected) in [
        (
            ScalarIntegral::A0,
            vec![C::new(4., -0.2), c(3.)],
            C::new(4., -0.2),
            c(0.),
        ),
        (
            ScalarIntegral::B0,
            vec![c(0.), c(0.), c(0.), c(3.)],
            c(1.),
            c(-1.),
        ),
        (
            ScalarIntegral::B0,
            vec![c(-11.), c(4.), c(9.), c(3.)],
            c(1.),
            c(0.),
        ),
        (
            ScalarIntegral::DB0,
            vec![c(4.), c(4.), c(0.), c(3.)],
            c(0.),
            c(-0.125),
        ),
        (
            ScalarIntegral::C0,
            vec![c(0.), c(0.), c(-11.), c(0.), c(0.), c(0.), c(3.)],
            c(0.),
            c((3_f64 / 11.).ln() / -11.),
        ),
    ] {
        let mut evaluator = NativeEvaluator::<f64>::new(family).unwrap();
        let mut total = [c(0.); 3];
        evaluator.evaluate(&input, &mut total).unwrap();
        let mut selected = [c(0.); 3];
        for selection in [
            PoleSelection::UV,
            PoleSelection::IR,
            PoleSelection::NONE,
            PoleSelection::BOTH,
        ] {
            evaluator.set_pole_selection(selection);
            assert_eq!(evaluator.pole_selection(), selection);
            evaluator.evaluate(&input, &mut selected).unwrap();
            assert_eq!(bits(&selected[..1]), bits(&total[..1]));
            close(
                selected[1],
                match selection {
                    PoleSelection::UV => uv_expected,
                    PoleSelection::IR => ir_expected,
                    PoleSelection::NONE => c(0.),
                    _ => total[1],
                },
            );
            close(selected[2], if selection.ir { total[2] } else { c(0.) });
            if selection == PoleSelection::BOTH {
                assert_eq!(bits(&selected), bits(&total));
            }
        }
        close(uv_expected + ir_expected, total[1]);
        let before = bits(&selected);
        assert!(
            evaluator
                .evaluate(&input[..input.len() - 1], &mut selected)
                .is_err()
        );
        assert_eq!(bits(&selected), before);
    }
}
