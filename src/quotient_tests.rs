//! Quotient and composed-log controls keep independent numerator/denominator jets.
use super::*;
use symbolica::{
    domains::{
        dual::HyperDual,
        float::{Complex, Float, FloatLike, Real},
    },
    evaluate::Dualizer,
};

const SHAPE: [[usize; 2]; 7] = [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [3, 0], [2, 1]];

fn c(re: f64, im: f64) -> Complex<f64> {
    Complex::new(re, im)
}

#[test]
fn equal_nonzero_sheet_arguments_keep_independent_mixed_jets() {
    let (x, y) = (
        symbolica::symbol!("quotient_test_x").to_atom(),
        symbolica::symbol!("quotient_test_y").to_atom(),
    );
    let ratio = &SheetAtom::lower(x.clone()) / &SheetAtom::lower(y.clone());
    let expressions = [ratio.value(), ratio.log(), ratio.log_over_one_minus()];
    let raw = &x / &y;
    let oracle = [
        raw.clone(),
        raw.log(),
        crate::removable_log_tests::oracle_polynomial(3, &(&raw - 1)),
    ];
    let shape = || HyperDual::<f64>::new(SHAPE.map(Vec::from).to_vec());
    let context = OneLoopExpressions::new();
    let exact = context
        .dual_evaluator(&expressions, &[x.clone(), y.clone()], shape())
        .unwrap();
    let oracle = Atom::evaluator_multiple(&oracle, &[x, y])
        .build()
        .unwrap()
        .vectorize(&Dualizer::new(shape(), vec![]))
        .unwrap();
    let map = |v: &Complex<symbolica::domains::rational::Rational>| c(v.re.to_f64(), v.im.to_f64());
    let mut actual = exact.clone().map_coeff(&map).into_batch_evaluator();
    let mut expected = oracle.clone().map_coeff(&map).into_batch_evaluator();
    let mut input = Vec::new();
    let primals = [(2., 0.), (-2., 0.), (3., 0.), (0.13, 0.71), (-0.13, -0.71)];
    let rows = primals.len();
    for (re, im) in primals {
        input.extend([
            c(re, im),
            c(0., 1.),
            c(0., 0.),
            c(1., 0.),
            c(0., 2.),
            c(-1., 1.),
            c(2., 1.),
        ]);
        input.extend([
            c(re, im),
            c(0., 0.),
            c(1., 0.),
            c(0., 1.),
            c(2., 0.),
            c(1., -1.),
            c(-1., 2.),
        ]);
    }
    let mut output = vec![c(0., 0.); rows * 3 * SHAPE.len()];
    let mut reference = output.clone();
    actual.evaluate_batch(&input, &mut output, rows).unwrap();
    expected
        .evaluate_batch(&input, &mut reference, rows)
        .unwrap();
    for row in output.chunks(3 * SHAPE.len()) {
        assert_eq!(row[0], c(1., 0.), "quotient primal must stay exactly one");
    }
    for (&actual, &expected) in output.iter().zip(&reference) {
        assert!(
            (actual.re - expected.re).hypot(actual.im - expected.im)
                < 5e-12 * (1. + expected.re.hypot(expected.im)),
            "{actual} != {expected}"
        );
    }
    // The first derivative of x/y is nonzero at equal primals; the old shortcut
    // returned literal one and erased this and every higher/mixed component.
    assert_eq!(output[1], c(0., 0.5));
    let map = |v: &Complex<symbolica::domains::rational::Rational>| {
        Complex::new(v.re.to_multi_prec_float(256), v.im.to_multi_prec_float(256))
    };
    let mut actual = exact.map_coeff_with_prec(&map, 256).into_batch_evaluator();
    let mut expected = oracle.map_coeff_with_prec(&map, 256).into_batch_evaluator();
    let inputs: Vec<_> = input
        .iter()
        .map(|v| Complex::new(Float::with_val(256, v.re), Float::with_val(256, v.im)))
        .collect();
    let zero = Complex::new(Float::new(256), Float::new(256));
    let mut output = vec![zero.clone(); reference.len()];
    let mut reference = vec![zero; reference.len()];
    actual.evaluate_batch(&inputs, &mut output, rows).unwrap();
    expected
        .evaluate_batch(&inputs, &mut reference, rows)
        .unwrap();
    let tolerance = Float::with_val(256, 2).pow(180).inv();
    for (actual, expected) in output.into_iter().zip(reference) {
        assert!((actual.re - expected.re).norm() < tolerance);
        assert!((actual.im - expected.im).norm() < tolerance);
    }
}

#[test]
fn derivative_constructor_does_not_invent_a_zero_over_zero_jet() {
    let (x, y) = (
        symbolica::symbol!("zero_quotient_x").to_atom(),
        symbolica::symbol!("zero_quotient_y").to_atom(),
    );
    let ratio = (&SheetAtom::lower(x.clone()) / &SheetAtom::lower(y.clone())).value();
    let context = OneLoopExpressions::new();
    let parameters = [x, y];
    let map = |v: &Complex<symbolica::domains::rational::Rational>| c(v.re.to_f64(), v.im.to_f64());
    let mut scalar = context
        .evaluator(std::slice::from_ref(&ratio), &parameters)
        .unwrap()
        .map_coeff(&map);
    assert_eq!(scalar.evaluate_single(&[c(0., 0.), c(0., 0.)]), c(1., 0.));
    let mut derivative = context
        .dual_evaluator(
            &[ratio],
            &parameters,
            HyperDual::<f64>::new(vec![vec![0], vec![1]]),
        )
        .unwrap()
        .map_coeff(&map);
    let mut output = [c(99., 99.); 2];
    derivative.evaluate(&[c(0., 0.), c(1., 0.), c(0., 0.), c(0., 1.)], &mut output);
    assert!(
        !output[0].re.is_finite() || !output[0].im.is_finite(),
        "undefined independent vanishing jets became {output:?}"
    );
}
