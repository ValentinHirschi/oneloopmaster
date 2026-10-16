//! Independent controls differentiate the defining logarithmic derivative.
use super::*;
use symbolica::{
    domains::{
        dual::HyperDual,
        float::{Complex, Float, FloatLike, Real, SingleFloat},
        integer::Integer,
    },
    evaluate::ExpressionEvaluator,
};

const SHAPE: [[usize; 2]; 11] = [
    [0, 0],
    [1, 0],
    [0, 1],
    [2, 0],
    [1, 1],
    [3, 0],
    [2, 1],
    [4, 0],
    [3, 1],
    [5, 0],
    [4, 1],
];

fn shape() -> HyperDual<f64> {
    HyperDual::new(SHAPE.map(Vec::from).to_vec())
}

fn c(re: f64, im: f64) -> Complex<f64> {
    Complex::new(re, im)
}

fn exact_inputs() -> [Atom; 2] {
    ["coincident_dilog_x", "coincident_dilog_y"].map(|s| symbolica::symbol!(s).to_atom())
}

fn row(first: Complex<f64>, second: Complex<f64>) -> Vec<Complex<f64>> {
    let mut row = vec![c(0., 0.); 2 * SHAPE.len()];
    row[0] = first;
    row[1] = c(0., 1.);
    row[SHAPE.len()] = second;
    row[SHAPE.len() + 2] = c(1., 0.);
    row
}

// D(x,y) = integral_0^1 F'(y+t(x-y)) dt gives the Taylor coefficient
// F^(a+b+1)(z)/(a+b+1)! at x=y=z for independent linear seeds. Symbolica
// differentiates F'(z)=log(z)/(1-z); this does not use the production series.
fn oracle(
    point: &Atom,
    negative_lip: i32,
) -> ExpressionEvaluator<Complex<symbolica::domains::rational::Rational>> {
    let z = symbolica::symbol!("coincident_dilog_oracle_z");
    let mut derivative = if point.is_one() {
        crate::removable_log_tests::oracle_polynomial(5, &z.to_atom())
    } else {
        let logarithm = if negative_lip != 0 {
            (-z.to_atom()).log() + negative_lip * sheet_exact::i() * Symbol::PI.to_atom()
        } else {
            z.to_atom().log()
        };
        logarithm / (Atom::num(1) - z.to_atom())
    };
    let center = if point.is_one() {
        Atom::num(0)
    } else {
        point.clone()
    };
    let mut derivatives = Vec::new();
    for degree in 0..=5 {
        derivatives.push(
            derivative.replace(z).with(center.clone()) / Atom::num(Integer::factorial(degree + 1)),
        );
        derivative = derivative.derivative(z);
    }
    let expressions: Vec<_> = SHAPE
        .iter()
        .map(|&[a, b]| &derivatives[a + b] * sheet_exact::i().pow(a as u64))
        .collect();
    Atom::evaluator_multiple(&expressions, &[] as &[Atom])
        .build()
        .unwrap()
}

#[test]
fn coincident_sheet_dilog_preserves_independent_jets_and_primal() {
    let [x, y] = exact_inputs();
    let expression =
        sheet_exact::divided_difference(&SheetAtom::lower(x.clone()), &SheetAtom::lower(y.clone()));
    let context = OneLoopExpressions::new();
    let exact = context
        .dual_evaluator(
            std::slice::from_ref(&expression),
            &[x.clone(), y.clone()],
            shape(),
        )
        .unwrap();
    let scalar = context.evaluator(&[expression], &[x, y]).unwrap();
    let f64_map =
        |v: &Complex<symbolica::domains::rational::Rational>| c(v.re.to_f64(), v.im.to_f64());
    let mut actual = exact.clone().map_coeff(&f64_map).into_batch_evaluator();
    let mut scalar = scalar.map_coeff(&f64_map);
    let points = [
        symbolica::parse!("2"),
        symbolica::parse!("1"),
        symbolica::parse!("-2"),
        symbolica::parse!("7/10+1i/5"),
        symbolica::parse!("1"),
    ];
    let mut inputs = Vec::new();
    let mut expected = Vec::new();
    let mut scalar_expected = Vec::new();
    for (index, point) in points.iter().enumerate() {
        let numeric = Complex::<f64>::try_from(point.to_float(16)).unwrap();
        inputs.extend(row(numeric, numeric));
        scalar_expected.push(scalar.evaluate_single(&[numeric, numeric]));
        let mut output = vec![c(0., 0.); SHAPE.len()];
        oracle(point, if index == 2 { -1 } else { 0 })
            .map_coeff(&f64_map)
            .evaluate(&[], &mut output);
        expected.extend(output);
    }
    let mut output = vec![c(0., 0.); expected.len()];
    actual
        .evaluate_batch(&inputs, &mut output, points.len())
        .unwrap();
    for (actual, expected) in output.chunks_exact(SHAPE.len()).zip(scalar_expected) {
        assert_eq!(
            actual[0], expected,
            "coincidence must preserve the scalar primal"
        );
    }
    for (index, (actual, expected)) in output.into_iter().zip(expected).enumerate() {
        assert!(
            (actual.re - expected.re).hypot(actual.im - expected.im)
                < 2e-10 * (1. + expected.re.hypot(expected.im)),
            "component {index}: {actual} != {expected}"
        );
    }
    let bits = 256;
    let precise_map = |v: &Complex<symbolica::domains::rational::Rational>| {
        Complex::new(
            v.re.to_multi_prec_float(bits),
            v.im.to_multi_prec_float(bits),
        )
    };
    let mut actual = exact
        .map_coeff_with_prec(&precise_map, bits)
        .into_batch_evaluator();
    let mut inputs = Vec::new();
    let mut expected = Vec::new();
    let zero = Complex::new(Float::new(bits), Float::new(bits));
    for (index, point) in points.iter().enumerate() {
        let numeric = Complex::<Float>::try_from(point.to_float(80)).unwrap();
        let mut input = vec![zero.clone(); 2 * SHAPE.len()];
        input[0] = numeric.clone();
        input[1] = Complex::new(Float::new(bits), Float::with_val(bits, 1));
        input[SHAPE.len()] = numeric;
        input[SHAPE.len() + 2] = Complex::new(Float::with_val(bits, 1), Float::new(bits));
        inputs.extend(input);
        let mut output = vec![zero.clone(); SHAPE.len()];
        oracle(point, if index == 2 { -1 } else { 0 })
            .map_coeff_with_prec(&precise_map, bits)
            .evaluate(&[], &mut output);
        expected.extend(output);
    }
    let mut output = vec![zero; expected.len()];
    actual
        .evaluate_batch(&inputs, &mut output, points.len())
        .unwrap();
    let tolerance = Float::with_val(bits, 2).pow(170).inv();
    for (actual, expected) in output.into_iter().zip(expected) {
        assert!((actual.re - expected.re).norm() < tolerance);
        assert!((actual.im - expected.im).norm() < tolerance);
    }
}

#[test]
fn incompatible_sheets_and_zero_do_not_get_removable_derivatives() {
    let parameters: Vec<_> = (0..6)
        .map(|i| symbolica::symbol!(format!("dilog_germ_{i}")).to_atom())
        .collect();
    let expression = symbolica::symbol!("__olo_sheet_dilog_difference").call(parameters.as_slice());
    let context = OneLoopExpressions::new();
    let mut actual = context
        .dual_evaluator(&[expression], &parameters, shape())
        .unwrap()
        .map_coeff(&|v| c(v.re.to_f64(), v.im.to_f64()))
        .into_batch_evaluator();
    let values = [
        [1., 0., 0., 1., 0., 0.],
        [0., 0., 0., 0., 0., 0.],
        [2., -1., 1., 2., 1., 1.],
        [1., 2., 0., 1., 2., 0.],
        [1., 0., 0., 1., 0., 0.],
        [2., -1., 1., 2., -1., 1.],
        [2., 1., 1., 2., 1., 1.],
    ];
    let mut input = vec![c(0., 0.); values.len() * 6 * SHAPE.len()];
    for (row, values) in input.chunks_exact_mut(6 * SHAPE.len()).zip(values) {
        for (components, value) in row.chunks_exact_mut(SHAPE.len()).zip(values) {
            components[0] = c(value, 0.);
        }
        row[1] = c(0., 1. - 2. * values[2]);
        row[3 * SHAPE.len() + 2] = c(1. - 2. * values[5], 0.);
    }
    let mut output = vec![c(0., 0.); values.len() * SHAPE.len()];
    actual
        .evaluate_batch(&input, &mut output, values.len())
        .unwrap();
    for (index, row) in output.chunks_exact(SHAPE.len()).enumerate() {
        if index == 0 || index >= 4 {
            assert!(row.iter().all(|v| v.re.is_finite() && v.im.is_finite()));
        } else {
            assert!(
                row.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()),
                "nonanalytic row {index}: {row:?}"
            );
        }
    }
    for (index, lip) in [(5, -1), (6, 1)] {
        let mut expected = vec![c(0., 0.); SHAPE.len()];
        oracle(&Atom::num(-2), lip)
            .map_coeff(&|v| c(v.re.to_f64(), v.im.to_f64()))
            .evaluate(&[], &mut expected);
        for (actual, expected) in output[index * SHAPE.len()..(index + 1) * SHAPE.len()]
            .iter()
            .zip(expected)
        {
            assert!((actual.re - expected.re).hypot(actual.im - expected.im) < 2e-10);
        }
    }
}

#[test]
fn coincident_and_unit_endpoint_jets_match_converged_float_differences() {
    let [x, y] = exact_inputs();
    let expression =
        sheet_exact::divided_difference(&SheetAtom::lower(x.clone()), &SheetAtom::lower(y.clone()));
    let context = OneLoopExpressions::new();
    let parameters = [x, y];
    let shape = HyperDual::<f64>::new(SHAPE[..5].iter().copied().map(Vec::from).collect());
    let exact = context
        .dual_evaluator(std::slice::from_ref(&expression), &parameters, shape)
        .unwrap();
    let bits = 384;
    let map = |v: &Complex<symbolica::domains::rational::Rational>| {
        Complex::new(
            v.re.to_multi_prec_float(bits),
            v.im.to_multi_prec_float(bits),
        )
    };
    let mut actual = exact.map_coeff_with_prec(&map, bits);
    let mut scalar = context
        .evaluator(&[expression], &parameters)
        .unwrap()
        .map_coeff_with_prec(&map, bits);
    let number = |re, im| Complex::new(Float::with_val(bits, re), Float::with_val(bits, im));
    let zero = number(0., 0.);
    for (first, second) in [
        (c(1., 0.), c(1., 0.)),
        (c(2., 0.), c(2., 0.)),
        (c(1., 0.), c(2., 0.)),
        (c(2., 0.), c(1., 0.)),
        (c(1.5, 0.25), c(1.5, 0.25)),
    ] {
        let first = number(first.re, first.im);
        let second = number(second.re, second.im);
        let mut input = vec![zero.clone(); 10];
        input[0] = first.clone();
        input[1] = number(0., 1.);
        input[5] = second.clone();
        input[7] = number(1., 0.);
        let mut output = vec![zero.clone(); 5];
        actual.evaluate(&input, &mut output);
        let at_center = scalar.evaluate_single(&[first.clone(), second.clone()]);
        let mut previous: Option<Float> = None;
        for exponent in [12, 22, 32] {
            let h = Float::with_val(bits, 2).pow(exponent).inv();
            let real = Complex::new(h.clone(), Float::new(bits));
            let imaginary = Complex::new(Float::new(bits), h.clone());
            let plus = scalar.evaluate_single(&[&first + &imaginary, second.clone()]);
            let minus = scalar.evaluate_single(&[&first - &imaginary, second.clone()]);
            let dx = (&plus - &minus) / (&real * number(2., 0.));
            let dxx =
                (&plus + &minus - &at_center * number(2., 0.)) / (&real * &real * number(2., 0.));
            let pp = scalar.evaluate_single(&[&first + &imaginary, &second + &real]);
            let pm = scalar.evaluate_single(&[&first + &imaginary, &second - &real]);
            let mp = scalar.evaluate_single(&[&first - &imaginary, &second + &real]);
            let mm = scalar.evaluate_single(&[&first - &imaginary, &second - &real]);
            let mixed = (pp - pm - mp + mm) / (&real * &real * number(4., 0.));
            let mut error = Float::new(bits);
            for (numerical, analytic) in [dx, dxx, mixed]
                .into_iter()
                .zip([&output[1], &output[3], &output[4]])
            {
                error +=
                    (numerical.re - &analytic.re).norm() + (numerical.im - &analytic.im).norm();
            }
            if let Some(previous) = previous {
                assert!(
                    error < previous / Float::with_val(bits, 100_000),
                    "differences did not converge at ({first},{second}), h=2^-{exponent}: {error}"
                );
            }
            previous = Some(error);
        }
        assert!(previous.unwrap() < Float::with_val(bits, 2).pow(60).inv());
    }
}

#[test]
fn nonlinear_mixed_unit_seeds_match_native_polynomial_integral() {
    let [x, y] = exact_inputs();
    let expression =
        sheet_exact::divided_difference(&SheetAtom::lower(x.clone()), &SheetAtom::lower(y.clone()));
    let context = OneLoopExpressions::new();
    let exact = context
        .dual_evaluator(&[expression], &[x, y], shape())
        .unwrap();
    let s = symbolica::symbol!("nonlinear_dilog_s");
    let t = symbolica::symbol!("nonlinear_dilog_t");
    let u = symbolica::symbol!("nonlinear_dilog_integration");
    let dx = sheet_exact::i() * s.to_atom()
        + s.to_atom().pow(2)
        + 2 * s.to_atom() * t.to_atom()
        + t.to_atom().pow(3);
    let dy = t.to_atom() + sheet_exact::i() * s.to_atom().pow(2) - s.to_atom() * t.to_atom()
        + s.to_atom().pow(3);
    // Independently integrate the logarithmic derivative along the straight
    // segment. Symbolica's polynomial integrator performs the entire integral;
    // no dilogarithm series or divided-difference cancellation is used here.
    let integrand =
        crate::removable_log_tests::oracle_polynomial(5, &(&dy + u.to_atom() * (&dx - &dy)));
    let primitive = integrand
        .to_polynomial_in_vars::<u16>(u)
        .integrate(0)
        .flatten(false);
    let oracle = primitive.replace(u).with(1) - primitive.replace(u).with(0);
    let components = |expression: &Atom| -> Vec<Atom> {
        SHAPE
            .iter()
            .map(|&[a, b]| {
                let mut derivative = expression.clone();
                for _ in 0..a {
                    derivative = derivative.derivative(s);
                }
                for _ in 0..b {
                    derivative = derivative.derivative(t);
                }
                derivative.replace(s).with(0).replace(t).with(0)
                    / Atom::num(Integer::factorial(a as u32) * Integer::factorial(b as u32))
            })
            .collect()
    };
    let input: Vec<_> = [Atom::num(1) + dx, Atom::num(1) + dy]
        .iter()
        .flat_map(components)
        .collect();
    let expected = components(&oracle);
    let input = Atom::evaluator_multiple(&input, &[] as &[Atom])
        .build()
        .unwrap();
    let oracle = Atom::evaluator_multiple(&expected, &[] as &[Atom])
        .build()
        .unwrap();
    let bits = 256;
    let map = |v: &Complex<symbolica::domains::rational::Rational>| {
        Complex::new(
            v.re.to_multi_prec_float(bits),
            v.im.to_multi_prec_float(bits),
        )
    };
    let zero = Complex::new(Float::new(bits), Float::new(bits));
    let mut inputs = vec![zero.clone(); 2 * SHAPE.len()];
    input
        .map_coeff_with_prec(&map, bits)
        .evaluate(&[], &mut inputs);
    let mut expected = vec![zero.clone(); SHAPE.len()];
    oracle
        .map_coeff_with_prec(&map, bits)
        .evaluate(&[], &mut expected);
    // Follow the nonlinear seeds by the same primal with all zero seeds.
    let mut unseeded = vec![zero.clone(); 2 * SHAPE.len()];
    unseeded[0] = Complex::new(Float::with_val(bits, 1), Float::new(bits));
    unseeded[SHAPE.len()] = unseeded[0].clone();
    inputs.extend(unseeded);
    let mut output = vec![zero; 2 * SHAPE.len()];
    exact
        .map_coeff_with_prec(&map, bits)
        .into_batch_evaluator()
        .evaluate_batch(&inputs, &mut output, 2)
        .unwrap();
    let tolerance = Float::with_val(bits, 2).pow(180).inv();
    for (actual, expected) in output[..SHAPE.len()].iter().zip(expected) {
        assert!((actual.re.clone() - expected.re).norm() < tolerance);
        assert!((actual.im.clone() - expected.im).norm() < tolerance);
    }
    assert_eq!(output[SHAPE.len()].re, Float::with_val(bits, -1));
    assert!(
        output[SHAPE.len() + 1..]
            .iter()
            .all(|v| v.re.is_zero() && v.im.is_zero())
    );
}
