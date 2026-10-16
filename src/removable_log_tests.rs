//! Independent Taylor controls use derivatives of the nonsingular numerator log.
#![allow(clippy::suspicious_arithmetic_impl, clippy::suspicious_op_assign_impl)]

use super::*;
use symbolica::{
    create_hyperdual_from_depths,
    domains::{
        dual::HyperDual,
        float::{Complex, Float, FloatLike, Real, SingleFloat},
    },
};

create_hyperdual_from_depths!(OrderEight, [8]);

fn c(re: f64, im: f64) -> Complex<f64> {
    Complex::new(re, im)
}

fn close(actual: Complex<f64>, expected: Complex<f64>) {
    assert!(
        (actual.re - expected.re).hypot(actual.im - expected.im)
            < 5e-12 * (1. + expected.re.hypot(expected.im)),
        "{actual} != {expected}"
    );
}

// If (1-x) f(x) = log(x), the coefficients of f(1+h) follow from
// derivatives of log(1+h) one order higher. Symbolica performs all derivatives.
// This oracle does not use the series operation used by production preparation.
pub(crate) fn oracle_polynomial(depth: usize, offset: &Atom) -> Atom {
    let h = symbolica::symbol!("removable_oracle_h");
    let mut derivative = (Atom::num(1) + h.to_atom()).log();
    let mut result = Atom::num(0);
    let mut factorial = 1usize;
    for degree in 0..=depth {
        derivative = derivative.derivative(h);
        factorial *= degree + 1;
        let coefficient = -derivative.replace(h).with(0) / Atom::num(factorial);
        result += coefficient * offset.pow(degree as i64);
    }
    result
}

#[test]
fn euclidean_c0_preserves_primal_and_zero_seeds() {
    let context = OneLoopExpressions::new();
    let parameters: Vec<_> = (0..7)
        .map(|i| symbolica::symbol!(format!("c0_primal_control_{i}")).to_atom())
        .collect();
    let series = context.c0(
        [&parameters[0], &parameters[1], &parameters[2]],
        [&parameters[3], &parameters[4], &parameters[5]],
        &parameters[6],
    );
    let scalar = context
        .evaluator(series.coefficients(), &parameters)
        .unwrap();
    let exact = context
        .dual_evaluator(
            series.coefficients(),
            &parameters,
            HyperDual::<f64>::new(vec![
                vec![0, 0],
                vec![1, 0],
                vec![0, 1],
                vec![2, 0],
                vec![1, 1],
            ]),
        )
        .unwrap();
    let values = [-1., -5., -4., 2., 3., 5., 1.].map(|v| c(v, 0.));
    let mut input = vec![c(0., 0.); 35];
    for (row, value) in input.chunks_exact_mut(5).zip(values) {
        row[0] = value;
    }
    let mut expected = [c(0., 0.); 3];
    let mut output = [c(0., 0.); 15];
    scalar
        .clone()
        .map_coeff(&|v| c(v.re.to_f64(), v.im.to_f64()))
        .evaluate(&values, &mut expected);
    exact
        .clone()
        .map_coeff(&|v| c(v.re.to_f64(), v.im.to_f64()))
        .evaluate(&input, &mut output);
    for (components, expected) in output.chunks_exact(5).zip(expected) {
        close(components[0], expected);
        assert!(components[1..].iter().all(|v| *v == c(0., 0.)));
    }
    let map = |v: &Complex<symbolica::domains::rational::Rational>| {
        Complex::new(v.re.to_multi_prec_float(256), v.im.to_multi_prec_float(256))
    };
    let precise =
        |v: &Complex<f64>| Complex::new(Float::with_val(256, v.re), Float::with_val(256, v.im));
    let zero = Complex::new(Float::new(256), Float::new(256));
    let mut expected = vec![zero.clone(); 3];
    let mut output = vec![zero; 15];
    scalar.map_coeff_with_prec(&map, 256).evaluate(
        &values.iter().map(precise).collect::<Vec<_>>(),
        &mut expected,
    );
    exact
        .map_coeff_with_prec(&map, 256)
        .evaluate(&input.iter().map(precise).collect::<Vec<_>>(), &mut output);
    let tolerance = Float::with_val(256, 2).pow(180).inv();
    for (components, expected) in output.chunks_exact(5).zip(expected) {
        assert!((components[0].re.clone() - expected.re).norm() < tolerance);
        assert!((components[0].im.clone() - expected.im).norm() < tolerance);
        assert!(
            components[1..]
                .iter()
                .all(|v| v.re.is_zero() && v.im.is_zero())
        );
    }
    // This checks evaluator semantics, not a general C0 derivative domain.
    // Coincident dilogarithm differences still need analytic neighborhoods.
}

#[test]
fn log_ratio_degree_eight_static_dynamic_and_precision() {
    let z = symbolica::symbol!("removable_z").to_atom();
    let expression = SheetAtom::lower(z.clone()).log_over_one_minus();
    let context = OneLoopExpressions::new();
    // Compile the zero-order map first: independently prepared higher-order
    // maps must not inherit its constant local branch through helper caching.
    let mut scalar = context
        .evaluator(std::slice::from_ref(&expression), std::slice::from_ref(&z))
        .unwrap()
        .map_coeff(&|v| c(v.re.to_f64(), v.im.to_f64()));
    close(scalar.evaluate_single(&[c(1., 0.)]), c(-1., 0.));
    let dynamic = context
        .dual_evaluator(
            std::slice::from_ref(&expression),
            std::slice::from_ref(&z),
            HyperDual::<f64>::new((0..=8).map(|i| vec![i]).collect()),
        )
        .unwrap();
    let static_shape = context
        .dual_evaluator(&[expression], &[z], OrderEight { values: [0.; 9] })
        .unwrap();
    let t = symbolica::symbol!("removable_t");
    let oracle = oracle_polynomial(8, &(sheet_exact::i() * t.to_atom()));
    let expected: Vec<_> = (0..=8)
        .map(|degree| {
            let mut value = oracle.clone();
            for _ in 0..degree {
                value = value.derivative(t);
            }
            value.replace(t).with(0) / Atom::num((1..=degree).product::<usize>())
        })
        .collect();
    let oracle = Atom::evaluator_multiple(&expected, &[] as &[Atom])
        .build()
        .unwrap();
    let input = [
        c(1., 0.),
        c(0., 1.),
        c(0., 0.),
        c(0., 0.),
        c(0., 0.),
        c(0., 0.),
        c(0., 0.),
        c(0., 0.),
        c(0., 0.),
    ];
    let mut expected = [c(0., 0.); 9];
    oracle
        .clone()
        .map_coeff(&|v| c(v.re.to_f64(), v.im.to_f64()))
        .evaluate(&[], &mut expected);
    for exact in [dynamic.clone(), static_shape] {
        let mut output = [c(0., 0.); 9];
        exact
            .map_coeff(&|v| c(v.re.to_f64(), v.im.to_f64()))
            .evaluate(&input, &mut output);
        for (actual, expected) in output.into_iter().zip(expected) {
            close(actual, expected);
        }
    }
    let map = |v: &Complex<symbolica::domains::rational::Rational>| {
        Complex::new(v.re.to_multi_prec_float(256), v.im.to_multi_prec_float(256))
    };
    let mut precise = dynamic.map_coeff_with_prec(&map, 256);
    let mut oracle = oracle.map_coeff_with_prec(&map, 256);
    let zero = Complex::new(Float::new(256), Float::new(256));
    let input: Vec<_> = input
        .iter()
        .map(|v| Complex::new(Float::with_val(256, v.re), Float::with_val(256, v.im)))
        .collect();
    let mut actual = vec![zero.clone(); 9];
    let mut expected = vec![zero; 9];
    precise.evaluate(&input, &mut actual);
    oracle.evaluate(&[], &mut expected);
    let tolerance = Float::with_val(256, 2).pow(180).inv();
    for (actual, expected) in actual.into_iter().zip(expected) {
        assert!((actual.re - expected.re).norm() < tolerance);
        assert!((actual.im - expected.im).norm() < tolerance);
    }
}

#[test]
fn log_ratio_mixed_transverse_lanes_and_lazy_branches() {
    let z = symbolica::symbol!("removable_mixed_z").to_atom();
    let shape = vec![
        vec![0, 0],
        vec![1, 0],
        vec![0, 1],
        vec![2, 0],
        vec![1, 1],
        vec![0, 2],
        vec![3, 0],
        vec![2, 1],
        vec![1, 2],
        vec![3, 1],
        vec![2, 2],
        vec![3, 2],
    ];
    let context = OneLoopExpressions::new();
    let expression = SheetAtom::lower(z.clone()).log_over_one_minus();
    let exact = context
        .dual_evaluator(
            std::slice::from_ref(&expression),
            std::slice::from_ref(&z),
            HyperDual::<f64>::new(shape.clone()),
        )
        .unwrap();
    let mut evaluator = exact
        .clone()
        .map_coeff(&|v| c(v.re.to_f64(), v.im.to_f64()))
        .into_batch_evaluator();
    let (t, u) = (
        symbolica::symbol!("removable_mixed_t"),
        symbolica::symbol!("removable_mixed_u"),
    );
    let coefficients = [
        (1, 0),
        (0, 1),
        (2, -1),
        (1, 2),
        (3, 1),
        (-1, 0),
        (0, 0),
        (2, -1),
        (0, 0),
        (0, 0),
        (0, 0),
        (0, 0),
    ];
    let offset =
        coefficients
            .iter()
            .zip(&shape)
            .skip(1)
            .fold(Atom::num(0), |sum, (&(re, im), powers)| {
                sum + (Atom::num(re) + sheet_exact::i() * im)
                    * t.to_atom().pow(powers[0] as i64)
                    * u.to_atom().pow(powers[1] as i64)
            });
    let polynomial = oracle_polynomial(5, &offset);
    let expected: Vec<_> = shape
        .iter()
        .map(|power| {
            let mut expression = polynomial.clone();
            for _ in 0..power[0] {
                expression = expression.derivative(t);
            }
            for _ in 0..power[1] {
                expression = expression.derivative(u);
            }
            expression.replace(t).with(0).replace(u).with(0)
                / Atom::num((1..=power[0]).product::<usize>() * (1..=power[1]).product::<usize>())
        })
        .collect();
    let mut oracle = Atom::evaluator_multiple(&expected, &[] as &[Atom])
        .build()
        .unwrap()
        .map_coeff(&|v| c(v.re.to_f64(), v.im.to_f64()));
    let mut expected = vec![c(0., 0.); shape.len()];
    oracle.evaluate(&[], &mut expected);
    let at_one: Vec<_> = coefficients
        .map(|(re, im)| c(f64::from(re), f64::from(im)))
        .to_vec();
    // A constant off-axis lane takes the quotient branch; alternating it with
    // the removable lane checks branch compaction, tail order and reactivation.
    let mut off_axis = vec![c(0., 0.); shape.len()];
    off_axis[0] = c(2., 1.);
    let mut scalar = context
        .evaluator(&[expression], &[z])
        .unwrap()
        .map_coeff(&|v| c(v.re.to_f64(), v.im.to_f64()));
    let off_axis_expected = scalar.evaluate_single(&off_axis[..1]);
    for rows in [0, 1, 9, 2, 17, 1] {
        let inputs: Vec<_> = (0..rows)
            .flat_map(|row| {
                if row % 2 == 0 {
                    at_one.clone()
                } else {
                    off_axis.clone()
                }
            })
            .collect();
        let mut outputs = vec![c(99., 99.); rows * shape.len()];
        evaluator
            .evaluate_batch(&inputs, &mut outputs, rows)
            .unwrap();
        for (row, output) in outputs.chunks(shape.len()).enumerate() {
            for (component, &actual) in output.iter().enumerate() {
                close(
                    actual,
                    if row % 2 == 0 {
                        expected[component]
                    } else if component == 0 {
                        off_axis_expected
                    } else {
                        c(0., 0.)
                    },
                );
            }
        }
    }
}

#[test]
fn log_ratio_unit_argument_on_wound_sheet_is_not_removable() {
    let value = symbolica::symbol!("__olo_sheet_value").to_atom();
    let phase = symbolica::symbol!("__olo_sheet_phase").to_atom();
    let odd = symbolica::symbol!("__olo_sheet_odd").to_atom();
    let expression =
        symbolica::symbol!("__olo_sheet_log_over_one_minus").call((&value, &phase, &odd));
    let mut evaluator = OneLoopExpressions::new()
        .evaluator(&[expression], &[value, phase, odd])
        .unwrap()
        .map_coeff(&|v| c(v.re.to_f64(), v.im.to_f64()));
    close(
        evaluator.evaluate_single(&[c(1., 0.), c(0., 0.), c(0., 0.)]),
        c(-1., 0.),
    );
    let pole = evaluator.evaluate_single(&[c(1., 0.), c(2., 0.), c(0., 0.)]);
    assert!(
        !pole.re.is_finite() || !pole.im.is_finite(),
        "a nonremovable sheet pole was replaced by {pole}"
    );
}

#[test]
fn log_ratio_mixed_near_unit_keeps_small_complex_displacement() {
    let z = symbolica::symbol!("removable_near_z").to_atom();
    let t = symbolica::symbol!("removable_near_t");
    let u = symbolica::symbol!("removable_near_u");
    let shape = vec![
        vec![0, 0],
        vec![1, 0],
        vec![0, 1],
        vec![2, 0],
        vec![1, 1],
        vec![0, 2],
        vec![2, 1],
    ];
    let i = sheet_exact::i();
    let input = vec![
        Atom::num(1) + Atom::num(1) / Atom::num(2).pow(20) - &i / Atom::num(2).pow(32),
        i.clone(),
        2 - &i,
        Atom::num(0),
        &i - 1,
        Atom::num(0),
        Atom::num(0),
    ];
    let polynomial = input
        .iter()
        .zip(&shape)
        .fold(Atom::num(0), |sum, (value, power)| {
            sum + value * t.to_atom().pow(power[0] as i64) * u.to_atom().pow(power[1] as i64)
        });
    let ordinary = polynomial.log() / (Atom::num(1) - &polynomial);
    let expected: Vec<_> = shape
        .iter()
        .map(|power| {
            let mut value = ordinary.clone();
            for _ in 0..power[0] {
                value = value.derivative(t);
            }
            for _ in 0..power[1] {
                value = value.derivative(u);
            }
            value.replace(t).with(0).replace(u).with(0)
                / Atom::num((1..=power[0]).product::<usize>() * (1..=power[1]).product::<usize>())
        })
        .collect();
    let map = |bits| {
        move |value: &Complex<symbolica::domains::rational::Rational>| {
            Complex::new(
                value.re.to_multi_prec_float(bits),
                value.im.to_multi_prec_float(bits),
            )
        }
    };
    let zero = |bits| Complex::new(Float::new(bits), Float::new(bits));
    let mut inputs = vec![zero(256); shape.len()];
    Atom::evaluator_multiple(&input, &[] as &[Atom])
        .build()
        .unwrap()
        .map_coeff_with_prec(&map(256), 256)
        .evaluate(&[], &mut inputs);
    let expression = SheetAtom::lower(z.clone()).log_over_one_minus();
    let mut evaluator = OneLoopExpressions::new()
        .dual_evaluator(&[expression], &[z], HyperDual::<f64>::new(shape.clone()))
        .unwrap()
        .map_coeff_with_prec(&map(256), 256);
    let mut actual = vec![zero(256); shape.len()];
    evaluator.evaluate(&inputs, &mut actual);
    let mut expected_values = vec![zero(384); shape.len()];
    Atom::evaluator_multiple(&expected, &[] as &[Atom])
        .build()
        .unwrap()
        .map_coeff_with_prec(&map(384), 384)
        .evaluate(&[], &mut expected_values);
    let tolerance = Float::with_val(256, 2).pow(140).inv();
    for (actual, expected) in actual.into_iter().zip(expected_values) {
        assert!((actual.re - expected.re).norm() < tolerance);
        assert!((actual.im - expected.im).norm() < tolerance);
    }
}
