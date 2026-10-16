//! End-to-end derivative controls at nonsingular, massive spacelike points.
//! Gram-zero, threshold and projective singular limits are deliberately absent.
//!
//! Algebra audit: Symbolica exposes native differentiation and Matrix::solve;
//! their implementations/tests (evaluate/dual.rs, numerica/tensors/matrix.rs,
//! examples/solve_linear_system.rs) cover these operations. No reusable finite
//! difference stencil generator was found in Symbolica or FeynKit. Exact stencil
//! weights below are therefore solved by Symbolica, never by a local solver.

use oneloop::OneLoopExpressions;
use symbolica::{
    atom::{Atom, AtomCore},
    domains::{
        dual::HyperDual,
        float::{Complex, Float, FloatLike, Real, RealLike, SingleFloat},
        integer::Integer,
        rational::{Q, Rational},
    },
    evaluate::ExpressionEvaluator,
    tensors::matrix::Matrix,
};

const BITS: u32 = 384;
const POINTS: [(&str, [i64; 7]); 4] = [
    ("negative Källén", [-1, -5, -4, 2, 3, 5, 1]),
    ("positive Källén", [-1, -2, -9, 2, 3, 5, 1]),
    ("negative Källén/equal masses", [-1, -5, -4, 3, 3, 3, 1]),
    ("positive Källén/equal masses", [-1, -2, -9, 3, 3, 3, 1]),
];
const SHAPE: [[usize; 2]; 7] = [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [3, 0], [2, 1]];
type C = Complex<Float>;

fn c(re: i64, im: i64) -> C {
    Complex::new(Float::with_val(BITS, re), Float::with_val(BITS, im))
}

fn map(value: &Complex<Rational>) -> C {
    Complex::new(
        value.re.to_multi_prec_float(BITS),
        value.im.to_multi_prec_float(BITS),
    )
}

fn size(value: &C) -> Float {
    value.re.norm() + value.im.norm()
}

fn close(actual: &C, expected: &C, bits: u64, label: &str) -> Float {
    assert!(
        actual.re.is_finite() && actual.im.is_finite(),
        "{label}: {actual}"
    );
    let error = size(&(actual - expected));
    let bound = (Float::with_val(BITS, 1) + size(expected)) / Float::with_val(BITS, 2).pow(bits);
    assert!(
        error < bound,
        "{label}: error={error:e}, bound={bound:e}, actual={actual}, expected={expected}"
    );
    error
}

struct Evaluators {
    scalar: ExpressionEvaluator<C>,
    dual: ExpressionEvaluator<C>,
    double: ExpressionEvaluator<Complex<f64>>,
}

impl Evaluators {
    fn new() -> Self {
        let context = OneLoopExpressions::new();
        let p: Vec<_> = (0..7)
            .map(|j| symbolica::symbol!(format!("c0_derivative_input_{j}")).to_atom())
            .collect();
        let series = context.c0([&p[0], &p[1], &p[2]], [&p[3], &p[4], &p[5]], &p[6]);
        let scalar = context.evaluator(series.coefficients(), &p).unwrap();
        let dual = context
            .dual_evaluator(
                series.coefficients(),
                &p,
                HyperDual::<f64>::new(SHAPE.map(Vec::from).to_vec()),
            )
            .unwrap();
        Self {
            scalar: scalar.map_coeff_with_prec(&map, BITS),
            double: dual
                .clone()
                .map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64())),
            dual: dual.map_coeff_with_prec(&map, BITS),
        }
    }

    fn scalar(&mut self, input: &[C]) -> C {
        let mut output = vec![c(0, 0); 3];
        self.scalar.evaluate(input, &mut output);
        assert!(
            output.iter().all(|v| v.re.is_finite() && v.im.is_finite()),
            "scalar input={input:?}, output={output:?}"
        );
        assert!(output[1..].iter().all(|v| v.re.is_zero() && v.im.is_zero()));
        output.remove(0)
    }

    fn dual(&mut self, input: &[C]) -> Vec<C> {
        let mut output = vec![c(0, 0); 3 * SHAPE.len()];
        self.dual.evaluate(input, &mut output);
        assert!(
            output.iter().all(|v| v.re.is_finite() && v.im.is_finite()),
            "nonfinite Dual output={output:?}"
        );
        assert!(
            output[SHAPE.len()..]
                .iter()
                .all(|v| v.re.is_zero() && v.im.is_zero())
        );
        let mut double = vec![Complex::new(0., 0.); output.len()];
        let input64: Vec<_> = input
            .iter()
            .map(|v| Complex::new(v.re.to_f64(), v.im.to_f64()))
            .collect();
        self.double.evaluate(&input64, &mut double);
        for (actual, expected) in double.iter().zip(&output) {
            let actual = Complex::new(
                Float::with_val(BITS, actual.re),
                Float::with_val(BITS, actual.im),
            );
            close(&actual, expected, 30, "f64/Float384 components");
        }
        output.truncate(SHAPE.len());
        output
    }
}

fn input(base: &[C], first: &[C], second: &[C]) -> Vec<C> {
    base.iter()
        .zip(first)
        .zip(second)
        .flat_map(|((primal, a), b)| {
            [
                primal.clone(),
                a.clone(),
                b.clone(),
                c(0, 0),
                c(0, 0),
                c(0, 0),
                c(0, 0),
            ]
        })
        .collect()
}

fn scaling_oracle() -> ExpressionEvaluator<C> {
    let t = symbolica::symbol!("c0_scaling_t");
    let u = symbolica::symbol!("c0_scaling_u");
    let f = symbolica::symbol!("c0_scaling_value");
    let i = Atom::num(-1).sqrt();
    let coefficients = [
        Atom::num(1),
        i.clone(),
        2 - &i,
        Atom::num(3) * &i,
        -2 + &i,
        1 + &i,
        -1 - Atom::num(2) * &i,
    ];
    let scale = coefficients
        .iter()
        .zip(SHAPE)
        .fold(Atom::num(0), |sum, (a, power)| {
            sum + a * t.to_atom().pow(power[0] as i64) * u.to_atom().pow(power[1] as i64)
        });
    let expression = f.to_atom() / scale;
    let derivatives: Vec<_> = SHAPE
        .iter()
        .map(|power| {
            let mut derivative = expression.clone();
            for _ in 0..power[0] {
                derivative = derivative.derivative(t);
            }
            for _ in 0..power[1] {
                derivative = derivative.derivative(u);
            }
            derivative
                / Atom::num(
                    Integer::factorial(power[0] as u32) * Integer::factorial(power[1] as u32),
                )
        })
        .collect();
    Atom::evaluator_multiple(&derivatives, &[f.to_atom(), t.to_atom(), u.to_atom()])
        .build()
        .unwrap()
        .map_coeff_with_prec(&map, BITS)
}

// The solution extracts the requested Taylor coefficient from a polynomial of
// degree <=6. Tensor products of these exact weights give mixed differences.
fn stencil(order: usize) -> Vec<Rational> {
    let nodes: Vec<i64> = (-3..=3).collect();
    let matrix = Matrix::from_nested_vec(
        (0..7)
            .map(|degree| {
                nodes
                    .iter()
                    .map(|x| Rational::from(x.pow(degree)))
                    .collect()
            })
            .collect(),
        Q,
    )
    .unwrap();
    let rhs = Matrix::new_vec(
        (0..7)
            .map(|degree| Rational::from(i64::from(degree == order)))
            .collect(),
        Q,
    );
    matrix.solve(&rhs).unwrap().into_vec()
}

fn check_differences(
    evaluators: &mut Evaluators,
    base: &[C],
    first: &[C],
    second: &[C],
    label: &str,
) {
    let actual = evaluators.dual(&input(base, first, second));
    // At a massive spacelike real point the Feynman parameter denominator is
    // strictly positive. Holomorphy gives complex directional coefficients by
    // multiplying the real-axis coefficient by the corresponding seed powers.
    // Scalar complex-invariant continuation is deliberately not assumed here.
    let first_index = first.iter().position(|v| !v.is_zero()).unwrap();
    let second_index = second.iter().position(|v| !v.is_zero()).unwrap();
    let mut real_first = vec![c(0, 0); 7];
    let mut real_second = real_first.clone();
    real_first[first_index] = c(1, 0);
    real_second[second_index] = c(1, 0);
    let primal = evaluators.scalar(base);
    close(&actual[0], &primal, 240, label);
    let weights: Vec<_> = (0..=3).map(stencil).collect();
    let mut previous = vec![Float::with_val(BITS, 1); SHAPE.len()];
    for exponent in [12, 18, 24] {
        let h = Float::with_val(BITS, 2).pow(exponent).inv();
        let mut grid = Vec::new();
        for a in -3..=3 {
            for b in -3..=3 {
                let t = Complex::new(h.clone() * Float::with_val(BITS, a), Float::new(BITS));
                let u = Complex::new(h.clone() * Float::with_val(BITS, b), Float::new(BITS));
                let row: Vec<_> = base
                    .iter()
                    .zip(&real_first)
                    .zip(&real_second)
                    .map(|((z, v), w)| z + v * &t + w * &u)
                    .collect();
                grid.push(evaluators.scalar(&row));
            }
        }
        for (component, &[a, b]) in SHAPE.iter().enumerate().skip(1) {
            let mut estimate = c(0, 0);
            for j in 0..7 {
                for k in 0..7 {
                    let weight = &weights[a][j] * &weights[b][k];
                    if !weight.is_zero() {
                        estimate += &grid[7 * j + k]
                            * Complex::new(weight.to_multi_prec_float(BITS), Float::new(BITS));
                    }
                }
            }
            estimate /= Complex::new(h.pow((a + b) as u64), Float::new(BITS));
            estimate *= first[first_index].pow(a as u64) * second[second_index].pow(b as u64);
            let error = size(&(&estimate - &actual[component]));
            assert!(
                error < previous[component].clone() / Float::with_val(BITS, 100),
                "{label}: component={a},{b} h=2^-{exponent} not converged: {error:e} vs {:e}",
                previous[component]
            );
            previous[component] = error;
            if exponent == 24 {
                close(&estimate, &actual[component], 80, label);
                eprintln!(
                    "{label}: directions {first_index}/{second_index}, component={a},{b}, final error={:e}",
                    previous[component]
                );
            }
        }
    }
}

#[test]
fn massive_spacelike_c0_higher_mixed_controls() {
    let mut evaluators = Evaluators::new();
    let mut oracle = scaling_oracle();
    let mut batch_inputs = Vec::new();
    let mut batch_expected = Vec::new();
    // Källén values are -16 and +28. Positive masses keep the Feynman
    // parameter denominator away from zero; neither point is Gram singular.
    for (label, values) in POINTS {
        let base = values.map(|v| c(v, 0));
        let primal = evaluators.scalar(&base);
        let coefficients = [
            c(1, 0),
            c(0, 1),
            c(2, -1),
            c(0, 3),
            c(-2, 1),
            c(1, 1),
            c(-1, -2),
        ];
        let mut scaled: Vec<_> = base
            .iter()
            .flat_map(|z| coefficients.iter().map(move |a| z * a))
            .collect();
        // A finite massive C0 is independent of the renormalization scale.
        scaled[6 * SHAPE.len()..].fill(c(0, 0));
        scaled[6 * SHAPE.len()] = c(1, 0);
        let actual = evaluators.dual(&scaled);
        let mut expected = vec![c(0, 0); SHAPE.len()];
        oracle.evaluate(&[primal.clone(), c(0, 0), c(0, 0)], &mut expected);
        for (a, b) in actual.iter().zip(&expected) {
            close(a, b, 220, label);
        }
        batch_inputs.push(scaled);
        batch_expected.push(actual);
        let mut first = vec![c(0, 0); 7];
        let mut second = first.clone();
        first[6] = c(0, 1);
        second[6] = c(2, -1);
        let scale_only = evaluators.dual(&input(&base, &first, &second));
        close(&scale_only[0], &primal, 240, "scale-only primal");
        assert!(scale_only[1..].iter().all(|v| v.is_zero()));
    }
    let mut batch = evaluators.dual.clone().into_batch_evaluator();
    batch.evaluate_batch(&[], &mut [], 0).unwrap();
    for order in [[0, 1, 2, 3], [3, 2, 1, 0], [0, 3, 1, 2]] {
        for indices in order.chunks(3) {
            let arguments: Vec<_> = indices
                .iter()
                .flat_map(|&row| batch_inputs[row].iter().cloned())
                .collect();
            let mut output = vec![c(0, 0); indices.len() * 3 * SHAPE.len()];
            batch
                .evaluate_batch(&arguments, &mut output, indices.len())
                .unwrap();
            for (&row, values) in indices.iter().zip(output.chunks_exact(3 * SHAPE.len())) {
                for (actual, expected) in values[..SHAPE.len()].iter().zip(&batch_expected[row]) {
                    close(actual, expected, 240, "mixed-branch batch permutation");
                }
                assert!(values[SHAPE.len()..].iter().all(|v| v.is_zero()));
            }
        }
    }
    for (label, values) in POINTS {
        let base = values.map(|v| c(v, 0));
        for (first_index, second_index) in [(0, 4), (3, 2)] {
            let mut first = vec![c(0, 0); 7];
            let mut second = first.clone();
            first[first_index] = c(0, 1);
            second[second_index] = c(2, -1);
            check_differences(&mut evaluators, &base, &first, &second, label);
        }
    }
}

#[test]
fn c0_degree_five_complex_scaling() {
    let context = OneLoopExpressions::new();
    let parameters: Vec<_> = (0..7)
        .map(|j| symbolica::symbol!(format!("c0_order_five_input_{j}")).to_atom())
        .collect();
    let series = context.c0(
        [&parameters[0], &parameters[1], &parameters[2]],
        [&parameters[3], &parameters[4], &parameters[5]],
        &parameters[6],
    );
    let mut scalar = context
        .evaluator(series.coefficients(), &parameters)
        .unwrap()
        .map_coeff_with_prec(&map, BITS);
    let exact = context
        .dual_evaluator(
            series.coefficients(),
            &parameters,
            HyperDual::<f64>::new((0..=5).map(|order| vec![order]).collect()),
        )
        .unwrap();
    let mut double = exact
        .clone()
        .map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64()));
    let mut dual = exact.map_coeff_with_prec(&map, BITS);
    let t = symbolica::symbol!("c0_order_five_t");
    let value = symbolica::symbol!("c0_order_five_value");
    let phase = Atom::num(1) + Atom::num(2) * Atom::num(-1).sqrt();
    let mut derivative = value.to_atom() / (Atom::num(1) + phase * t);
    let mut coefficients = Vec::new();
    for order in 0..=5 {
        coefficients.push(&derivative / Atom::num(Integer::factorial(order)));
        derivative = derivative.derivative(t);
    }
    let mut oracle = Atom::evaluator_multiple(&coefficients, &[value.to_atom(), t.to_atom()])
        .build()
        .unwrap()
        .map_coeff_with_prec(&map, BITS);
    for (label, values) in POINTS {
        let base = values.map(|v| c(v, 0));
        let mut primal = vec![c(0, 0); 3];
        scalar.evaluate(&base, &mut primal);
        let input: Vec<_> = base
            .iter()
            .enumerate()
            .flat_map(|(j, z)| {
                [
                    z.clone(),
                    if j == 6 { c(0, 0) } else { z * c(1, 2) },
                    c(0, 0),
                    c(0, 0),
                    c(0, 0),
                    c(0, 0),
                ]
            })
            .collect();
        let mut actual = vec![c(0, 0); 18];
        dual.evaluate(&input, &mut actual);
        let mut expected = vec![c(0, 0); 6];
        oracle.evaluate(&[primal[0].clone(), c(0, 0)], &mut expected);
        for (a, b) in actual[..6].iter().zip(&expected) {
            close(a, b, 200, label);
        }
        assert!(actual[6..].iter().all(|v| v.is_zero()));
        let mut output64 = vec![Complex::new(0., 0.); 18];
        let input64: Vec<_> = input
            .iter()
            .map(|z| Complex::new(z.re.to_f64(), z.im.to_f64()))
            .collect();
        double.evaluate(&input64, &mut output64);
        for (a, b) in output64.iter().zip(&actual) {
            close(
                &Complex::new(Float::with_val(BITS, a.re), Float::with_val(BITS, a.im)),
                b,
                30,
                "degree-five f64/Float384",
            );
        }
    }
}
