//! Direct Rust tensor coefficients against original OneLOop, higher precision,
//! the exact symbolic expressions, and Ninja's rank-three bubble formulas.
use oneloop::{LaurentSeries, NativeTensorEvaluator, TensorIntegral};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::float::{Complex, DoubleFloat, Float, RealLike},
    parse,
};

struct Row {
    name: String,
    family: TensorIntegral,
    input: Vec<Complex<f64>>,
    expected: Vec<Complex<f64>>,
}

fn rows() -> Vec<Row> {
    let mut name = String::new();
    let mut rows = Vec::new();
    for line in include_str!("data/tensor_parity.txt").lines() {
        if let Some(rest) = line.strip_prefix("# name ") {
            name = rest.trim().into();
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let values: Vec<f64> = line
            .split_whitespace()
            .map(|v| v.parse().unwrap())
            .collect();
        let family = match values[0] {
            5. => TensorIntegral::An,
            6. => TensorIntegral::B11,
            7. => TensorIntegral::B111,
            _ => panic!("invalid tensor fixture kind"),
        };
        let mu_squared = Complex::new(values[1], 0.);
        let (input, rest) = match family {
            TensorIntegral::An => (
                vec![Complex::new(values[2], values[3]), mu_squared],
                &values[4..],
            ),
            _ => (
                vec![
                    Complex::new(values[2], 0.),
                    Complex::new(values[3], values[4]),
                    Complex::new(values[5], values[6]),
                    mu_squared,
                ],
                &values[7..],
            ),
        };
        assert_eq!(rest.len(), 2 * family.outputs());
        rows.push(Row {
            name: std::mem::take(&mut name),
            family,
            input,
            expected: rest
                .chunks_exact(2)
                .map(|v| Complex::new(v[0], v[1]))
                .collect(),
        });
    }
    rows
}

/// Binary64 accuracy target against 320 bits, relative to the coefficient plus
/// 1e-3 of its natural scale: the root-pair representation cancels, e.g. two
/// complex-conjugate bnlog terms of size 0.08 summing to 5e-4.
const BINARY64_TOLERANCE: f64 = 5e-14;

fn norm(z: Complex<f64>) -> f64 {
    z.re.hypot(z.im)
}

/// Natural size of a coefficient: B's are dimensionless; A's scale with m^(2k+2).
fn magnitude(row: &Row, coefficient: usize) -> f64 {
    match row.family {
        TensorIntegral::An if norm(row.input[0]) > 0. => {
            norm(row.input[0]).powi(coefficient as i32 + 1)
        }
        _ => 1.0,
    }
}

fn to_float(z: Complex<f64>, bits: u32) -> Complex<Float> {
    Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im))
}

fn evaluate_float(row: &Row, bits: u32) -> Vec<Complex<Float>> {
    let mut evaluator =
        NativeTensorEvaluator::<Float>::with_binary_precision(row.family, bits).unwrap();
    let input: Vec<_> = row.input.iter().map(|&z| to_float(z, bits)).collect();
    let mut output = vec![to_float(Complex::new(0., 0.), bits); row.family.outputs()];
    evaluator.evaluate(&input, &mut output).unwrap();
    output
}

fn evaluate_f64(row: &Row) -> Vec<Complex<f64>> {
    let mut evaluator = NativeTensorEvaluator::<f64>::new(row.family).unwrap();
    let mut output = vec![Complex::new(f64::NAN, f64::NAN); row.family.outputs()];
    evaluator.evaluate(&row.input, &mut output).unwrap();
    output
}

fn difference(a: &Complex<Float>, b: &Complex<Float>) -> f64 {
    let re = (a.re.clone() - &b.re).to_f64();
    let im = (a.im.clone() - &b.im).to_f64();
    re.hypot(im)
}

fn float_norm(z: &Complex<Float>) -> f64 {
    z.re.to_f64().hypot(z.im.to_f64())
}

#[test]
fn binary64_matches_original_oneloop() {
    let rows = rows();
    assert_eq!(rows.len(), 366);
    let mut worst = (0.0_f64, String::new());
    let mut failures = Vec::new();
    let mut attributed = Vec::new();
    for row in &rows {
        let actual = evaluate_f64(row);
        let reference = evaluate_float(row, 320);
        for (index, (value, expected)) in actual.iter().zip(&row.expected).enumerate() {
            let coefficient = index / 3;
            let scale = magnitude(row, coefficient);
            let error = norm(*value - *expected) / (norm(*expected) + 1e-3 * scale);
            // Fortran itself loses a few digits in places (cancellation in
            // its direct bnlog formula, or its finite regulators); attribute a
            // disagreement to Fortran only when this port matches 320 bits.
            let accurate = difference(&to_float(*value, 320), &reference[index])
                <= BINARY64_TOLERANCE * (float_norm(&reference[index]) + 1e-3 * scale);
            if error > 1e-12 && accurate && error < 1e-9 {
                attributed.push(format!("{} output {index} ({error:.1e})", row.name));
            }
            if error > 1e-12 && !(accurate && error < 1e-9) {
                failures.push(format!(
                    "{} output {index}: {value:?} vs Fortran {expected:?} (relative {error:e})",
                    row.name
                ));
            }
            if error > worst.0 {
                worst = (error, format!("{} output {index}", row.name));
            }
        }
    }
    eprintln!(
        "worst relative difference to Fortran: {:e} at {}",
        worst.0, worst.1
    );
    eprintln!(
        "{} outputs differ from Fortran by more than 1e-12 where this port matches 320 bits: {}",
        attributed.len(),
        attributed.join(", ")
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// Change of each output at 320 bits when p² moves by a relative `2^-bits`:
/// threshold rows are ill-conditioned there (`Im B ~ sqrt(p² - threshold)`),
/// so no binary64 algorithm can do better for an input rounded to binary64.
fn momentum_sensitivity(row: &Row, reference: &[Complex<Float>], bits: i32) -> Vec<f64> {
    if row.family == TensorIntegral::An || row.input[0].re == 0. {
        return vec![0.; reference.len()];
    }
    let mut sensitivity = vec![0.0_f64; reference.len()];
    for sign in [-1.0, 1.0] {
        let mut shifted = Row {
            name: String::new(),
            family: row.family,
            input: row.input.clone(),
            expected: vec![],
        };
        let mut evaluator =
            NativeTensorEvaluator::<Float>::with_binary_precision(row.family, 320).unwrap();
        let mut input: Vec<_> = shifted.input.drain(..).map(|z| to_float(z, 320)).collect();
        let factor = Float::with_val(320, 1) + Float::with_val(320, sign * 2f64.powi(-bits));
        input[0].re = input[0].re.clone() * factor;
        let mut output = vec![to_float(Complex::new(0., 0.), 320); row.family.outputs()];
        evaluator.evaluate(&input, &mut output).unwrap();
        for (value, (shifted, original)) in sensitivity.iter_mut().zip(output.iter().zip(reference))
        {
            *value = value.max(difference(shifted, original));
        }
    }
    sensitivity
}

#[test]
fn binary64_and_double_float_match_higher_precision() {
    let mut worst = 0.0_f64;
    for row in rows() {
        let reference = evaluate_float(&row, 320);
        let binary64_sensitivity = momentum_sensitivity(&row, &reference, 53);
        let double_sensitivity = momentum_sensitivity(&row, &reference, 106);
        let actual = evaluate_f64(&row);
        let mut double = NativeTensorEvaluator::<DoubleFloat>::new(row.family).unwrap();
        let input: Vec<_> = row
            .input
            .iter()
            .map(|z| Complex::new(DoubleFloat::from(z.re), DoubleFloat::from(z.im)))
            .collect();
        let zero = Complex::new(DoubleFloat::from(0.), DoubleFloat::from(0.));
        let mut double_output = vec![zero; row.family.outputs()];
        double.evaluate(&input, &mut double_output).unwrap();
        for index in 0..row.family.outputs() {
            let scale = magnitude(&row, index / 3);
            let size = float_norm(&reference[index]) + 1e-3 * scale;
            let error = (difference(&to_float(actual[index], 320), &reference[index])
                - 4. * binary64_sensitivity[index])
                .max(0.)
                / size;
            worst = worst.max(error);
            assert!(
                error < BINARY64_TOLERANCE,
                "{} output {index}: f64 relative error {error:e}",
                row.name
            );
            let double_value = Complex::new(
                Float::with_val(320, double_output[index].re.into_inner().hi())
                    + Float::with_val(320, double_output[index].re.into_inner().lo()),
                Float::with_val(320, double_output[index].im.into_inner().hi())
                    + Float::with_val(320, double_output[index].im.into_inner().lo()),
            );
            let error = (difference(&double_value, &reference[index])
                - 4. * double_sensitivity[index])
                .max(0.)
                / size;
            assert!(
                error < 1e-28,
                "{} output {index}: DoubleFloat relative error {error:e}",
                row.name
            );
        }
    }
    eprintln!("worst binary64 relative error against 320 bits: {worst:e}");
}

#[test]
fn arbitrary_precision_converges() {
    for row in rows().iter().step_by(7) {
        let coarse = evaluate_float(row, 200);
        let fine = evaluate_float(row, 400);
        for index in 0..row.family.outputs() {
            let size = float_norm(&fine[index]) + 1e-3 * magnitude(row, index / 3);
            let error = difference(&coarse[index], &fine[index]) / size;
            // Symbolica's Float tracks cancellation of its inputs: near the
            // on-shell point (p²-m² = 1e-12) the 200-bit result keeps ~160 bits.
            assert!(error < 1e-45, "{} output {index}: {error:e}", row.name);
        }
    }
}

fn symbolic_coefficients(family: TensorIntegral, parameters: &[Atom]) -> Vec<Atom> {
    let series: Vec<LaurentSeries> = match family {
        TensorIntegral::An => oneloop::an(4, &parameters[0], &parameters[1]).unwrap(),
        _ => {
            let [p, m0, m1, mu] = [0, 1, 2, 3].map(|i| &parameters[i]);
            let mut series = vec![
                oneloop::b0(p, m0, m1, mu),
                oneloop::b1(p, m0, m1, mu),
                oneloop::b11(p, m0, m1, mu),
            ];
            if family == TensorIntegral::B111 {
                series.push(oneloop::b111(p, m0, m1, mu));
            }
            series
        }
    };
    series
        .iter()
        .flat_map(|s| s.coefficients().iter().cloned())
        .collect()
}

#[test]
fn arbitrary_precision_matches_symbolic_expressions() {
    // The Passarino-Veltman expressions cancel like (m²/p²)^rank; evaluate
    // them with enough guard bits to keep 1e-60 at the smallest momenta.
    const BITS: u32 = 256;
    const SYMBOLIC_BITS: u32 = 640;
    for family in [TensorIntegral::An, TensorIntegral::B111] {
        let parameters: Vec<Atom> = (0..family.arity())
            .map(|i| parse!(format!("tensor_parameter_{}_{i}", family.name())))
            .collect();
        let coefficients = symbolic_coefficients(family, &parameters);
        let views: Vec<_> = coefficients.iter().map(Atom::as_view).collect();
        let mut evaluator = Atom::evaluator_multiple(&views, &parameters)
            .build()
            .unwrap()
            .map_coeff_with_prec(
                &|v| {
                    Complex::new(
                        v.re.to_multi_prec_float(SYMBOLIC_BITS),
                        v.im.to_multi_prec_float(SYMBOLIC_BITS),
                    )
                },
                SYMBOLIC_BITS,
            );
        let mut checked = 0;
        for row in rows().iter().filter(|row| row.family == family) {
            let input: Vec<_> = row
                .input
                .iter()
                .map(|&z| to_float(z, SYMBOLIC_BITS))
                .collect();
            let mut expected =
                vec![to_float(Complex::new(0., 0.), SYMBOLIC_BITS); family.outputs()];
            evaluator.evaluate(&input, &mut expected);
            let actual = evaluate_float(row, BITS);
            for index in 0..family.outputs() {
                // The symbolic zero-momentum B1/B11/B111 anchors put real
                // opposite-sign masses on the wrong lip; see the analytic test.
                if row.name.contains("p0_opposite_signs") && index >= 3 {
                    continue;
                }
                let size = float_norm(&expected[index]) + 1e-3 * magnitude(row, index / 3);
                let error = difference(&actual[index], &expected[index]) / size;
                assert!(
                    error < 1e-60,
                    "{} output {index}: native {} vs symbolic {} ({error:e})",
                    row.name,
                    actual[index],
                    expected[index]
                );
            }
            checked += 1;
        }
        assert!(
            checked
                >= if family == TensorIntegral::An {
                    18
                } else {
                    100
                }
        );
    }
}

/// At p² = 0 with real masses, `B_n = (-1)^(n+1) int_0^1 x^n log(D(x) - i0) dx`
/// with `D(x) = (1-x) m0² + x m1²`, so `Im B_n = (-1)^n pi int_(D<0) x^n dx`.
#[test]
fn zero_momentum_imaginary_parts_match_feynman_parameters() {
    let mut checked = 0;
    for row in rows()
        .iter()
        .filter(|row| row.family == TensorIntegral::B111)
    {
        let (p, m0, m1) = (row.input[0], row.input[1], row.input[2]);
        if p.re != 0. || m0.im != 0. || m1.im != 0. || (m0.re == 0. && m1.re == 0.) {
            continue;
        }
        // D(x) < 0 on x in (a, b).
        let (a, b) = match (m0.re < 0., m1.re < 0.) {
            (false, false) => (0., 0.),
            (true, true) => (0., 1.),
            _ => {
                let root = m0.re / (m0.re - m1.re);
                if m1.re < 0. { (root, 1.) } else { (0., root) }
            }
        };
        let actual = evaluate_f64(row);
        for n in 0..4 {
            let integral = (b.powi(n as i32 + 1) - a.powi(n as i32 + 1)) / (n as f64 + 1.);
            let expected = if n % 2 == 0 { 1. } else { -1. } * std::f64::consts::PI * integral;
            let value = actual[3 * n].im;
            assert!(
                (value - expected).abs() < 1e-14,
                "{} B rank {n}: Im {value} vs {expected}",
                row.name
            );
        }
        checked += 1;
    }
    assert!(checked >= 20, "{checked}");
}

/// Ninja's IntegralLibrary::getRank3BubbleIntegral, from B0/B1/B11 and A0.
fn ninja_b111(row: &Row, bits: u32) -> [Complex<Float>; 3] {
    let float = |x: f64| Float::with_val(bits, x);
    let complex = |z: Complex<f64>| to_float(z, bits);
    let tensors = evaluate_float(row, bits);
    let tadpole = |mass: Complex<f64>| {
        let input = Row {
            name: String::new(),
            family: TensorIntegral::An,
            input: vec![mass, row.input[3]],
            expected: vec![],
        };
        evaluate_float(&input, bits)
    };
    let p2 = complex(row.input[0]);
    let (m1, m2) = (complex(row.input[1]), complex(row.input[2]));
    let (b0, b1, b11) = (&tensors[0..3], &tensors[3..6], &tensors[6..9]);
    if row.input[0].re.abs() > 1e-8 {
        let a1 = tadpole(row.input[1]);
        let a2 = tadpole(row.input[2]);
        let six = complex(Complex::new(6., 0.));
        let three = complex(Complex::new(3., 0.));
        let quarter = float(0.25);
        let rational_b00 = (p2.clone() - three * (m1.clone() + &m2)) / six;
        let b00: Vec<_> = (0..3)
            .map(|i| {
                let mut i2q2 = a2[i].clone() + m1.clone() * &b0[i];
                if i == 0 {
                    i2q2 -= &rational_b00;
                }
                (i2q2 - p2.clone() * &b11[i]) * &quarter
            })
            .collect();
        let rational = -(p2.clone() * &p2)
            * (p2.clone() - m1.clone() * float(2.) - m2.clone() * float(4.))
            / float(3.)
            + p2.clone() * (m1.clone() * &m1 - m2.clone() * &m2) * float(0.5);
        let c1 = -(p2.clone() * &p2 * &m1) * float(4.);
        let c2 = (m2.clone() - &m1 - &p2) * float(4.);
        let denominator = p2.clone() * &p2 * &p2 * float(4.);
        std::array::from_fn(|i| {
            let mut value = c1.clone() * &b1[i]
                + c2.clone() * &p2 * (b00[i].clone() + p2.clone() * &b11[i])
                + p2.clone() * (m1.clone() * &a1[i] - m2.clone() * &a2[i]);
            if i == 0 {
                value += &rational;
            }
            value / &denominator
        })
    } else {
        let bubble = |mass: Complex<f64>| {
            let input = Row {
                name: String::new(),
                family: TensorIntegral::B11,
                input: vec![Complex::new(0., 0.), mass, mass, row.input[3]],
                expected: vec![],
            };
            evaluate_float(&input, bits)
        };
        let z11 = bubble(row.input[1]);
        if row.input[1] != row.input[2] {
            let z22 = bubble(row.input[2]);
            let d = m1.clone() - &m2;
            let d3 = d.clone() * &d * &d;
            let inv4 = d3.clone() * &d;
            let m1q4 = m1.clone() * &m1 * &m1 * &m1;
            let coefficient_1 = -(m1q4.clone() / &inv4) * float(0.25);
            let coefficient_2 = (m1q4 / &inv4 - float(1.)) * float(0.25);
            let rational = (m1.clone() * &m1 * &m1 * float(-25.)
                + m1.clone() * &m1 * &m2 * float(23.)
                - m1.clone() * &m2 * &m2 * float(13.)
                + m2.clone() * &m2 * &m2 * float(3.))
                / (d3 * float(48.));
            std::array::from_fn(|i| {
                let mut value = coefficient_1.clone() * &z11[i] + coefficient_2.clone() * &z22[i];
                if i == 0 {
                    value += &rational;
                }
                value
            })
        } else {
            std::array::from_fn(|i| -(z11[i].clone() * float(0.25)))
        }
    }
}

#[test]
fn rank_three_bubble_matches_ninja_formulas() {
    const BITS: u32 = 400;
    let mut checked = 0;
    for row in rows()
        .iter()
        .filter(|row| row.family == TensorIntegral::B111)
    {
        // Ninja replaces 0 < |p²| <= 1e-8 by its p²=0 form, which has no
        // scaleless limit; compare where it is exact.
        let p = row.input[0].re;
        let scaleless = row.input[1..3].iter().all(|m| m.re == 0. && m.im == 0.);
        if (p != 0. && p.abs() <= 1e-8) || (p == 0. && scaleless) {
            continue;
        }
        let native = evaluate_float(row, BITS);
        let ninja = ninja_b111(row, BITS);
        for (laurent, value) in ninja.iter().enumerate() {
            let error = difference(&native[9 + laurent], value) / (float_norm(value) + 1e-3);
            // PV division by p⁶ cancels like (m²/p²)^3 even at 400 bits.
            assert!(
                error < 1e-60,
                "{} order {laurent}: native {} vs Ninja {value} ({error:e})",
                row.name,
                native[9 + laurent]
            );
        }
        checked += 1;
    }
    assert!(checked >= 100, "{checked}");
}

#[test]
fn batches_validate_rows_and_match_single_rows() {
    let rows = rows();
    for family in [
        TensorIntegral::An,
        TensorIntegral::B11,
        TensorIntegral::B111,
    ] {
        let selected: Vec<_> = rows.iter().filter(|row| row.family == family).collect();
        let input: Vec<_> = selected
            .iter()
            .flat_map(|row| row.input.iter().copied())
            .collect();
        let mut evaluator = NativeTensorEvaluator::<f64>::new(family).unwrap();
        let mut output = vec![Complex::new(f64::NAN, f64::NAN); selected.len() * family.outputs()];
        evaluator
            .evaluate_batch(&input, &mut output, selected.len())
            .unwrap();
        for (row, values) in selected.iter().zip(output.chunks_exact(family.outputs())) {
            assert_eq!(values, evaluate_f64(row).as_slice(), "{}", row.name);
        }
        evaluator.evaluate_batch(&[], &mut [], 0).unwrap();
        // One invalid row rejects the whole batch without writing any output.
        let mut invalid = input.clone();
        let last = invalid.len() - 1;
        invalid[last] = Complex::new(-1.0, 0.0);
        let mut untouched = vec![Complex::new(7.0, 7.0); output.len()];
        assert!(
            evaluator
                .evaluate_batch(&invalid, &mut untouched, selected.len())
                .is_err()
        );
        assert!(untouched.iter().all(|z| *z == Complex::new(7.0, 7.0)));
        assert!(
            evaluator
                .evaluate_batch(&input, &mut output[1..], selected.len())
                .is_err()
        );
        let mut positive_width = selected[0].input.clone();
        positive_width[family.arity() - 2].im = 0.5;
        assert!(
            evaluator
                .evaluate(&positive_width, &mut output[..family.outputs()])
                .is_err()
        );
    }
    assert!(NativeTensorEvaluator::<f64>::with_binary_precision(TensorIntegral::B11, 64).is_err());
    assert!(NativeTensorEvaluator::<Float>::with_binary_precision(TensorIntegral::B11, 0).is_err());
}
