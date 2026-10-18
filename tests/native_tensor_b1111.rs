//! The rank-four bubble family `TensorIntegral::B1111` (`[B0, B1, B11, B111,
//! B1111]`, OneLOop's `bub1111` with the corrected mass-swap sign) on every
//! `NativeFloat` route: the binary64 kernels, the generic code in binary64
//! (`set_binary64(false)`), `DoubleFloat` and `Float`.
//!
//! - Against the exact symbolic `b0`..`b1111` expressions at 320 bits, or at
//!   1280/5120 bits where their Passarino-Veltman forms cancel beyond that.
//! - Against direct Feynman-parameter integrals (`tests/data/tensor_feynman.txt`,
//!   `tests/generate_tensor_feynman_fixtures.py`), per sector, including every
//!   branch of the mass swap and a dyadic p² ladder 2^-4..2^-24.
//! - The mass-swap identity `B_n(p; m0, m1) = (-1)^n sum_k C(n,k) B_k(p; m1, m0)`.
//!   OneLOop 3.4 and 3.6 negate the rank-four sum in all three swap branches;
//!   that sign violates the identity by `2 |B1111|`.
//! - `B0`..`B111` of the family are bitwise those of `B111` on every route.
//!
//! Errors of `B0`..`B111` are scaled as in `native_tensor.rs`, by the
//! coefficient plus 1e-3. `B1111` is scaled by the largest coefficient of its
//! row plus 1e-3: after a mass swap, `B1111 = B0' + 4B1' + 6B11' + 4B111' +
//! B1111'` cancels to a value up to about ten times smaller than its terms,
//! so its rounding error is a few ulps of `max_k |B_k|`, not of `|B1111|`.
use oneloop::{LaurentSeries, NativeTensorEvaluator, TensorIntegral};
use std::collections::BTreeMap;
use symbolica::{
    atom::{Atom, AtomCore},
    domains::float::{Complex, DoubleFloat, Float, RealLike},
    parse,
};

type C = Complex<f64>;
type D = Complex<DoubleFloat>;

const FAMILY: TensorIntegral = TensorIntegral::B1111;
/// Binary64 accuracy target against 320 bits, relative to the coefficient plus
/// 1e-3 of its natural scale, as in `native_tensor.rs`.
const BINARY64_TOLERANCE: f64 = 5e-14;
/// `DoubleFloat` target, as in `native_tensor.rs`.
const DOUBLE_TOLERANCE: f64 = 1e-28;
/// Native `Float` at 256 bits against the expressions, as for `B111`.
const FLOAT_TOLERANCE: f64 = 1e-60;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn uniform(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn log_uniform(&mut self, low: f64, high: f64) -> f64 {
        (low.ln() + (high.ln() - low.ln()) * self.uniform()).exp()
    }
}

/// A squared mass: massless, real (positive or negative) or complex.
fn mass(rng: &mut Rng, scale: f64) -> C {
    match rng.below(8) {
        0 | 1 => C::new(0., 0.),
        2 => C::new(-rng.log_uniform(1e-3, 1e3) * scale, 0.),
        3 => C::new(
            rng.log_uniform(1e-3, 1e3) * scale,
            -rng.log_uniform(1e-4, 1.) * scale,
        ),
        4 => C::new(173.0f64.powi(2), 0.),
        _ => C::new(rng.log_uniform(1e-4, 1e4) * scale, 0.),
    }
}

/// `[p², m0², m1², mu²]` over every sector (as `native_tensor_binary64.rs`).
fn bubble_point(rng: &mut Rng) -> Vec<C> {
    let scale = rng.log_uniform(1e-3, 1e5);
    let m0 = mass(rng, scale);
    let mut m1 = mass(rng, scale);
    if rng.below(6) == 0 {
        m1 = m0;
    }
    let real = |z: C| z.re;
    let p = match rng.below(12) {
        0 => 0.,
        1 => real(m0),
        2 => real(m1),
        3 | 4 => {
            let (a, b) = (real(m0).abs().sqrt(), real(m1).abs().sqrt());
            let sign = if rng.below(2) == 0 { 1. } else { -1. };
            let threshold = (a + sign * b).powi(2);
            threshold * (1. + (rng.uniform() - 0.5) * 10f64.powi(-(rng.below(12) as i32)))
        }
        5 => (rng.uniform() - 0.5) * scale * 10f64.powi(-(rng.below(16) as i32)),
        6 => -rng.log_uniform(1e-4, 1e6) * scale,
        _ => rng.log_uniform(1e-4, 1e6) * scale,
    };
    let (m0, m1) = if rng.below(2) == 0 {
        (m0, m1)
    } else {
        (m1, m0)
    };
    let mu = C::new(rng.log_uniform(1e-2, 1e2) * scale, 0.);
    vec![C::new(p, 0.), m0, m1, mu]
}

fn random_points(seed: u64, count: usize) -> Vec<(String, Vec<C>)> {
    let mut rng = Rng(seed);
    (0..count)
        .map(|i| (format!("random_{i}"), bubble_point(&mut rng)))
        .collect()
}

fn norm(z: C) -> f64 {
    z.re.hypot(z.im)
}

fn to_float(z: C, bits: u32) -> Complex<Float> {
    Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im))
}

fn double_to_float(z: &D) -> Complex<Float> {
    let part = |x: &DoubleFloat| {
        let x = (*x).into_inner();
        Float::with_val(320, x.hi()) + Float::with_val(320, x.lo())
    };
    Complex::new(part(&z.re), part(&z.im))
}

fn difference(a: &Complex<Float>, b: &Complex<Float>) -> f64 {
    let re = (a.re.clone() - &b.re).to_f64();
    let im = (a.im.clone() - &b.im).to_f64();
    re.hypot(im)
}

fn float_norm(z: &Complex<Float>) -> f64 {
    z.re.to_f64().hypot(z.im.to_f64())
}

fn evaluate_f64(family: TensorIntegral, input: &[C], binary64: bool) -> Vec<C> {
    let mut evaluator = NativeTensorEvaluator::<f64>::new(family).unwrap();
    evaluator.set_binary64(binary64).unwrap();
    let mut output = vec![C::new(f64::NAN, f64::NAN); family.outputs()];
    evaluator.evaluate(input, &mut output).unwrap();
    output
}

fn evaluate_double(family: TensorIntegral, input: &[C]) -> Vec<D> {
    let mut evaluator = NativeTensorEvaluator::<DoubleFloat>::new(family).unwrap();
    let input: Vec<_> = input
        .iter()
        .map(|z| Complex::new(DoubleFloat::from(z.re), DoubleFloat::from(z.im)))
        .collect();
    let zero = Complex::new(DoubleFloat::from(0.), DoubleFloat::from(0.));
    let mut output = vec![zero; family.outputs()];
    evaluator.evaluate(&input, &mut output).unwrap();
    output
}

fn evaluate_float(family: TensorIntegral, input: &[C], bits: u32) -> Vec<Complex<Float>> {
    let mut evaluator =
        NativeTensorEvaluator::<Float>::with_binary_precision(family, bits).unwrap();
    let input: Vec<_> = input.iter().map(|&z| to_float(z, bits)).collect();
    let mut output = vec![to_float(C::new(0., 0.), bits); family.outputs()];
    evaluator.evaluate(&input, &mut output).unwrap();
    output
}

/// Largest change of each output at 320 bits when p² moves by a relative
/// `2^-bits`: threshold rows are ill-conditioned (`Im B ~ sqrt(p² - threshold)`),
/// so no algorithm at that precision can do better for its rounded input.
fn momentum_sensitivity(input: &[C], reference: &[Complex<Float>], bits: i32) -> Vec<f64> {
    if input[0].re == 0. {
        return vec![0.; reference.len()];
    }
    let mut sensitivity = vec![0.0_f64; reference.len()];
    let mut evaluator = NativeTensorEvaluator::<Float>::with_binary_precision(FAMILY, 320).unwrap();
    for sign in [-1.0, 1.0] {
        let mut shifted: Vec<_> = input.iter().map(|&z| to_float(z, 320)).collect();
        let factor = Float::with_val(320, 1) + Float::with_val(320, sign * 2f64.powi(-bits));
        shifted[0].re = shifted[0].re.clone() * factor;
        let mut output = vec![to_float(C::new(0., 0.), 320); reference.len()];
        evaluator.evaluate(&shifted, &mut output).unwrap();
        for (value, (shifted, original)) in sensitivity.iter_mut().zip(output.iter().zip(reference))
        {
            *value = value.max(difference(shifted, original));
        }
    }
    sensitivity
}

struct FeynmanRow {
    name: String,
    input: Vec<C>,
    /// Finite parts of B0, B1, B11, B111, B1111 at 40 digits.
    expected: Vec<Complex<Float>>,
}

fn feynman_rows() -> Vec<FeynmanRow> {
    let mut name = String::new();
    let mut rows = Vec::new();
    for line in include_str!("data/tensor_feynman.txt").lines() {
        if let Some(rest) = line.strip_prefix("# name ") {
            name = rest.trim().into();
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let words: Vec<&str> = line.split_whitespace().collect();
        assert_eq!(words.len(), 16, "{line}");
        let value = |i: usize| words[i].parse::<f64>().unwrap();
        let input = vec![
            C::new(value(0), 0.),
            C::new(value(1), value(2)),
            C::new(value(3), value(4)),
            C::new(value(5), 0.),
        ];
        let parse_float = |s: &str| Float::parse(s, Some(320)).unwrap();
        let expected = words[6..]
            .chunks_exact(2)
            .map(|pair| Complex::new(parse_float(pair[0]), parse_float(pair[1])))
            .collect();
        rows.push(FeynmanRow {
            name: std::mem::take(&mut name),
            input,
            expected,
        });
    }
    rows
}

/// The bubble rows of the OneLOop fixtures (B11 and B111 kinematics).
fn parity_points() -> Vec<(String, Vec<C>)> {
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
        if values[0] == 6. || values[0] == 7. {
            rows.push((
                std::mem::take(&mut name),
                vec![
                    C::new(values[2], 0.),
                    C::new(values[3], values[4]),
                    C::new(values[5], values[6]),
                    C::new(values[1], 0.),
                ],
            ));
        }
    }
    rows
}

fn all_points(random: usize) -> Vec<(String, Vec<C>)> {
    let mut points: Vec<_> = feynman_rows()
        .into_iter()
        .map(|row| (row.name, row.input))
        .collect();
    points.extend(parity_points());
    points.extend(random_points(0x3c6e_f372_fe94_f82b, random));
    points
}

/// The branch of OneLOop's bubble algorithm a row takes, and whether the
/// result is reconstructed from the swapped mass order.
fn sector(input: &[C]) -> (&'static str, bool) {
    let size = |z: C| z.re.abs() + z.im.abs();
    let is_zero = |z: C| z.re == 0. && z.im == 0.;
    let (p, m0, m1) = (input[0].re, input[1], input[2]);
    if p == 0. && is_zero(m0) && is_zero(m1) {
        return ("scaleless", false);
    }
    let switch = size(m1) < size(m0);
    let (a0, a1) = if switch { (m1, m0) } else { (m0, m1) };
    if p == 0. {
        let name = if a0 == a1 {
            "p0_equal"
        } else if is_zero(a0) {
            "p0_massless"
        } else {
            "p0_unequal"
        };
        return (name, !switch);
    }
    let name = if is_zero(a0) {
        if a1.im == 0. && a1.re == p {
            "massless_on_shell"
        } else {
            "massless"
        }
    } else if a0.im != 0. || a1.im != 0. {
        "complex_masses"
    } else {
        let linear = a1.re - a0.re - p;
        if linear * linear - 4. * p * a0.re < 0. {
            "conjugate_roots"
        } else {
            "real_roots"
        }
    };
    (name, switch)
}

#[test]
fn lower_ranks_are_those_of_b111_bitwise() {
    let points = all_points(2000);
    for (index, (name, input)) in points.iter().enumerate() {
        for binary64 in [true, false] {
            let rank_four = evaluate_f64(FAMILY, input, binary64);
            assert_eq!(
                evaluate_f64(TensorIntegral::B111, input, binary64),
                rank_four[..12],
                "{name} binary64 kernels {binary64}"
            );
            assert_eq!(
                (rank_four[13], rank_four[14]),
                if sector(input).0 == "scaleless" {
                    (C::new(0., 0.), C::new(0., 0.))
                } else {
                    (C::new(0.2, 0.), C::new(0., 0.))
                },
                "{name}: B1111 poles"
            );
        }
        let rank_four = evaluate_double(FAMILY, input);
        let rank_three = evaluate_double(TensorIntegral::B111, input);
        let bits = |z: &D| {
            [&z.re, &z.im].map(|x| {
                let x = (*x).into_inner();
                (x.hi().to_bits(), x.lo().to_bits())
            })
        };
        for (a, b) in rank_three.iter().zip(&rank_four) {
            assert_eq!(bits(a), bits(b), "{name} DoubleFloat");
        }
        if index % 4 == 0 {
            let rank_four = evaluate_float(FAMILY, input, 256);
            let rank_three = evaluate_float(TensorIntegral::B111, input, 256);
            for (a, b) in rank_three.iter().zip(&rank_four) {
                assert!(a.re == b.re && a.im == b.im, "{name} Float: {a} vs {b}");
            }
        }
    }
}

fn symbolic_coefficients(parameters: &[Atom]) -> Vec<Atom> {
    let [p, m0, m1, mu] = [0, 1, 2, 3].map(|i| &parameters[i]);
    let series: [LaurentSeries; 5] = [
        oneloop::b0(p, m0, m1, mu),
        oneloop::b1(p, m0, m1, mu),
        oneloop::b11(p, m0, m1, mu),
        oneloop::b111(p, m0, m1, mu),
        oneloop::b1111(p, m0, m1, mu),
    ];
    series
        .iter()
        .flat_map(|s| s.coefficients().iter().cloned())
        .collect()
}

type Symbolic = Box<dyn FnMut(&[C]) -> Vec<Complex<Float>>>;

/// The exact expressions `b0`..`b1111`, evaluated at `bits`.
fn symbolic(bits: u32) -> Symbolic {
    let parameters: Vec<Atom> = (0..4)
        .map(|i| parse!(format!("b1111_parameter_{i}")))
        .collect();
    let coefficients = symbolic_coefficients(&parameters);
    let views: Vec<_> = coefficients.iter().map(Atom::as_view).collect();
    let mut evaluator = Atom::evaluator_multiple(&views, &parameters)
        .build()
        .unwrap()
        .map_coeff_with_prec(
            &|v| {
                Complex::new(
                    v.re.to_multi_prec_float(bits),
                    v.im.to_multi_prec_float(bits),
                )
            },
            bits,
        );
    Box::new(move |input: &[C]| {
        let input: Vec<_> = input.iter().map(|&z| to_float(z, bits)).collect();
        let mut output = vec![to_float(C::new(0., 0.), bits); FAMILY.outputs()];
        evaluator.evaluate(&input, &mut output);
        output
    })
}

type SymbolicF64 = Box<dyn FnMut(&[C]) -> Vec<C>>;

/// The exact expressions in binary64 (Symbolica's evaluator over `f64`).
fn symbolic_f64() -> SymbolicF64 {
    let parameters: Vec<Atom> = (0..4)
        .map(|i| parse!(format!("b1111_f64_parameter_{i}")))
        .collect();
    let coefficients = symbolic_coefficients(&parameters);
    let views: Vec<_> = coefficients.iter().map(Atom::as_view).collect();
    let mut evaluator = Atom::evaluator_multiple(&views, &parameters)
        .build()
        .unwrap()
        .map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64()));
    Box::new(move |input: &[C]| {
        let mut output = vec![C::new(0., 0.); FAMILY.outputs()];
        evaluator.evaluate(input, &mut output);
        output
    })
}

/// The error scale of each output of `reference` (see the module notes).
fn scales(reference: &[Complex<Float>]) -> Vec<f64> {
    (0..FAMILY.outputs())
        .map(|index| {
            let (rank, order) = (index / 3, index % 3);
            if rank < 4 {
                float_norm(&reference[index]) + 1e-3
            } else {
                (0..5)
                    .map(|k| float_norm(&reference[3 * k + order]))
                    .fold(0.0, f64::max)
                    + 1e-3
            }
        })
        .collect()
}

/// Whether a row has an input component beyond `1e±100` (`DoubleFloat`'s
/// logarithm keeps about 90 bits there, already for `B0`).
fn extreme(input: &[C]) -> bool {
    input
        .iter()
        .flat_map(|z| [z.re, z.im])
        .any(|x| x != 0. && !(1e-100..=1e100).contains(&x.abs()))
}

fn finite(z: &Complex<Float>) -> bool {
    z.re.to_f64().is_finite() && z.im.to_f64().is_finite()
}

/// Every row in every route against the exact expressions. The expressions'
/// Passarino-Veltman forms divide by `2p²` per rank above one and cancel like
/// `(m²/p²)^4` at rank four; each row uses the lowest of 320, 1280 and 5120
/// bits that agrees with the next to `1e-70`, and is skipped only if none does
/// (then `|p²|` is below `1e-200` of the masses, where no fixed precision of
/// these forms is meaningful).
#[test]
fn every_route_matches_the_exact_expressions() {
    let points = all_points(1500);
    let mut at_320 = symbolic(320);
    let mut at_1280 = symbolic(1280);
    let mut at_5120 = symbolic(5120);
    let (mut worst_f64, mut worst_generic, mut worst_double, mut worst_float) =
        (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    let (mut beyond_320, mut skipped) = (0, Vec::new());
    let (mut limited, mut failures, mut nonfinite) = (Vec::new(), Vec::new(), Vec::new());
    for (name, input) in &points {
        let coarse = at_320(input);
        let fine = at_1280(input);
        let size = scales(&fine);
        let agrees = |a: &[Complex<Float>], b: &[Complex<Float>]| {
            a.iter()
                .zip(b)
                .zip(&size)
                .all(|((a, b), size)| finite(a) && finite(b) && difference(a, b) <= 1e-70 * size)
        };
        let expression = if agrees(&coarse, &fine) {
            coarse
        } else {
            let finest = at_5120(input);
            if agrees(&fine, &finest) {
                beyond_320 += 1;
                fine
            } else {
                let scale = [input[1], input[2], input[3]]
                    .iter()
                    .map(|&z| norm(z))
                    .fold(0.0, f64::max);
                if input[0].re.abs() <= 1e-200 * scale {
                    skipped.push(name.clone());
                } else {
                    nonfinite.push(format!(
                        "{name} {input:?}: the expressions do not converge with precision"
                    ));
                }
                continue;
            }
        };
        let reference = evaluate_float(FAMILY, input, 320);
        let sensitivity_53 = momentum_sensitivity(input, &reference, 53);
        let sensitivity_106 = momentum_sensitivity(input, &reference, 106);
        let fast = evaluate_f64(FAMILY, input, true);
        let generic = evaluate_f64(FAMILY, input, false);
        let double = evaluate_double(FAMILY, input);
        let float = evaluate_float(FAMILY, input, 256);
        let double_errors: Vec<f64> = (0..FAMILY.outputs())
            .map(|index| {
                (difference(&double_to_float(&double[index]), &expression[index])
                    - 4. * sensitivity_106[index])
                    .max(0.)
                    / size[index]
            })
            .collect();
        let lower_double = double_errors[..12].iter().copied().fold(0.0, f64::max);
        for index in 0..FAMILY.outputs() {
            let expected = &expression[index];
            let size = size[index];
            let scaled = |value: &Complex<Float>, sensitivity: f64| {
                (difference(value, expected) - 4. * sensitivity).max(0.) / size
            };
            let fast_error = scaled(&to_float(fast[index], 320), sensitivity_53[index]);
            let generic_error = scaled(&to_float(generic[index], 320), sensitivity_53[index]);
            let double_error = double_errors[index];
            let float_error = difference(&float[index], expected) / size;
            worst_generic = worst_generic.max(generic_error);
            worst_float = worst_float.max(float_error);
            let message = format!(
                "{name} {input:?} output {index}: expression {expected}, binary64 kernels {:?} \
                 ({fast_error:.1e}), generic binary64 {:?} ({generic_error:.1e}), \
                 DoubleFloat {double_error:.1e}, Float {float_error:.1e}",
                fast[index], generic[index]
            );
            // As in native_tensor_binary64.rs: the binary64 target, or the
            // generic code's own accuracy where the algorithm loses digits.
            if fast_error < BINARY64_TOLERANCE {
                worst_f64 = worst_f64.max(fast_error);
            } else if fast_error <= 4. * generic_error {
                limited.push(message.clone());
            } else {
                failures.push(message.clone());
            }
            // DoubleFloat: the target, or for B1111 the accuracy of the lower
            // ranks of the row; B0..B111 only beyond 1e±100.
            if double_error < DOUBLE_TOLERANCE {
                worst_double = worst_double.max(double_error);
            } else if (index < 12 && extreme(input))
                || (index >= 12 && double_error <= 4. * lower_double)
            {
                limited.push(message.clone());
            } else {
                failures.push(message.clone());
            }
            if float_error >= FLOAT_TOLERANCE {
                failures.push(message);
            }
        }
    }
    eprintln!(
        "{} rows against the exact b0..b1111: worst binary64 kernels {worst_f64:e} below the \
         target, generic binary64 {worst_generic:e}, DoubleFloat {worst_double:e}, Float(256) \
         {worst_float:e}; {beyond_320} rows need 1280 bits of the expressions, {} rows skipped \
         ({}); {} outputs above the target as for the generic code (binary64) or the lower \
         ranks (DoubleFloat):\n{}",
        points.len(),
        skipped.len(),
        skipped.join(", "),
        limited.len(),
        limited.join("\n")
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert!(nonfinite.is_empty(), "{}", nonfinite.join("\n"));
    assert!(skipped.len() <= 24, "{}", skipped.len());
}

/// The coarse sector of a fixture row, from its name.
fn fixture_sector(name: &str) -> String {
    let words: Vec<&str> = name.split('_').collect();
    match words[0] {
        "rank4" | "small" | "complex" | "normal" | "pseudo" => words[..2].join("_"),
        "m0" | "m1" | "p0" | "top" => words[..2].join("_"),
        _ => words[0].into(),
    }
}

/// Every route against the Feynman-parameter integrals, worst errors reported
/// per sector. The integrals carry 40 digits and agree with the 320-bit
/// algorithm to about 1e-27, which bounds the `DoubleFloat` check here (it is
/// compared with the exact expressions at 1e-28 above).
#[test]
fn every_route_matches_feynman_parameter_integrals() {
    let rows = feynman_rows();
    assert!(rows.len() >= 720, "{}", rows.len());
    let mut worst: BTreeMap<String, [f64; 4]> = BTreeMap::new();
    let mut failures = Vec::new();
    for row in &rows {
        let reference = evaluate_float(FAMILY, &row.input, 320);
        let sensitivity_53 = momentum_sensitivity(&row.input, &reference, 53);
        let sensitivity_106 = momentum_sensitivity(&row.input, &reference, 106);
        let fast = evaluate_f64(FAMILY, &row.input, true);
        let generic = evaluate_f64(FAMILY, &row.input, false);
        let double = evaluate_double(FAMILY, &row.input);
        let [p, m0, m1, _] = [0, 1, 2, 3].map(|i| row.input[i]);
        // As native_tensor.rs: absolute accuracy next to a massless on-shell point.
        let near_on_shell = [(m0, m1), (m1, m0)].iter().any(|(zero, mass)| {
            *zero == C::new(0., 0.)
                && mass.im == 0.
                && p.re != mass.re
                && (p.re - mass.re).abs() <= 1e-9 * mass.re.abs()
        });
        let entry = worst.entry(fixture_sector(&row.name)).or_default();
        let largest = row.expected.iter().map(float_norm).fold(0.0, f64::max);
        for (n, expected) in row.expected.iter().enumerate() {
            let size = if near_on_shell {
                1.0
            } else if n == 4 {
                largest + 1e-3
            } else {
                float_norm(expected) + 1e-3
            };
            let index = 3 * n;
            let float_error = difference(&reference[index], expected) / size;
            let scaled = |value: &Complex<Float>, sensitivity: f64| {
                (difference(value, expected) - 4. * sensitivity).max(0.) / size
            };
            let fast_error = scaled(&to_float(fast[index], 320), sensitivity_53[index]);
            let generic_error = scaled(&to_float(generic[index], 320), sensitivity_53[index]);
            let double_error = scaled(&double_to_float(&double[index]), sensitivity_106[index]);
            if n == 4 {
                for (slot, error) in [fast_error, generic_error, double_error, float_error]
                    .into_iter()
                    .enumerate()
                {
                    entry[slot] = entry[slot].max(error);
                }
            }
            if float_error >= 1e-24
                || fast_error >= BINARY64_TOLERANCE
                || generic_error >= BINARY64_TOLERANCE
                || double_error >= 1e-24
            {
                failures.push(format!(
                    "{} B rank {n}: integral {expected}, binary64 {:?} ({fast_error:.1e}), generic \
                     {:?} ({generic_error:.1e}), DoubleFloat {double_error:.1e}, 320 bits \
                     {float_error:.1e}",
                    row.name, fast[index], generic[index]
                ));
            }
        }
    }
    eprintln!(
        "B1111 against the integrals, worst scaled error per sector (binary64 kernels / generic binary64 / DoubleFloat / 320 bits):"
    );
    for (sector, errors) in &worst {
        eprintln!(
            "  {sector:<24} {}",
            errors.map(|e| format!("{e:.1e}")).join(" / ")
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// The native B1111 on the dyadic ladder p² = ±2^-k (k = 4..24, masses of
/// order one) against the integrals, next to the exact expression in
/// binary64, whose Passarino-Veltman form loses about `4 log10(m²/p²)` digits
/// (reported, not asserted: it is the reason fixed-precision callers use the
/// native family). Errors relative to `max_k |B_k|` of the row.
#[test]
fn small_momentum_ladder_keeps_binary64_accuracy() {
    let mut exact = symbolic_f64();
    let mut lines = Vec::new();
    let mut checked = 0;
    for row in feynman_rows()
        .iter()
        .filter(|row| row.name.starts_with("rank4_ladder_"))
    {
        let expected = &row.expected[4];
        let size = row.expected.iter().map(float_norm).fold(0.0, f64::max);
        let native = evaluate_f64(FAMILY, &row.input, true);
        let error = difference(&to_float(native[12], 320), expected) / size;
        assert!(
            error < 32. * f64::EPSILON,
            "{}: native {:?} vs integral {expected} ({error:e})",
            row.name,
            native[12]
        );
        assert_eq!(native[13], C::new(0.2, 0.), "{}", row.name);
        let symbolic = exact(&row.input);
        let symbolic_error = difference(&to_float(symbolic[12], 320), expected) / size;
        lines.push(format!(
            "  {:<32} native {error:.1e}, exact expression in binary64 {symbolic_error:.1e} \
             (its 1/eps coefficient {:.3e})",
            row.name, symbolic[13].re
        ));
        checked += 1;
    }
    eprintln!(
        "B1111 relative error on the p² ladder:\n{}",
        lines.join("\n")
    );
    assert!(checked >= 80, "{checked}");
}

/// `B_n(p; m0, m1) = (-1)^n sum_k C(n,k) B_k(p; m1, m0)` at every rank, on
/// every route, in every branch of the algorithm and in both mass orders (one
/// of which takes the reconstruction).
#[test]
fn mass_swap_identity_holds_in_every_branch() {
    const BINOMIAL: [[f64; 5]; 5] = [
        [1., 0., 0., 0., 0.],
        [1., 1., 0., 0., 0.],
        [1., 2., 1., 0., 0.],
        [1., 3., 3., 1., 0.],
        [1., 4., 6., 4., 1.],
    ];
    let points = all_points(4000);
    let mut coverage: BTreeMap<(&str, bool), usize> = BTreeMap::new();
    let mut worst = [0.0_f64; 4];
    let mut failures = Vec::new();
    // The identity on values `direct` (order m0, m1) and `swapped` (m1, m0),
    // normalized by the size of its terms.
    let check = |direct: &[Complex<Float>], swapped: &[Complex<Float>]| -> f64 {
        let mut error = 0.0_f64;
        for n in 0..5 {
            let mut sum = to_float(C::new(0., 0.), 320);
            // Each term with the 1e-3 floor of the other tests: values far
            // below one (equal masses at small p²) carry absolute errors.
            let mut terms = float_norm(&direct[3 * n]) + 1e-3;
            for k in 0..=n {
                let weight =
                    Float::with_val(320, BINOMIAL[n][k] * if n % 2 == 0 { 1. } else { -1. });
                sum = Complex::new(
                    sum.re + swapped[3 * k].re.clone() * &weight,
                    sum.im + swapped[3 * k].im.clone() * &weight,
                );
                terms += BINOMIAL[n][k] * (float_norm(&swapped[3 * k]) + 1e-3);
            }
            error = error.max(difference(&direct[3 * n], &sum) / terms);
        }
        error
    };
    for (index, (name, input)) in points.iter().enumerate() {
        let swapped_input = vec![input[0], input[2], input[1], input[3]];
        *coverage.entry(sector(input)).or_default() += 1;
        let as_float = |values: Vec<C>| -> Vec<Complex<Float>> {
            values.into_iter().map(|z| to_float(z, 320)).collect()
        };
        let as_float_double = |values: Vec<D>| -> Vec<Complex<Float>> {
            values.iter().map(double_to_float).collect()
        };
        let mut errors = [
            check(
                &as_float(evaluate_f64(FAMILY, input, true)),
                &as_float(evaluate_f64(FAMILY, &swapped_input, true)),
            ),
            check(
                &as_float(evaluate_f64(FAMILY, input, false)),
                &as_float(evaluate_f64(FAMILY, &swapped_input, false)),
            ),
            check(
                &as_float_double(evaluate_double(FAMILY, input)),
                &as_float_double(evaluate_double(FAMILY, &swapped_input)),
            ),
            0.,
        ];
        if index % 4 == 0 {
            errors[3] = check(
                &evaluate_float(FAMILY, input, 320),
                &evaluate_float(FAMILY, &swapped_input, 320),
            );
        }
        // Tens of ulps of the terms; DoubleFloat's logarithm keeps about 90
        // bits beyond 1e±100, and Float tracks the cancellation next to a
        // massless on-shell point (about 270 of 320 bits remain).
        let double_tolerance = if extreme(input) { 1e-24 } else { 1e-29 };
        for (slot, (error, tolerance)) in errors
            .iter()
            .zip([1e-14, 1e-14, double_tolerance, 1e-80])
            .enumerate()
        {
            worst[slot] = worst[slot].max(*error);
            if *error > tolerance {
                failures.push(format!(
                    "{name} {input:?} ({:?}): route {slot} violates the identity by {error:e}",
                    sector(input)
                ));
            }
        }
    }
    eprintln!(
        "mass-swap identity through rank four: worst {} (binary64 kernels / generic binary64 / \
         DoubleFloat / 320 bits); rows per (branch, reconstructed): {coverage:?}",
        worst.map(|e| format!("{e:.1e}")).join(" / ")
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    for branch in [
        "p0_unequal",
        "p0_massless",
        "massless",
        "massless_on_shell",
        "real_roots",
        "conjugate_roots",
        "complex_masses",
    ] {
        for reconstructed in [false, true] {
            let count = coverage.get(&(branch, reconstructed)).copied().unwrap_or(0);
            assert!(
                count >= 5,
                "{branch} reconstructed={reconstructed}: {count} rows"
            );
        }
    }
}

/// Reference values of the exact-expression tests (`src/tests.rs`) and of the
/// masters study, on every route. `B1111(0; 2, 5)` takes the zero-momentum
/// reconstruction, where OneLOop returns `+0.2998...`.
#[test]
fn reference_values() {
    let cases = [
        (
            [0., 0., 2., 0., 5., 0.],
            C::new(-0.299_847_502_997_085_9, 0.),
        ),
        (
            [0., 0., 5., 0., 2., 0.],
            C::new(-0.180_652_923_442_634_8, 0.),
        ),
        (
            [3., 0., 0.7, -0.03, 1.4, -0.08],
            C::new(0.029_939_862_301_517_71, 0.017_531_313_325_083_214),
        ),
        (
            [3., 0., 1.5, 0., 0.75, 0.],
            C::new(0.138_383_602_547_043_8, 0.),
        ),
        (
            [2.5, 0., 1.5, 0., 0., 0.],
            C::new(0.591_880_873_305_368_5, 0.579_460_481_769_330_2),
        ),
    ];
    for (values, expected) in cases {
        let input = [
            C::new(values[0], values[1]),
            C::new(values[2], values[3]),
            C::new(values[4], values[5]),
            C::new(1., 0.),
        ];
        let tolerance = 64. * f64::EPSILON * norm(expected);
        for binary64 in [true, false] {
            let output = evaluate_f64(FAMILY, &input, binary64);
            assert!(
                norm(output[12] - expected) < tolerance.max(1e-15),
                "{input:?} binary64 kernels {binary64}: {:?} vs {expected:?}",
                output[12]
            );
            assert_eq!(output[13], C::new(0.2, 0.));
        }
        let double = double_to_float(&evaluate_double(FAMILY, &input)[12]);
        let float = evaluate_float(FAMILY, &input, 256);
        assert!(difference(&double, &float[12]) < 1e-30, "{input:?}");
        assert!(
            difference(&float[12], &to_float(expected, 256)) < 1e-15,
            "{input:?}"
        );
    }
}

#[test]
fn family_shape() {
    assert_eq!(FAMILY.name(), "B1111");
    assert_eq!(FAMILY.arity(), 4);
    assert_eq!(FAMILY.coefficients(), 5);
    assert_eq!(FAMILY.outputs(), 15);
    assert_eq!(
        FAMILY.coefficient_names(),
        &["B0", "B1", "B11", "B111", "B1111"]
    );
    let mut evaluator = NativeTensorEvaluator::<f64>::new(FAMILY).unwrap();
    assert_eq!(evaluator.family(), FAMILY);
    let input = [
        C::new(3., 0.),
        C::new(1., 0.),
        C::new(2., 0.),
        C::new(1., 0.),
    ];
    let mut short = [C::new(0., 0.); 12];
    assert!(evaluator.evaluate(&input, &mut short).is_err());
    let mut positive_width = input;
    positive_width[2].im = 0.5;
    let mut output = [C::new(7., 7.); 15];
    assert!(evaluator.evaluate(&positive_width, &mut output).is_err());
    assert!(output.iter().all(|z| *z == C::new(7., 7.)));
    assert!(evaluator.evaluate(&input, &mut output).is_ok());
}
