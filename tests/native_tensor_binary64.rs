//! The binary64 tensor kernels (`NativeTensorEvaluator<f64>`) against the
//! generic code of the same algorithm in binary64 (`set_binary64(false)`), at
//! 320 bits and against the exact symbolic expressions (`b0`, `b1`, `b11`,
//! `b111`, `an`) at 256 bits: randomized kinematics over every sector, and
//! every distinct bubble/tadpole call of Ninja's recorded master streams,
//! also against OneLOop.
use oneloop::{LaurentSeries, NativeTensorEvaluator, TensorIntegral};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::float::{Complex, Float, RealLike},
    parse,
};

type C = Complex<f64>;

/// Binary64 accuracy target against 320 bits, as in `native_tensor.rs`:
/// relative to the coefficient plus 1e-3 of its natural scale.
const BINARY64_TOLERANCE: f64 = 5e-14;

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

fn norm(z: C) -> f64 {
    z.re.hypot(z.im)
}

fn to_float(z: C, bits: u32) -> Complex<Float> {
    Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im))
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

fn evaluate_float(family: TensorIntegral, input: &[C], bits: u32) -> Vec<Complex<Float>> {
    let mut evaluator =
        NativeTensorEvaluator::<Float>::with_binary_precision(family, bits).unwrap();
    let input: Vec<_> = input.iter().map(|&z| to_float(z, bits)).collect();
    let mut output = vec![to_float(C::new(0., 0.), bits); family.outputs()];
    evaluator.evaluate(&input, &mut output).unwrap();
    output
}

/// Largest change of each output at 320 bits when p² moves by a relative
/// 2^-53: threshold rows are ill-conditioned (`Im B ~ sqrt(p² - threshold)`),
/// so no binary64 algorithm can do better for an input rounded to binary64.
fn momentum_sensitivity(input: &[C], reference: &[Complex<Float>]) -> Vec<f64> {
    if input[0].re == 0. {
        return vec![0.; reference.len()];
    }
    let mut sensitivity = vec![0.0_f64; reference.len()];
    for sign in [-1.0, 1.0] {
        let mut evaluator =
            NativeTensorEvaluator::<Float>::with_binary_precision(TensorIntegral::B111, 320)
                .unwrap();
        let mut shifted: Vec<_> = input.iter().map(|&z| to_float(z, 320)).collect();
        let factor = Float::with_val(320, 1) + Float::with_val(320, sign * 2f64.powi(-53));
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

/// `[p², m0², m1², mu²]` covering every sector: zero momentum, a massless
/// propagator (on and off shell), equal masses, normal and pseudo
/// thresholds, small and large momenta (series and direct `bnlog`), real
/// negative and complex masses.
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

/// Errors of the binary64 kernels and of the generic binary64 code, scaled
/// like `native_tensor.rs`, after removing the input conditioning.
fn errors(input: &[C]) -> (Vec<f64>, Vec<f64>, Vec<C>, Vec<C>) {
    let family = TensorIntegral::B111;
    let fast = evaluate_f64(family, input, true);
    let generic = evaluate_f64(family, input, false);
    let reference = evaluate_float(family, input, 320);
    let sensitivity = momentum_sensitivity(input, &reference);
    let scaled = |values: &[C]| -> Vec<f64> {
        values
            .iter()
            .zip(&reference)
            .zip(&sensitivity)
            .map(|((value, reference), sensitivity)| {
                let size = float_norm(reference) + 1e-3;
                ((difference(&to_float(*value, 320), reference) - 4. * sensitivity).max(0.)) / size
            })
            .collect()
    };
    (scaled(&fast), scaled(&generic), fast, generic)
}

/// Rounding noise of a scaled error: both codes round differently, and a
/// cancelling reconstruction (B111 after a mass swap) amplifies an ulp.
const NOISE: f64 = 1e-14;

#[test]
fn binary64_kernels_match_generic_code_and_higher_precision() {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let (mut all_fast, mut all_generic) = (Vec::new(), Vec::new());
    let (mut fast_better, mut generic_better) = (0, 0);
    let mut failures = Vec::new();
    for _ in 0..8000 {
        let input = bubble_point(&mut rng);
        let (fast, generic, fast_values, _) = errors(&input);
        // B11 is the leading part of B111, bitwise.
        assert_eq!(
            evaluate_f64(TensorIntegral::B11, &input, true),
            fast_values[..9],
            "{input:?}"
        );
        for (index, (&fast, &generic)) in fast.iter().zip(&generic).enumerate() {
            if fast > 16. * generic.max(f64::EPSILON) && fast > NOISE {
                generic_better += 1;
            }
            if generic > 16. * fast.max(f64::EPSILON) && generic > NOISE {
                fast_better += 1;
            }
            // The generic code's own accuracy, or the binary64 target.
            if fast > (4. * generic).max(BINARY64_TOLERANCE) {
                failures.push(format!(
                    "{input:?} output {index}: binary64 kernels {fast:e}, generic {generic:e}"
                ));
            }
            all_fast.push(fast);
            all_generic.push(generic);
        }
    }
    let quantiles = |values: &mut Vec<f64>| {
        values.sort_by(f64::total_cmp);
        [0.5, 0.9, 0.99, 0.999, 1.0].map(|q| values[((values.len() - 1) as f64 * q) as usize])
    };
    let (fast, generic) = (quantiles(&mut all_fast), quantiles(&mut all_generic));
    let show = |values: &[f64; 5]| values.map(|v| format!("{v:.1e}")).join("/");
    eprintln!(
        "{} outputs, scaled error quantiles 50/90/99/99.9/100%: binary64 kernels {}, generic {}; \
         above {NOISE:e}, kernels 16x more accurate on {fast_better} outputs, 16x less on {generic_better}",
        all_fast.len(),
        show(&fast),
        show(&generic)
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert_eq!(generic_better, 0);
    // No systematic loss: every quantile within twice the generic code's.
    for (fast, generic) in fast.iter().zip(&generic) {
        assert!(
            *fast <= 2. * generic.max(f64::EPSILON),
            "{fast:e} vs {generic:e}"
        );
    }
}

#[test]
fn binary64_tadpoles_match_generic_code() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    for _ in 0..2000 {
        let scale = rng.log_uniform(1e-6, 1e6);
        let input = [
            mass(&mut rng, scale),
            C::new(rng.log_uniform(1e-2, 1e2) * scale, 0.),
        ];
        let fast = evaluate_f64(TensorIntegral::An, &input, true);
        let generic = evaluate_f64(TensorIntegral::An, &input, false);
        assert_eq!(fast, generic, "{input:?}");
    }
}

struct StreamRow {
    case: String,
    family: TensorIntegral,
    input: Vec<C>,
    /// OneLOop: A0, or B0, B1, B11 (olo_b11) as Laurent triples.
    oneloop: Vec<C>,
}

/// Every distinct bubble/tadpole call of Ninja's recorded master streams
/// (tests/generate_ninja_tensor_rows.py).
fn stream_rows() -> Vec<StreamRow> {
    let mut rows = Vec::new();
    for line in include_str!("data/ninja_tensor_rows.txt").lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let (values, case) = line.split_once('#').unwrap();
        let values: Vec<f64> = values
            .split_whitespace()
            .map(|v| v.parse().unwrap())
            .collect();
        let mu = C::new(values[1], 0.);
        let pairs = |v: &[f64]| {
            v.chunks_exact(2)
                .map(|p| C::new(p[0], p[1]))
                .collect::<Vec<_>>()
        };
        let (family, input, oneloop) = match values[0] {
            5. => (
                TensorIntegral::An,
                vec![C::new(values[2], values[3]), mu],
                pairs(&values[4..]),
            ),
            6. | 7. => (
                if values[0] == 6. {
                    TensorIntegral::B11
                } else {
                    TensorIntegral::B111
                },
                vec![
                    C::new(values[2], 0.),
                    C::new(values[3], values[4]),
                    C::new(values[5], values[6]),
                    mu,
                ],
                pairs(&values[7..]),
            ),
            _ => panic!("invalid stream row kind"),
        };
        rows.push(StreamRow {
            case: case.trim().into(),
            family,
            input,
            oneloop,
        });
    }
    rows
}

#[test]
fn ninja_stream_rows_match_oneloop_and_higher_precision() {
    let rows = stream_rows();
    assert_eq!(rows.len(), 290);
    let (mut worst_oneloop, mut worst_reference, mut worst_generic) = (0.0_f64, 0.0_f64, 0.0_f64);
    let mut attributed = Vec::new();
    for row in &rows {
        let actual = evaluate_f64(row.family, &row.input, true);
        let reference = evaluate_float(row.family, &row.input, 320);
        let generic = evaluate_f64(row.family, &row.input, false);
        let magnitude = |index: usize| {
            if row.family == TensorIntegral::An {
                norm(row.input[0]).powi(index as i32 / 3 + 1)
            } else {
                1.0
            }
        };
        for index in 0..row.family.outputs() {
            let size = float_norm(&reference[index]) + 1e-3 * magnitude(index);
            let error = difference(&to_float(actual[index], 320), &reference[index]) / size;
            let generic_error =
                difference(&to_float(generic[index], 320), &reference[index]) / size;
            worst_reference = worst_reference.max(error);
            worst_generic = worst_generic.max(generic_error);
            assert!(
                error < BINARY64_TOLERANCE && error <= (4. * generic_error).max(NOISE),
                "{} {:?} output {index}: {error:e} against 320 bits (generic {generic_error:e})",
                row.case,
                row.input
            );
            // OneLOop's A0 or B0, B1, B11; its direct bnlog formula cancels
            // near |y| = 1/(8(n+1)), so attribute a difference to OneLOop
            // where this kernel matches 320 bits.
            if let Some(expected) = row.oneloop.get(index) {
                let error =
                    norm(actual[index] - *expected) / (norm(*expected) + 1e-3 * magnitude(index));
                worst_oneloop = worst_oneloop.max(error);
                if error > 1e-13 {
                    attributed.push(format!("{} output {index} ({error:.1e})", row.case));
                    assert!(
                        error < 1e-10,
                        "{} {:?} output {index}: {error:e} against OneLOop",
                        row.case,
                        row.input
                    );
                }
            }
        }
    }
    eprintln!(
        "Ninja stream rows: worst against 320 bits {worst_reference:e} (generic code {worst_generic:e}), \
         against OneLOop {worst_oneloop:e}; \
         OneLOop differs by more than 1e-13 where this kernel matches 320 bits: {}",
        attributed.join(", ")
    );
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

type Symbolic = Box<dyn FnMut(&[C]) -> Vec<Complex<Float>>>;

/// The exact expressions of `family`, evaluated at `bits`.
fn symbolic(family: TensorIntegral, bits: u32) -> Symbolic {
    let parameters: Vec<Atom> = (0..family.arity())
        .map(|i| parse!(format!("binary64_tensor_parameter_{}_{i}", family.name())))
        .collect();
    let coefficients = symbolic_coefficients(family, &parameters);
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
        let mut output = vec![to_float(C::new(0., 0.), bits); family.outputs()];
        evaluator.evaluate(&input, &mut output);
        output
    })
}

/// Every regression row (Ninja's streams, random bubbles and tadpoles) in
/// binary64 against the symbolic expressions at 256 bits. Where the
/// expressions themselves cancel beyond 256 bits (Passarino-Veltman forms at
/// tiny p², like `(m²/p²)^rank`) the 1024-bit evaluation is the reference.
///
/// This also checks two repaired expression defects of `d10eaa0`: the
/// zero-momentum B1/B11/B111 anchors with real masses of opposite sign took
/// the wrong lip, and the two-mass B0 formula at `m1² = 0` was not finite
/// where its root `x = 1` rounded to `1 - ulp`, so that `log(1 - 1/x) =
/// log(0)` (e.g. `B(5937.799487939787; 29929, 0)` at 256 bits).
#[test]
fn binary64_kernels_match_symbolic_expressions() {
    let mut rows: Vec<(TensorIntegral, Vec<C>)> = stream_rows()
        .into_iter()
        .map(|row| (row.family, row.input))
        .collect();
    let mut rng = Rng(0x6a09_e667_f3bc_c908);
    for _ in 0..1500 {
        rows.push((TensorIntegral::B111, bubble_point(&mut rng)));
    }
    for _ in 0..300 {
        let scale = rng.log_uniform(1e-6, 1e6);
        let mass = mass(&mut rng, scale);
        rows.push((
            TensorIntegral::An,
            vec![mass, C::new(rng.log_uniform(1e-2, 1e2) * scale, 0.)],
        ));
    }
    let families = [
        TensorIntegral::An,
        TensorIntegral::B11,
        TensorIntegral::B111,
    ];
    let mut at_256: Vec<Symbolic> = families.iter().map(|&f| symbolic(f, 256)).collect();
    let mut at_1024: Vec<Symbolic> = families.iter().map(|&f| symbolic(f, 1024)).collect();
    let (mut worst, mut worst_generic, mut beyond_256) = (0.0_f64, 0.0_f64, 0);
    let (mut nonfinite, mut limited, mut failures) = (Vec::new(), Vec::new(), Vec::new());
    for (family, input) in &rows {
        let slot = families.iter().position(|f| f == family).unwrap();
        let actual = evaluate_f64(*family, input, true);
        let generic = evaluate_f64(*family, input, false);
        let coarse = at_256[slot](input);
        let fine = at_1024[slot](input);
        let sensitivity = if *family == TensorIntegral::An {
            vec![0.; 12]
        } else {
            let reference = evaluate_float(TensorIntegral::B111, input, 320);
            momentum_sensitivity(input, &reference)
        };
        for index in 0..family.outputs() {
            let magnitude = if *family == TensorIntegral::An && norm(input[0]) > 0. {
                norm(input[0]).powi(index as i32 / 3 + 1)
            } else {
                1.0
            };
            let finite =
                |z: &Complex<Float>| z.re.to_f64().is_finite() && z.im.to_f64().is_finite();
            if !finite(&coarse[index]) || !finite(&fine[index]) {
                nonfinite.push(format!(
                    "{} {input:?} output {index}: binary64 {:?}, expression {} (256 bits) / {} (1024 bits)",
                    family.name(),
                    actual[index],
                    coarse[index],
                    fine[index]
                ));
                continue;
            }
            let size = float_norm(&fine[index]) + 1e-3 * magnitude;
            let expression = if difference(&coarse[index], &fine[index]) <= 1e-40 * size {
                &coarse[index]
            } else {
                beyond_256 += 1;
                &fine[index]
            };
            let scaled = |value: C| {
                (difference(&to_float(value, 320), expression) - 4. * sensitivity[index]).max(0.)
                    / size
            };
            let (error, generic_error) = (scaled(actual[index]), scaled(generic[index]));
            worst_generic = worst_generic.max(generic_error);
            let message = format!(
                "{} {input:?} output {index}: binary64 {:?} vs expression {} ({error:.1e}, \
                 generic binary64 code {generic_error:.1e})",
                family.name(),
                actual[index],
                expression
            );
            // The binary64 target, or the generic code's own accuracy where the
            // algorithm loses digits (e.g. the root pair at p² << m²).
            if error < BINARY64_TOLERANCE {
                worst = worst.max(error);
            } else if error <= 4. * generic_error {
                limited.push(message);
            } else {
                failures.push(message);
            }
        }
    }
    eprintln!(
        "{} rows against the symbolic expressions: worst binary64 error {worst:e} \
         below the target; {} outputs above it as for the generic binary64 code \
         (worst generic error {worst_generic:e}):\n{}\n\
         {beyond_256} outputs need more than 256 bits of the expressions",
        rows.len(),
        limited.len(),
        limited.join("\n"),
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert!(
        nonfinite.is_empty(),
        "{} outputs where the expression is not finite:\n{}",
        nonfinite.len(),
        nonfinite.join("\n")
    );
}
