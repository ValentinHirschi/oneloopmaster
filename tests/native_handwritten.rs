//! The hand-written native C0/D0 sector kernels against the generated code of
//! the exact expressions, in binary64 and at 256 bits, on randomized kinematics
//! covering the dispatched sectors, exact-zero sub-sectors and thresholds.
#![cfg(feature = "generated-evaluators")]

use oneloop::{NativeEvaluator, ScalarIntegral};
use symbolica::domains::float::{Complex, Float};
use symbolica::domains::float::{FloatLike, RealLike};

type C = Complex<f64>;

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
    /// A "nice" value: small integers make exact zeros, equalities and
    /// thresholds frequent, as in physical S-matrices.
    fn nice(&mut self, scale: f64) -> f64 {
        match self.below(4) {
            0 => (self.below(9) as f64 - 4.0) * scale / 4.0,
            _ => (2.0 * self.uniform() - 1.0) * scale,
        }
    }
}

fn real(x: f64) -> C {
    C::new(x, 0.0)
}

/// Random triangle inputs [p1,p2,p3,m1,m2,m3,mu2] with three or two nonzero masses.
fn triangle_point(rng: &mut Rng) -> Vec<C> {
    let mass = |rng: &mut Rng, common: f64| -> C {
        let re = match rng.below(3) {
            0 => common,
            1 => (1 + rng.below(4)) as f64 * 0.25,
            _ => 0.05 + rng.uniform(),
        };
        let im = if rng.below(4) == 0 { -0.1 * rng.uniform() } else { 0.0 };
        C::new(re, im)
    };
    let common = 0.05 + rng.uniform();
    let mut m: Vec<C> = (0..3).map(|_| mass(rng, common)).collect();
    if rng.below(3) == 0 {
        // Two-mass sectors: one massless propagator.
        m[rng.below(3)] = C::new(0.0, 0.0);
    }
    let momentum = |rng: &mut Rng, m: &[C]| -> f64 {
        match rng.below(6) {
            0 => 0.0,
            1 => m[rng.below(3)].re,
            2 => {
                let (a, b) = (m[rng.below(3)].re.sqrt(), m[rng.below(3)].re.sqrt());
                (a + b) * (a + b)
            }
            3 => -rng.uniform() * 4.0,
            _ => rng.nice(6.0),
        }
    };
    let mut p: Vec<C> = (0..3).map(|_| real(momentum(rng, &m))).collect();
    if rng.below(5) == 0 {
        p[rng.below(3)] = p[rng.below(3)];
    }
    p.extend(m);
    p.push(real(0.3 + rng.uniform()));
    p
}

/// Random box inputs [p1..p4,s12,s23,m1..m4,mu2] with four or three nonzero masses.
fn box_point(rng: &mut Rng, complex: bool) -> Vec<C> {
    let common = 0.05 + rng.uniform();
    let m: Vec<C> = (0..4)
        .map(|_| {
            let re = match rng.below(3) {
                0 | 1 => common,
                _ => 0.05 + rng.uniform(),
            };
            let im = if complex && rng.below(2) == 0 {
                -0.05 * rng.uniform()
            } else {
                0.0
            };
            C::new(re, im)
        })
        .collect();
    let mut m = m;
    if rng.below(3) == 0 {
        // Three-mass sectors: one massless propagator.
        m[rng.below(4)] = C::new(0.0, 0.0);
    }
    let mut p: Vec<C> = (0..6)
        .map(|_| {
            real(match rng.below(6) {
                0 => 0.0,
                1 => m[rng.below(4)].re,
                2 => -rng.uniform() * 3.0,
                _ => rng.nice(8.0),
            })
        })
        .collect();
    if rng.below(4) == 0 {
        p[rng.below(6)] = p[rng.below(6)];
    }
    p.extend(m);
    p.push(real(0.3 + rng.uniform()));
    p
}

fn evaluate_f64(evaluator: &mut NativeEvaluator<f64>, hand: bool, input: &[C]) -> [C; 3] {
    evaluator.set_handwritten(hand);
    let mut out = [C::new(0.0, 0.0); 3];
    evaluator.evaluate(input, &mut out).unwrap();
    out
}

fn evaluate_float(
    evaluator: &mut NativeEvaluator<Float>,
    hand: bool,
    input: &[C],
    bits: u32,
) -> Vec<Complex<Float>> {
    evaluator.set_handwritten(hand);
    let input: Vec<Complex<Float>> = input
        .iter()
        .map(|z| Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im)))
        .collect();
    let mut out: Vec<Complex<Float>> = (0..3)
        .map(|_| Complex::new(Float::new(bits), Float::new(bits)))
        .collect();
    evaluator.evaluate(&input, &mut out).unwrap();
    out
}

/// The smallest tracked precision of the nonzero output components.
fn tracked_bits(values: &[Complex<Float>]) -> u32 {
    values
        .iter()
        .flat_map(|z| [&z.re, &z.im])
        .filter(|x| !x.is_fully_zero())
        .map(|x| x.prec())
        .min()
        .unwrap_or(u32::MAX)
}

fn to_f64(values: &[Complex<Float>]) -> Vec<C> {
    values.iter().map(|z| C::new(z.re.to_f64(), z.im.to_f64())).collect()
}

/// max_k |a_k - b_k| / max_k |b_k|, in the precision of the inputs.
fn float_difference(a: &[Complex<Float>], b: &[Complex<Float>]) -> f64 {
    let scale = to_f64(b).iter().map(|z| z.re.hypot(z.im)).fold(0.0, f64::max);
    a.iter()
        .zip(b)
        .map(|(x, y)| {
            let re = (x.re.clone() - &y.re).to_f64();
            let im = (x.im.clone() - &y.im).to_f64();
            re.hypot(im) / scale.max(1e-300)
        })
        .fold(0.0, |w: f64, d| if d.is_nan() { f64::INFINITY } else { w.max(d) })
}

fn scaled_difference(a: &[C], b: &[C]) -> f64 {
    let scale = b.iter().map(|z| z.re.hypot(z.im)).fold(0.0, f64::max).max(1e-300);
    a.iter()
        .zip(b)
        .map(|(x, y)| (x.re - y.re).hypot(x.im - y.im) / scale)
        .fold(0.0, |w: f64, d| if d.is_nan() { f64::INFINITY } else { w.max(d) })
}

/// Binary64 accuracy of the hand-written kernels against the generated
/// expressions at 256 (or 512) bits, compared with the generated binary64 code;
/// and the kernels' logic at 256 bits against the same reference wherever both
/// retain at least 120 tracked bits.
fn run(family: ScalarIntegral, points: &[Vec<C>]) {
    let mut fast = NativeEvaluator::<f64>::new(family).unwrap();
    let mut precise = NativeEvaluator::<Float>::with_binary_precision(family, 256).unwrap();
    let mut very_precise = NativeEvaluator::<Float>::with_binary_precision(family, 512).unwrap();
    let mut failures = vec![];
    let (mut worst_hand, mut worst_generated, mut worst_logic) = (0.0f64, 0.0f64, 0.0f64);
    let (mut unresolved, mut logic_checked, mut nonfinite) = (0, 0, 0);
    let (mut worse, mut better) = (0, 0);
    for input in points {
        let mut truth = evaluate_float(&mut precise, false, input, 256);
        if tracked_bits(&truth) < 80 {
            truth = evaluate_float(&mut very_precise, false, input, 512);
        }
        let truth_bits = tracked_bits(&truth);
        let reference = to_f64(&truth);
        if truth_bits < 60 || reference.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
            unresolved += 1;
            continue;
        }
        let hand = evaluate_f64(&mut fast, true, input);
        let generated = evaluate_f64(&mut fast, false, input);
        let error_hand = scaled_difference(&hand, &reference);
        let error_generated = scaled_difference(&generated, &reference);
        if !error_generated.is_finite() {
            nonfinite += 1;
        }
        worst_hand = worst_hand.max(error_hand);
        if error_generated.is_finite() {
            worst_generated = worst_generated.max(error_generated);
        }
        // Near thresholds both codes lose digits to the same conditioning,
        // and which one is luckier varies point by point: count those cases
        // both ways. A gross deviation (a wrong sheet or branch) fails.
        if error_hand > 1e-12 && error_hand > 16.0 * error_generated {
            worse += 1;
        }
        if error_generated > 1e-12 && error_generated > 16.0 * error_hand {
            better += 1;
        }
        if error_hand > 1e-7 && !(error_hand <= 16.0 * error_generated) {
            failures.push(format!(
                "binary64 error {error_hand:.2e} (generated {error_generated:.2e}): {input:?}\n  hand={hand:?}\n  generated={generated:?}\n  reference={reference:?}"
            ));
        }
        let logic = evaluate_float(&mut precise, true, input, 256);
        let bits = tracked_bits(&logic).min(truth_bits);
        if bits >= 120 {
            logic_checked += 1;
            let d = float_difference(&logic, &truth);
            worst_logic = worst_logic.max(d);
            if !(d <= 2f64.powi(12 - bits as i32)) {
                failures.push(format!(
                    "256-bit difference {d:.2e} ({bits} bits): {input:?}\n  hand={:?}\n  reference={reference:?}",
                    to_f64(&logic)
                ));
            }
        }
    }
    eprintln!(
        "{family:?}: {} points ({unresolved} unresolved, {nonfinite} nonfinite generated binary64): worst binary64 error hand {worst_hand:.2e}, generated {worst_generated:.2e}; hand-written more than 16x less accurate on {worse}, more than 16x more accurate on {better}; {logic_checked} 256-bit logic checks, worst {worst_logic:.2e}",
        points.len()
    );
    assert!(
        worse <= better + points.len() / 500,
        "hand-written kernels systematically less accurate: {worse} against {better}"
    );
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.iter().take(count(20, "HANDWRITTEN_SHOW")).cloned().collect::<Vec<_>>().join("\n")
    );
}

fn count(default: usize, variable: &str) -> usize {
    std::env::var(variable)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

#[test]
fn handwritten_c0_matches_generated_expressions() {
    let mut rng = Rng(0x9e3779b97f4a7c15);
    let n = count(3000, "HANDWRITTEN_POINTS");
    let points: Vec<_> = (0..n).map(|_| triangle_point(&mut rng)).collect();
    run(ScalarIntegral::C0, &points);
}

#[test]
fn handwritten_d0_matches_generated_expressions() {
    let mut rng = Rng(0x2545f4914f6cdd1d);
    let n = count(600, "HANDWRITTEN_BOX_POINTS");
    let points: Vec<_> = (0..n).map(|i| box_point(&mut rng, i % 3 == 0)).collect();
    run(ScalarIntegral::D0, &points);
}

/// Every distinct C0/D0 call of the recorded Ninja benchmark streams
/// (tests/data/ninja_master_rows.txt, with the OneLOop 3.6 results). In the
/// hand-written sectors the binary64 kernels are compared with OneLOop and,
/// where the generated expressions keep at least 60 tracked bits at 256 or
/// 1024 bits, with that reference; each against the generated binary64 code.
#[test]
fn handwritten_kernels_on_ninja_streams() {
    let mut failures = vec![];
    for (family, label) in [(ScalarIntegral::C0, "3"), (ScalarIntegral::D0, "4")] {
        let rows: Vec<(Vec<C>, Vec<C>)> = include_str!("data/ninja_master_rows.txt")
            .lines()
            .filter(|line| line.split_whitespace().next() == Some(label))
            .map(|line| {
                let values: Vec<f64> = line
                    .split('#')
                    .next()
                    .unwrap()
                    .split_whitespace()
                    .skip(1)
                    .map(|v| v.parse().unwrap())
                    .collect();
                let momenta = if label == "3" { 3 } else { 6 };
                let masses = if label == "3" { 3 } else { 4 };
                let mut input: Vec<C> = values[1..=momenta].iter().map(|&x| real(x)).collect();
                let offset = 1 + momenta;
                input.extend((0..masses).map(|i| C::new(values[offset + 2 * i], values[offset + 2 * i + 1])));
                input.push(real(values[0]));
                let offset = offset + 2 * masses;
                let oneloop = (0..3).map(|i| C::new(values[offset + 2 * i], values[offset + 2 * i + 1])).collect();
                (input, oneloop)
            })
            .collect();
        let mut fast = NativeEvaluator::<f64>::new(family).unwrap();
        let mut precise = NativeEvaluator::<Float>::with_binary_precision(family, 256).unwrap();
        let mut very_precise = NativeEvaluator::<Float>::with_binary_precision(family, 1024).unwrap();
        let median = |mut v: Vec<f64>| {
            v.sort_by(|a, b| a.total_cmp(b));
            v.get(v.len() / 2).copied().unwrap_or(0.0)
        };
        let worst = |v: &[f64]| v.iter().cloned().fold(0.0, f64::max);
        let (mut covered, mut resolved) = (0, 0);
        let (mut hand_oneloop, mut generated_oneloop, mut hand_truth, mut generated_truth) =
            (vec![], vec![], vec![], vec![]);
        for (input, oneloop) in &rows {
            let hand = evaluate_f64(&mut fast, true, input);
            let generated = evaluate_f64(&mut fast, false, input);
            let bits = |v: &[C; 3]| v.map(|z| (z.re.to_bits(), z.im.to_bits()));
            if bits(&hand) == bits(&generated) {
                continue; // a sector of the generated code
            }
            covered += 1;
            let (h, g) = (scaled_difference(&hand, oneloop), scaled_difference(&generated, oneloop));
            hand_oneloop.push(h);
            generated_oneloop.push(g);
            if !(h <= 1e-12 || h <= 16.0 * g) {
                failures.push(format!("{family:?} vs OneLOop: {h:.2e} (generated {g:.2e}): {input:?}"));
            }
            // A reference keeps at least 60 tracked bits (rounding of an exact
            // zero coefficient can cost all of them, e.g. c in trif3).
            let mut truth = evaluate_float(&mut precise, false, input, 256);
            if tracked_bits(&truth) < 80 {
                truth = evaluate_float(&mut very_precise, false, input, 1024);
            }
            if tracked_bits(&truth) < 60 {
                continue;
            }
            resolved += 1;
            let truth = to_f64(&truth);
            let (h, g) = (scaled_difference(&hand, &truth), scaled_difference(&generated, &truth));
            hand_truth.push(h);
            generated_truth.push(g);
            if !(h <= 1e-12 || h <= 16.0 * g) {
                failures.push(format!("{family:?} vs high precision: {h:.2e} (generated {g:.2e}): {input:?}"));
            }
        }
        eprintln!(
            "{family:?}: {} distinct rows, {covered} in hand-written sectors. Against OneLOop 3.6, median / worst: hand-written {:.1e} / {:.2e}, generated {:.1e} / {:.2e}. Against the expressions at high precision ({resolved} resolved rows): hand-written {:.1e} / {:.2e}, generated {:.1e} / {:.2e}",
            rows.len(),
            median(hand_oneloop.clone()),
            worst(&hand_oneloop),
            median(generated_oneloop.clone()),
            worst(&generated_oneloop),
            median(hand_truth.clone()),
            worst(&hand_truth),
            median(generated_truth.clone()),
            worst(&generated_truth),
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
