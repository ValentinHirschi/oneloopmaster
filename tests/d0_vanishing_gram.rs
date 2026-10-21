//! Finite D0 near a vanishing Gram determinant, which the contour
//! representation (boxc) evaluates: two external momenta nearly parallel, the
//! invariant of the pair next to its threshold.
//!
//! The reported box (NinjaRS gg -> HHH, p² = (202500, 15625, 15625, 15625),
//! s23 = 62501.31 next to (2 x 125)², four masses 29929) lost two defects:
//!
//! - Its residue term `r0` takes `log(q)/(1 - q)` of quotients of nearly
//!   coincident poles lying on opposite lips of a cut: the sheet winding
//!   contributes `i pi n/(1 - q)`, large and cancelling against the other
//!   terms, and `1 - q` formed from the rounded quotient cost binary64 five
//!   digits (6.8e-10 relative; OneLOop 7.9e-15). The winding is added as
//!   `i pi n/(y1 - y2)` with the difference of the poles.
//! - That left 3.4e-11 to 3.6e-11 (native, SymJIT) and 7.5e-12 (interpreted).
//!   In the canonical (sorted-lambda) relabelling, two of the three `tfun`
//!   take a contour direction `(1, beta)` nearly parallel to their line:
//!   `a1 = g + beta h` is 1.4e-5 of `g` (the Gram determinant), and one pole
//!   lies within 5e-6 to 7e-5 of the root of each line edge, whose `s3` terms
//!   then took `y - z` from rounded values, rounded differently on each edge.
//!   The offsets now come from Vieta's product with the poles' own `a1`
//!   (`box_complex::line_product`): 1.9e-14 to 4.1e-14 on every route.
//!
//! References: OneLOop 3.6 in quadruple precision.
#![cfg(feature = "generated-evaluators")]
use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend,
};
use symbolica::domains::float::{Complex, Float, RealLike};

type C = Complex<f64>;

/// `[p1², p2², p3², p4², s12, s23, m0², m1², m2², m3², mu²]`.
const ARGUMENTS: [f64; 11] = [
    202500.0,
    15625.0,
    15625.0,
    15625.0,
    93305.84790585532,
    62501.31001345166,
    29929.0,
    29929.0,
    29929.0,
    29929.0,
    1e4,
];
const EXPECTED: [f64; 2] = [-1.907_030_378_722_972_01e-10, 8.834_687_875_568_863_24e-10];

/// Relative tolerances in binary64 and at 256 bits.
const BINARY64: f64 = 1e-13;
const PRECISE: f64 = 1e-15;

fn on_stack(test: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(test)
        .unwrap()
        .join()
        .unwrap();
}

/// Largest Laurent coefficient error relative to the largest reference
/// coefficient (infinite for a non-finite value).
fn error(actual: &[C], expected: &[C]) -> f64 {
    let scale = expected
        .iter()
        .fold(0.0_f64, |s, z| s.max(z.re.hypot(z.im)));
    let mut worst = 0.0_f64;
    for (a, e) in actual.iter().zip(expected) {
        if !a.re.is_finite() || !a.im.is_finite() {
            return f64::INFINITY;
        }
        worst = worst.max((a.re - e.re).hypot(a.im - e.im));
    }
    worst / scale
}

/// The binary64 routes: native code (hand-written kernels and generated
/// code), the interpreted expressions and SymJIT.
struct Binary64 {
    native: [NativeEvaluator<f64>; 2],
}

impl Binary64 {
    fn new() -> Self {
        Self {
            native: [true, false].map(|handwritten| {
                let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
                evaluator.set_handwritten(handwritten);
                evaluator
            }),
        }
    }

    fn evaluate(&mut self, input: &[C]) -> Vec<(&'static str, [C; 3])> {
        let mut values = Vec::new();
        for (route, evaluator) in ["native (hand-written)", "native (generated)"]
            .into_iter()
            .zip(&mut self.native)
        {
            let mut output = [C::new(f64::NAN, f64::NAN); 3];
            evaluator.evaluate(input, &mut output).unwrap();
            values.push((route, output));
        }
        for backend in [EvaluationBackend::Expression, EvaluationBackend::SymJit] {
            let mut output = [C::new(f64::NAN, f64::NAN); 3];
            evaluate_with_backend(ScalarIntegral::D0, input, &mut output, backend).unwrap();
            values.push((backend.name(), output));
        }
        values
    }
}

/// Native code and the interpreted expressions at `bits`.
struct Precise {
    bits: u32,
    native: NativeEvaluator<Float>,
    expression: PrecisionEvaluator,
}

impl Precise {
    fn new(bits: u32) -> Self {
        Self {
            bits,
            native: NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::D0, bits)
                .unwrap(),
            expression: PrecisionEvaluator::with_binary_precision_and_backend(
                ScalarIntegral::D0,
                bits,
                EvaluationBackend::Expression,
            )
            .unwrap(),
        }
    }

    fn evaluate(&mut self, input: &[C]) -> [(&'static str, [C; 3]); 2] {
        let bits = self.bits;
        let input: Vec<_> = input
            .iter()
            .map(|z| Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im)))
            .collect();
        let fresh = || {
            core::array::from_fn::<_, 3, _>(|_| Complex::new(Float::new(bits), Float::new(bits)))
        };
        let rounded =
            |output: [Complex<Float>; 3]| output.map(|z| C::new(z.re.to_f64(), z.im.to_f64()));
        let mut native = fresh();
        self.native.evaluate(&input, &mut native).unwrap();
        let mut expression = fresh();
        self.expression.evaluate(&input, &mut expression).unwrap();
        [
            ("native 256 bits", rounded(native)),
            ("expression 256 bits", rounded(expression)),
        ]
    }
}

/// Every evaluation route: native binary64 (hand-written kernels and generated
/// code), the interpreted expressions and SymJIT in binary64, and 256-bit
/// native code and interpreted expressions.
#[test]
fn box_near_a_vanishing_gram_determinant_matches_oneloop_on_every_route() {
    on_stack(|| {
        let input = ARGUMENTS.map(|x| C::new(x, 0.0));
        let expected = [
            C::new(EXPECTED[0], EXPECTED[1]),
            C::new(0.0, 0.0),
            C::new(0.0, 0.0),
        ];
        let mut failures = Vec::new();
        let binary64 = Binary64::new().evaluate(&input);
        let precise = Precise::new(256).evaluate(&input);
        for (route, value, tolerance) in binary64
            .into_iter()
            .map(|(route, value)| (route, value, BINARY64))
            .chain(
                precise
                    .into_iter()
                    .map(|(route, value)| (route, value, PRECISE)),
            )
        {
            let error = error(&value, &expected);
            if !(error <= tolerance) {
                failures.push(format!("{route}: {value:?} vs {EXPECTED:?} ({error:e})"));
            }
        }
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    });
}

/// A box next to an exactly vanishing Gram determinant in binary64 (p1², p2²
/// and p4² equal to 2e-16, s12 and s23 2e-6 apart), where the default contour
/// of the first `tfun` has `K = b (g + beta' h) = 0` exactly: its chart
/// divides by K. Native binary64 code (hand-written and generated) returned
/// NaN here, and every binary64 route (interpreted and SymJIT too) in 28% to
/// 40% of its 1-ulp perturbations. That chart now puts its far pole 2^52
/// beyond the near one in binary64 (`box_complex::chart_k`; the second
/// contour, taken there before, has its `a1` at the rounding level and was
/// wrong by O(1) on other such boxes, tests/data/contour_gram_rounding.txt),
/// and keeps it: moved to K's rounding level, the chart keeps 13 bits.
/// OneLOop 3.6 in quadruple precision; its double precision keeps 6e-14.
#[test]
fn box_with_a_vanishing_contour_coefficient_is_finite_on_every_route() {
    on_stack(|| {
        let input = [
            2.430909632840068,
            2.4309096328400637,
            144.9999453227875,
            2.4309096328400654,
            71.2844473512537,
            71.2845883386827,
        ]
        .map(|x| C::new(x, 0.0))
        .into_iter()
        .chain([
            C::new(3.2296223697605373, 0.0),
            C::new(3.5171272384943344, -0.007868476275933284),
            C::new(2.436298702662891, -0.00013942610297326737),
            C::new(1.5947319955211854, 0.0),
            C::new(1.0, 0.0),
        ])
        .collect::<Vec<_>>();
        let expected = [
            C::new(-6.78208635248455749e-4, -2.92577907244204226e-3),
            C::new(0.0, 0.0),
            C::new(0.0, 0.0),
        ];
        let mut failures = Vec::new();
        for (route, value) in Binary64::new().evaluate(&input) {
            let error = error(&value, &expected);
            if !(error <= 1e-13) {
                failures.push(format!("{route}: {value:?} ({error:e})"));
            }
        }
        for (route, value) in Precise::new(256).evaluate(&input) {
            let error = error(&value, &expected);
            if !(error <= PRECISE) {
                failures.push(format!("{route}: {value:?} ({error:e})"));
            }
        }
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    });
}

/// `kind mu² p.. (re im)..` fixture rows with references, as evaluator inputs
/// (mu² last) and the three reference coefficients.
fn near_threshold_rows() -> Vec<(usize, Vec<C>, [C; 3])> {
    include_str!("data/near_threshold_boxes.txt")
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.starts_with('#') && !line.trim().is_empty())
        .map(|(index, line)| {
            let v: Vec<f64> = line
                .split_whitespace()
                .map(|x| x.parse().unwrap())
                .collect();
            assert_eq!(v.len(), 22, "line {}", index + 1);
            assert_eq!(v[0], 4.0);
            let mut input: Vec<C> = v[2..8].iter().map(|&x| C::new(x, 0.0)).collect();
            input.extend(v[8..16].chunks_exact(2).map(|m| C::new(m[0], m[1])));
            input.push(C::new(v[1], 0.0));
            let expected = core::array::from_fn(|i| C::new(v[16 + 2 * i], v[17 + 2 * i]));
            (index + 1, input, expected)
        })
        .collect()
}

/// Binary64 tolerance of the near-threshold scan, relative to the largest
/// reference coefficient. Every route keeps 3.9e-12 or better on every row
/// but line 73 (see [NEAR_THRESHOLD_RESIDUAL]). Before the poles' offsets
/// from the line edges' roots (`box_complex::line_product`), 34 of the 200
/// boxes were beyond it on some route and 22 beyond 1e-9 (up to 5e-3 native,
/// 2e-2 interpreted and SymJIT); OneLOop's double precision loses more than
/// 5e-12 on 32 and more than 1e-9 on 16.
const NEAR_THRESHOLD: f64 = 5e-12;

/// Line 73 (three masses and a massless line, the pair's Källén function
/// 1.5e-11 of the largest invariant squared) keeps
/// 1.35e-11 interpreted (09b2c74 1.16e-10) to 3.2e-10 hand-written
/// (OneLOop's double precision 5.5e-11). There the canonical relabelling's first `tfun` has a nearly
/// double contour direction `beta` and the other two `tfun` nearly
/// coincident poles next to the end x = 0 of an edge, whose quadratic `s3`
/// terms (5e3 each, before the division by K = 5e13) cancel to 3e2 and take
/// the poles' rounding with them.
const NEAR_THRESHOLD_RESIDUAL: (usize, f64) = (73, 1e-9);

/// 200 random physical boxes next to a vanishing Gram determinant
/// (tests/data/near_threshold_boxes.txt, from
/// tests/generate_near_threshold_rows.py): two external momenta nearly
/// parallel, their invariant 1e-12 to 1e-2 (relative) above its threshold;
/// equal, random, complex masses and a massless line; scales 1e-2 to 1e4. Every
/// binary64 route against OneLOop 3.6 in quadruple precision, and native code
/// and the interpreted expressions at 256 bits (within 4e-31 here).
#[test]
fn near_threshold_scan_matches_oneloop_on_every_route() {
    on_stack(|| {
        let rows = near_threshold_rows();
        assert_eq!(rows.len(), 200);
        let mut binary64 = Binary64::new();
        let mut precise = Precise::new(256);
        let mut failures = Vec::new();
        for (line, input, expected) in &rows {
            let tolerance = match NEAR_THRESHOLD_RESIDUAL {
                (residual, tolerance) if residual == *line => tolerance,
                _ => NEAR_THRESHOLD,
            };
            for (route, value) in binary64.evaluate(input) {
                let error = error(&value, expected);
                if !(error <= tolerance) {
                    failures.push(format!(
                        "line {line} {route}: {value:?} vs {expected:?} ({error:.2e})"
                    ));
                }
            }
            for (route, value) in precise.evaluate(input) {
                let error = error(&value, expected);
                if !(error <= 1e-14) {
                    failures.push(format!(
                        "line {line} {route}: {value:?} vs {expected:?} ({error:.2e})"
                    ));
                }
            }
        }
        assert!(
            failures.is_empty(),
            "{} failures:\n{}",
            failures.len(),
            failures.join("\n")
        );
    });
}
