//! A four-mass D0 with a tiny but nonzero invariant (NinjaRS, a light-like
//! photon invariant rounded to 2^-114 at 128 bits): p = [eps, 40000, 0,
//! -100000/3, eps, 0] with eps = 2^-114, m² = 29929 x4, mu² = 10^4, exact
//! inputs. After the channel exchange its externals are (eps, 0, 0, eps), and
//! the contour box's zero-external chart formed coefficients such as
//! `m2 - m3 - p2 = -eps` and `a + b + c = eps` by cancellation. The values were
//! exact, but `Float` tracks a difference with only the bits left after the
//! cancellation and computes everything downstream at that precision: from
//! 132 to 184 bits the result kept 7-55 bits and was wrong beyond them (6e-2
//! at 136 bits). Mass differences are now subtracted first and the coefficient
//! sums that cancel use the invariants they equal (`box_complex::in_order`,
//! `coefficient_sum`). Reference: 512 bits (383 tracked bits), which agrees
//! with the four-mass Denner-Nierste-Scharf representation at 512 bits and,
//! to 2^-129, with the eps -> 0 limit.
#![cfg(feature = "generated-evaluators")]
use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend,
};
use symbolica::domains::{
    float::{Complex, Float, FloatLike, RealLike, SingleFloat},
    rational::Rational,
};

/// 100 digits of the 512-bit value (imaginary part zero).
const EXPECTED: &str = "1.970543895668874901801963663849211004697908797491396868467915495952052342393457601063637686226326854e-10";

fn arguments(bits: u32) -> Vec<Complex<Float>> {
    let prototype = Float::new(bits);
    let real = |value: Rational| {
        let value = prototype.from_rational(&value);
        Complex::new(value.clone(), value.zero())
    };
    let eps = Rational::from((1, 2)).pow(114);
    let integer = |n: i64| Rational::from(n);
    vec![
        real(eps.clone()),
        real(integer(40000)),
        real(integer(0)),
        real(Rational::from((-100000, 3))),
        real(eps),
        real(integer(0)),
        real(integer(29929)),
        real(integer(29929)),
        real(integer(29929)),
        real(integer(29929)),
        real(integer(10000)),
    ]
}

/// Native code and interpreted expressions in `Float` from 128 to 1024 bits:
/// up to 256 bits at most 16 bits are lost (the chart's own cancellation costs
/// about 13, as in binary64), and the tracked precision is within 10 bits of
/// the actual error (Numerica's tracking is approximate; it was off by 2-4 bits
/// and the value off by up to 6e-2 before).
#[test]
fn box_with_a_tiny_invariant_is_accurate_at_every_precision() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let mut failures = Vec::new();
            let precisions = (128..=200).step_by(4).chain([256, 384, 512, 1024]);
            for bits in precisions {
                let input = arguments(bits);
                let mut native =
                    NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::D0, bits)
                        .unwrap();
                let mut expression = PrecisionEvaluator::with_binary_precision_and_backend(
                    ScalarIntegral::D0,
                    bits,
                    EvaluationBackend::Expression,
                )
                .unwrap();
                for route in ["native", "expression"] {
                    let mut output = core::array::from_fn::<_, 3, _>(|_| {
                        Complex::new(Float::new(bits), Float::new(bits))
                    });
                    if route == "native" {
                        native.evaluate(&input, &mut output).unwrap();
                    } else {
                        expression.evaluate(&input, &mut output).unwrap();
                    }
                    let value = &output[0];
                    let expected = Float::parse(EXPECTED, Some(bits.max(384))).unwrap();
                    let scale = expected.to_f64();
                    let error = (value.re.clone() - &expected)
                        .to_f64()
                        .hypot(value.im.to_f64())
                        / scale;
                    let tracked = value.re.get_precision().min(bits) as i32;
                    // The reference itself has 100 digits.
                    let honest = error <= 2f64.powi(10 - tracked).max(1e-99);
                    let accurate = bits > 256 || error <= 2f64.powi(16 - bits as i32);
                    if !honest || !accurate {
                        failures.push(format!(
                            "{route} {bits} bits: {} {} ({tracked} tracked bits, error {error:e})",
                            value.re, value.im
                        ));
                    }
                }
            }
            assert!(failures.is_empty(), "{}", failures.join("\n"));
        })
        .unwrap()
        .join()
        .unwrap();
}

/// Binary64, where eps is far below the resolution of the other invariants:
/// every route stays within 1e-12 (the chart cancels by about 10^3 there).
#[test]
fn box_with_a_tiny_invariant_in_binary64() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let eps = 0.5f64.powi(114);
            let input = [eps, 40000.0, 0.0, -100000.0 / 3.0, eps, 0.0]
                .into_iter()
                .chain([29929.0; 4])
                .chain([1e4])
                .map(|x| Complex::new(x, 0.0))
                .collect::<Vec<_>>();
            let expected = 1.970_543_895_668_874_901_801_963_663_849_211e-10_f64;
            let mut failures = Vec::new();
            let mut check = |route: String, value: Complex<f64>| {
                let error = (value.re - expected).hypot(value.im) / expected;
                if !(error <= 1e-12) {
                    failures.push(format!("{route}: {value:?} ({error:e})"));
                }
            };
            for handwritten in [true, false] {
                let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
                evaluator.set_handwritten(handwritten);
                let mut output = [Complex::new(0.0, 0.0); 3];
                evaluator.evaluate(&input, &mut output).unwrap();
                check(format!("native (hand-written: {handwritten})"), output[0]);
            }
            for backend in [EvaluationBackend::Expression, EvaluationBackend::SymJit] {
                let mut output = [Complex::new(0.0, 0.0); 3];
                evaluate_with_backend(ScalarIntegral::D0, &input, &mut output, backend).unwrap();
                check(backend.name().to_string(), output[0]);
            }
            assert!(failures.is_empty(), "{}", failures.join("\n"));
        })
        .unwrap()
        .join()
        .unwrap();
}
