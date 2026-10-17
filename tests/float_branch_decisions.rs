//! Exact branch decisions of the interpreted expressions in arbitrary
//! precision (`Complex<Float>`, as NinjaRS's Symbolica masters evaluate them).
use oneloop::OneLoopExpressions;
use symbolica::{
    atom::Atom,
    domains::{
        float::{Complex, Float, FloatLike, RealLike, SingleFloat},
        rational::Rational,
    },
    evaluate::ExpressionEvaluator,
};

type C = Complex<Float>;

fn c0_evaluator() -> ExpressionEvaluator<Complex<Rational>> {
    let x: Vec<Atom> = (0..7)
        .map(|i| symbolica::symbol!(format!("float_branch::arg_{i}")).to_atom())
        .collect();
    let e = OneLoopExpressions::new();
    let series = e.c0([&x[0], &x[1], &x[2]], [&x[3], &x[4], &x[5]], &x[6]);
    e.evaluator(series.coefficients(), &x).unwrap()
}

fn real(bits: u32, n: i64) -> C {
    let v = Float::new(bits).from_rational(&Rational::from(n));
    Complex::new(v.clone(), v.zero())
}

/// C0(-4000, 0, 0; 29929 x3; mu² = 10^4) with exact inputs. With two vanishing
/// invariants and equal masses the trif3 coefficient `c` is analytically zero;
/// before the exact zero-momentum r-function roots, Float evaluated it as
/// `173 (1/173) - 1`, a few-bit residue that selected the quadratic formula
/// and returned NaN at 116 bits and a 5-bit value at 108 bits.
#[test]
fn equal_mass_triangle_with_two_vanishing_invariants_at_every_precision() {
    std::thread::Builder::new()
        .stack_size(512 << 20)
        .spawn(|| {
            let evaluator = c0_evaluator();
            let expected = -1.652_338_593_294_892_214_409_220_250_488_900_9e-5_f64;
            let mut failures = Vec::new();
            for bits in (100..=260).step_by(4).chain([108, 116]) {
                let p = Float::new(bits);
                let mut e = evaluator.clone().map_coeff_with_prec(
                    &|v: &Complex<Rational>| {
                        Complex::new(p.from_rational(&v.re), p.from_rational(&v.im))
                    },
                    bits,
                );
                let args = [-4000, 0, 0, 29929, 29929, 29929, 10000].map(|v| real(bits, v));
                let mut out = vec![args[0].zero(); 3];
                e.evaluate(&args, &mut out);
                let mut values = vec![("expression", out[0].clone())];
                #[cfg(feature = "generated-evaluators")]
                {
                    let mut native = oneloop::NativeEvaluator::<Float>::with_binary_precision(
                        oneloop::ScalarIntegral::C0,
                        bits,
                    )
                    .unwrap();
                    let mut out = vec![args[0].zero(); 3];
                    native.evaluate(&args, &mut out).unwrap();
                    values.push(("native", out[0].clone()));
                }
                for (route, value) in values {
                    let error = (value.re.to_f64() - expected).hypot(value.im.to_f64());
                    if !(error <= 1e-15 * expected.abs()) || value.re.get_precision() < bits - 16 {
                        failures.push(format!(
                            "{route} {bits} bits: {} {} ({} tracked bits)",
                            value.re,
                            value.im,
                            value.re.get_precision()
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

fn complex(bits: u32, re: f64, im: f64) -> C {
    Complex::new(Float::with_val(bits, re), Float::with_val(bits, im))
}

/// Two-mass C0 with a vanishing invariant between unequal masses, one of them
/// complex (NinjaRS syn-n6-r6): the coefficient `c` of trif2 is analytically
/// zero there. Computed as `s4/s2 - r24` it was a few-bit residue in `Float`,
/// the second root became noise and the root term divided by it, so that
/// results from 248 bits up (to 2048 interpreted) had 2-8 tracked bits and
/// wrong leading digits. References: OneLOop 3.6 (double and quadruple).
#[test]
fn two_mass_triangle_with_a_vanishing_invariant_at_every_precision() {
    std::thread::Builder::new()
        .stack_size(512 << 20)
        .spawn(|| {
            let evaluator = c0_evaluator();
            let calls = [
                (
                    [58911.168392294385, 0.0, -236743.32762192388],
                    (-2.523_145_603_457_950_70e-5, -8.171_478_148_049_019_03e-6),
                ),
                (
                    [70987.92439374398, 0.0, -110372.46037321811],
                    (-3.420_847_151_881_713_67e-5, -1.901_034_314_146_606_00e-5),
                ),
            ];
            let mut failures = Vec::new();
            for bits in [128, 248, 256, 264, 300, 512, 1024, 2048] {
                let p = Float::new(bits);
                let mut interpreted = evaluator.clone().map_coeff_with_prec(
                    &|v: &Complex<Rational>| {
                        Complex::new(p.from_rational(&v.re), p.from_rational(&v.im))
                    },
                    bits,
                );
                #[cfg(feature = "generated-evaluators")]
                let mut native = oneloop::NativeEvaluator::<Float>::with_binary_precision(
                    oneloop::ScalarIntegral::C0,
                    bits,
                )
                .unwrap();
                for (momenta, expected) in calls {
                    let args = vec![
                        complex(bits, momenta[0], 0.0),
                        complex(bits, momenta[1], 0.0),
                        complex(bits, momenta[2], 0.0),
                        complex(bits, 0.0, 0.0),
                        complex(bits, 29929.0, 0.0),
                        complex(bits, 6464.160000000001, -168.84000000000003),
                        complex(bits, 1e4, 0.0),
                    ];
                    let mut values = Vec::new();
                    let mut out = vec![args[0].zero(); 3];
                    interpreted.evaluate(&args, &mut out);
                    values.push(("expression", out[0].clone()));
                    #[cfg(feature = "generated-evaluators")]
                    {
                        let mut out = vec![args[0].zero(); 3];
                        native.evaluate(&args, &mut out).unwrap();
                        values.push(("native", out[0].clone()));
                    }
                    for (route, value) in values {
                        let error = (value.re.to_f64() - expected.0)
                            .hypot(value.im.to_f64() - expected.1)
                            / expected.0.hypot(expected.1);
                        if !(error <= 1e-15) || value.re.get_precision() < 64 {
                            failures.push(format!(
                                "{route} {bits} bits {momenta:?}: {} {} ({} tracked bits)",
                                value.re,
                                value.im,
                                value.re.get_precision()
                            ));
                        }
                    }
                }
            }
            assert!(failures.is_empty(), "{}", failures.join("\n"));
        })
        .unwrap()
        .join()
        .unwrap();
}
