//! Exact branch decisions of the interpreted expressions in arbitrary
//! precision (`Complex<Float>`, as NinjaRS's Symbolica masters evaluate them).
use oneloop::OneLoopExpressions;
use symbolica::{
    atom::{Atom, AtomCore},
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
