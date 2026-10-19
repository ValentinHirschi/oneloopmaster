//! A sheet sign must use one evaluated argument, including tracked Float inputs.
use oneloop::OneLoopExpressions;
use symbolica::{
    atom::{Atom, AtomCore},
    domains::{
        float::{Complex, DoubleFloat, Float, FloatLike, Real, SingleFloat},
        rational::Rational,
    },
};

type C = Complex<Float>;
const BITS: u32 = 266;

fn real(value: &str) -> C {
    C::new(Float::parse(value, Some(BITS)).unwrap(), Float::new(BITS))
}

fn masses() -> [C; 2] {
    let mass = real("91.188");
    let pole = (&mass * &mass).sqrt();
    let squared = &pole * &pole;
    let reparsed = real(&format!("{:.100e}", squared.re));
    assert_eq!(squared, reparsed);
    assert_eq!(squared.re.get_precision(), BITS + 1);
    assert_eq!(reparsed.re.get_precision(), BITS);
    [squared, reparsed]
}

fn on_stack(test: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(test)
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn compound_discriminant_keeps_its_real_sign_at_mixed_precision() {
    on_stack(|| {
        let context = OneLoopExpressions::new();
        let x = symbolica::symbol!("mixed_sheet::m0").to_atom();
        let y = symbolica::symbol!("mixed_sheet::m1").to_atom();
        let sqrt = symbolica::get_symbol!("oneloop::__olo_sqrt_lower").unwrap();
        let sign = symbolica::get_symbol!("oneloop::__olo_sign_nonnegative").unwrap();
        let k = (&x + &y - 1000000) / (sqrt.call((&x,)) * sqrt.call((&y,)));
        let discriminant = k.pow(2) - 4;
        let outputs = [
            discriminant.clone(),
            sign.call((&discriminant,)),
            sign.call((-&discriminant,)),
            sign.call((Atom::num(0),)),
        ];
        let exact = context.evaluator(&outputs, &[x, y]).unwrap();
        macro_rules! check_fixed {
            ($prototype:expr) => {{
                let prototype = $prototype;
                let mut evaluator = exact.clone().map_coeff(&|value| {
                    Complex::new(
                        prototype.from_rational(&value.re),
                        prototype.from_rational(&value.im),
                    )
                });
                let mass = Complex::new(
                    prototype.from_rational(&Rational::from((8315251344i64, 1000000i64))),
                    prototype.zero(),
                );
                let zero = Complex::new(prototype.zero(), prototype.zero());
                let mut result = [zero; 4];
                evaluator.evaluate(&[mass, mass], &mut result);
                for (value, expected) in result[1..].iter().zip([1, -1, 1]) {
                    assert_eq!(
                        *value,
                        Complex::new(prototype.from_i64(expected), prototype.zero())
                    );
                }
            }};
        }
        check_fixed!(0.0_f64);
        check_fixed!(DoubleFloat::from(0.0));
        let prototype = Float::new(BITS);
        let mut evaluator = exact.map_coeff_with_prec(
            &|value| {
                C::new(
                    prototype.from_rational(&value.re),
                    prototype.from_rational(&value.im),
                )
            },
            BITS,
        );
        let bound = symbolica::symbol!("mixed_sheet::bound_sign_argument").to_atom();
        let mut bound_sign = context
            .evaluator(&[sign.call((&bound,))], &[bound])
            .unwrap()
            .map_coeff_with_prec(
                &|value| {
                    C::new(
                        prototype.from_rational(&value.re),
                        prototype.from_rational(&value.im),
                    )
                },
                BITS,
            );
        for mass in masses() {
            let input_bits = mass.re.get_precision();
            let mut result = std::array::from_fn::<_, 4, _>(|_| real("0"));
            evaluator.evaluate(&[mass.clone(), mass], &mut result);
            assert!(result[0].re > real("10000").re);
            assert_eq!(result[1], real("1"));
            assert_eq!(result[2], real("-1"));
            assert_eq!(result[3], real("1"));
            // Unary native negation retains the stronger precision. A symbolic
            // -discriminant would multiply by the prepared -1 before entry.
            assert_eq!(result[0].re.get_precision(), input_bits);
            let negative = -result[0].clone();
            assert_eq!(negative.re.get_precision(), input_bits);
            for (argument, expected) in [
                (result[0].clone(), real("1")),
                (negative, real("-1")),
                (real("0"), real("1")),
            ] {
                let mut output = [real("0")];
                bound_sign.evaluate(&[argument], &mut output);
                assert_eq!(output[0], expected);
            }
        }
    });
}

#[test]
fn equal_masses_with_different_tracked_precision_keep_the_feynman_sheet() {
    on_stack(|| {
        let context = OneLoopExpressions::new();
        let arguments = (0..11)
            .map(|i| symbolica::symbol!(format!("mixed_sheet::arg_{i}")).to_atom())
            .collect::<Vec<_>>();
        let series = context.d0(
            arguments[..6]
                .iter()
                .collect::<Vec<_>>()
                .try_into()
                .unwrap(),
            arguments[6..10]
                .iter()
                .collect::<Vec<_>>()
                .try_into()
                .unwrap(),
            &arguments[10],
        );
        let prototype = Float::new(BITS);
        let mut evaluator = context
            .evaluator(series.coefficients(), &arguments)
            .unwrap()
            .map_coeff_with_prec(
                &|value| {
                    C::new(
                        prototype.from_rational(&value.re),
                        prototype.from_rational(&value.im),
                    )
                },
                BITS,
            );
        // Independent OneLOop3.6 real128 result, with mu_r=1000 in its interface.
        let expected = C::new(
            real("-3.43176590414174697596155714935708795e-10").re,
            real("4.62278850551111554358291148483041637e-10").re,
        );
        let mut results = Vec::new();
        for mass in masses() {
            let row = vec![
                real("0"),
                real("1000000"),
                real("0"),
                real("-100000"),
                real("0"),
                real("0"),
                real("0"),
                mass.clone(),
                mass,
                real("0"),
                real("1000000"),
            ];
            let mut output = std::array::from_fn::<_, 3, _>(|_| real("0"));
            evaluator.evaluate(&row, &mut output);
            assert!(
                (&output[0] - &expected).norm().re < real("1e-40").re,
                "D0 must match the independently established Feynman sheet: {output:?}"
            );
            assert_eq!(output[1], real("0"));
            assert_eq!(output[2], real("0"));
            #[cfg(feature = "generated-evaluators")]
            {
                let mut native = oneloop::NativeEvaluator::<Float>::with_binary_precision(
                    oneloop::ScalarIntegral::D0,
                    BITS,
                )
                .unwrap();
                let mut native_output = std::array::from_fn::<_, 3, _>(|_| real("0"));
                native.evaluate(&row, &mut native_output).unwrap();
                assert!((&output[0] - &native_output[0]).norm().re < real("1e-80").re);
            }
            results.push(output);
        }
        assert!((&results[0][0] - &results[1][0]).norm().re < real("1e-80").re);
    });
}
