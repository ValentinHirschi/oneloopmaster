//! Exact one-mass IR boxes whose equivalent denominator chart reaches box09.
//!
//! The original second row returned NaN and a zero double pole in eight of
//! its 24 labellings. The box10 chart exchange made p4=m, but did not repeat
//! the infrared classification. References are unchanged OneLOop3.6 real128
//! values from non-error orderings, independently checked by Collier1.2.9
//! COLI and DD (all flags zero; finite conversion +pi²/6*double). OneLOop's
//! original order emits an explicit error: its returned zeros are not targets.
#![cfg(feature = "generated-evaluators")]
use oneloop::{NativeEvaluator, NativeFloat, OneLoopExpressions, ScalarIntegral};
use symbolica::{
    domains::float::{Complex, DoubleFloat, Float, Real, SingleFloat},
    evaluate::EvaluationDomain,
};

const BITS: u32 = 256;
const REFERENCES: [[[&str; 2]; 3]; 2] = [
    [
        [
            "-2.60132556713013326849351586964136613E-13",
            "3.43770874061117486475847201381883162E-13",
        ],
        [
            "2.71393732633611655399202998979820782E-13",
            "0.00000000000000000000000000000000000E+00",
        ],
        [
            "-1.91999999999999999999999999999999997E-13",
            "-0.00000000000000000000000000000000000E+00",
        ],
    ],
    [
        [
            "2.98569165331053244434628966883798097E-13",
            "1.43570153780154044618338274416424854E-12",
        ],
        [
            "7.30683574087754623175189149732696656E-13",
            "-1.60849543863797413809287341223910543E-12",
        ],
        [
            "-2.55999999999999999999999999999999997E-13",
            "-0.00000000000000000000000000000000000E+00",
        ],
    ],
];
// The exact pairwise shift invariants and squared masses, independent of any
// chosen cyclic ordering. All input integers are exactly representable.
const S: [[i64; 5]; 5] = [
    [0, 0, 0, -937500, -1000000],
    [0, 0, 6250000, 1562500, 250000],
    [0, 6250000, 0, 0, 2593750],
    [-937500, 1562500, 0, 0, 250000],
    [-1000000, 250000, 2593750, 250000, 0],
];
const MASSES: [i64; 5] = [0, 0, 0, 0, 250000];

fn permutations(subset: [usize; 4]) -> Vec<[usize; 4]> {
    let mut result = Vec::new();
    for a in subset {
        for b in subset {
            for c in subset {
                for d in subset {
                    if a != b && a != c && a != d && b != c && b != d && c != d {
                        result.push([a, b, c, d]);
                    }
                }
            }
        }
    }
    assert_eq!(result.len(), 24);
    result
}

fn close<T: NativeFloat>(
    values: &[Complex<T>],
    expected: &[[&str; 2]; 3],
    scale: i64,
    convert: &impl Fn(&T) -> Float,
    tolerance: &str,
    label: &str,
) {
    let tolerance = Float::parse(tolerance, Some(BITS)).unwrap();
    for (coefficient, (value, expected)) in values.iter().zip(expected).enumerate() {
        let expected = expected
            .map(|s| Float::parse(s, Some(BITS)).unwrap() / Float::with_val(BITS, scale * scale));
        let norm = (expected[0].clone() * &expected[0] + expected[1].clone() * &expected[1]).sqrt();
        assert!(!norm.is_zero(), "these three coefficients are all nonzero");
        for (part, (actual, expected)) in
            [&value.re, &value.im].into_iter().zip(expected).enumerate()
        {
            assert!(
                actual.is_finite(),
                "{label}, coefficient{coefficient},part{part}: nonfinite"
            );
            // A zero imaginary part uses the same coefficient's nonzero norm,
            // not a unit floor that could hide a lost 1e-13 physical result.
            let error = (convert(actual) - expected).norm();
            assert!(
                error <= norm.clone() * &tolerance,
                "{label},coefficient{coefficient},part{part}: error={error:e},norm={norm:e}"
            );
        }
    }
}

fn check<T: NativeFloat>(prototype: T, bits: u32, convert: impl Fn(&T) -> Float, tolerance: &str)
where
    Complex<T>: EvaluationDomain,
{
    let context = OneLoopExpressions::new();
    let variables = (0..11)
        .map(|i| symbolica::symbol!(format!("box10_chart::x{i}")).to_atom())
        .collect::<Vec<_>>();
    let series = context.d0(
        variables[..6]
            .iter()
            .collect::<Vec<_>>()
            .try_into()
            .unwrap(),
        variables[6..10]
            .iter()
            .collect::<Vec<_>>()
            .try_into()
            .unwrap(),
        &variables[10],
    );
    let mut expression = context
        .evaluator(series.coefficients(), &variables)
        .unwrap()
        .map_coeff_with_prec(
            &|value| {
                Complex::new(
                    prototype.from_rational(&value.re),
                    prototype.from_rational(&value.im),
                )
            },
            bits,
        );
    let mut native = NativeEvaluator::<T>::with_binary_precision(ScalarIntegral::D0, bits).unwrap();
    for (case, subset) in [[0, 1, 2, 4], [0, 1, 3, 4]].into_iter().enumerate() {
        for order in permutations(subset) {
            let [a, b, c, d] = order;
            for scale in [1, 4] {
                let row = [
                    S[a][b], S[b][c], S[c][d], S[d][a], S[a][c], S[b][d], MASSES[a], MASSES[b],
                    MASSES[c], MASSES[d], 1000000,
                ];
                let input =
                    row.map(|v| Complex::new(prototype.from_i64(v * scale), prototype.zero()));
                for generated in [false, true] {
                    let mut output = core::array::from_fn::<_, 3, _>(|_| {
                        Complex::new(prototype.zero(), prototype.zero())
                    });
                    if generated {
                        native.evaluate(&input, &mut output).unwrap();
                    } else {
                        expression.evaluate(&input, &mut output);
                    }
                    close(
                        &output,
                        &REFERENCES[case],
                        scale,
                        &convert,
                        tolerance,
                        &format!("bits{bits},native{generated},order{order:?},scale{scale}"),
                    );
                }
            }
        }
    }
}

#[test]
fn chart_exchange_reclassifies_the_exact_mass_shell_at_every_precision() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            check(0.0_f64, 53, |x| Float::with_val(BITS, *x), "2e-12");
            check(
                DoubleFloat::from(0.0),
                106,
                |x| {
                    let pair = x.into_inner();
                    Float::with_val(BITS, pair.hi()) + Float::with_val(BITS, pair.lo())
                },
                "1e-27",
            );
            for bits in [128, 256] {
                check(Float::new(bits), bits, Clone::clone, "1e-29");
            }
        })
        .unwrap()
        .join()
        .unwrap();
}
