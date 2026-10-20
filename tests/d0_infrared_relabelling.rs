//! Soft-divergent boxes in all 24 labellings of their propagators.
//!
//! The scalar box depends only on its set of propagators, so the 24
//! relabellings (invariants permuted with the lines) are one integral. OneLOop's
//! dispatch finds the soft limits on the four external invariants only: a
//! massless line whose neighbours are on their mass shells,
//! (q_i - q_k)^2 = m_k. The 16 labellings that are not rotations or
//! reflections of a cyclic one put one such neighbour pair on a channel (s12
//! or s23) instead, and the finite formulas were taken: binary64 returned a
//! finite value without the pole, DoubleFloat and Float a wrong finite value
//! or NaN (e+ e- -> t t~ h, pyAmpliCol). Every numerical D0 entry point now
//! evaluates these points in the labelling of the same box with both pairs on
//! the legs, also instead of the vanishing-channel choice when both apply.
//! The symbolic expressions do not relabel: their numerical users apply
//! `oneloop::regular_d0_labelling` first, as the transparent-expression routes
//! below do.
//!
//! Cases: the e+ e- -> t t~ h box (m, m, 0, m) of pyAmpliCol (its exact
//! rational invariants rounded to binary64), three massive lines with
//! different masses, the same with a vanishing invariant between two massive
//! lines (8 labellings with a vanishing channel), two opposite massless lines
//! with both or one of them soft, and two adjacent massless lines with a
//! collinear leg and both soft (its double pole; PR #10 covered its
//! labellings, now relabelled up front). Every route: hand-written, generated
//! and interpreted binary64, SymJIT, DoubleFloat (generated and interpreted)
//! and Float at 128 and 256 bits (generated and interpreted), and the public
//! D0 master's numeric callbacks in binary64, DoubleFloat and Float. References:
//! OneLOop 3.6 in quadruple precision for the identical binary64 inputs in
//! the cyclic labelling 0123, with all three Laurent coefficients [finite,
//! 1/eps, 1/eps^2]. COLLIER 1.2.9 in quadruple precision (COLI and DD,
//! DeltaIR1 and DeltaIR2 = 0, 1, finite part + pi^2/6 times the double pole)
//! agrees in all 24 labellings (DD to 1e-18 and below, COLI to 5e-17, 2.4e-16
//! for the double-pole case).
#![cfg(feature = "generated-evaluators")]
use oneloop::{
    EvaluationBackend, NativeEvaluator, NativeFloat, OneLoopExpressions, PrecisionEvaluator,
    ScalarIntegral, evaluate_with_backend,
};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::float::{Complex, DoubleFloat, Float, Real, SingleFloat},
    domains::rational::Rational,
    evaluate::{EvaluationDomain, ExpressionEvaluator},
};

const BITS: u32 = 256;
const MU_SQUARED: f64 = 1e6;

struct Case {
    name: &'static str,
    /// [p1², p2², p3², p4², s12, s23] of the labelling 0123.
    invariants: [f64; 6],
    masses: [f64; 4],
    /// OneLOop quadruple precision, [finite, 1/eps, 1/eps²] as [re, im].
    reference: [[&'static str; 2]; 3],
}

const CASES: [Case; 6] = [
    Case {
        // gamma* -> H t t~: legs (H, t, t~, gamma*), lines (t, t, g, t). The
        // inputs 1299773111931/1562500, 1756052/5 and 976955119487/3125000
        // rounded to binary64.
        name: "tth",
        invariants: [
            15625.0,
            29929.0,
            29929.0,
            831854.79163584,
            351210.4,
            312625.63823584,
        ],
        masses: [29929.0, 29929.0, 0.0, 29929.0],
        reference: [
            [
                "-2.06880701099902689869386588799378576E-10",
                "4.12462200360029074629706555698828346E-11",
            ],
            [
                "-2.68609520425235211222314830975903066E-11",
                "3.98175954777224141582499983994458669E-11",
            ],
            ["0", "0"],
        ],
    },
    Case {
        name: "three-masses",
        invariants: [20000.0, 29929.0, 6464.0, 300000.0, 90000.0, 500000.0],
        masses: [8315.0, 29929.0, 0.0, 6464.0],
        reference: [
            [
                "-6.59506278274097905169153338634493660E-10",
                "4.40299782438432406530167571520429442E-11",
            ],
            [
                "-9.27369949332159886784082971613908701E-11",
                "8.31076127376110490777699176942716560E-11",
            ],
            ["0", "0"],
        ],
    },
    Case {
        // Lines 2 and 4 on shell next to the massless line 1, and p3² = 0
        // between the massive lines 3 and 4: in the 8 labellings where that
        // is a channel, the vanishing-channel choice took a labelling with an
        // on-shell pair of line 1 on the other channel.
        name: "three-masses-vanishing-invariant",
        invariants: [29929.0, -20000.0, 0.0, 6464.0, 2000.0, 100000.0],
        masses: [0.0, 29929.0, 8315.0, 6464.0],
        reference: [
            [
                "4.59583269127477715652231483581665722E-8",
                "-5.86249578599173189390604877698850303E-8",
            ],
            [
                "4.06536678549715165979961602076657314E-9",
                "-8.69699974176854477071953708593246945E-9",
            ],
            ["0", "0"],
        ],
    },
    Case {
        name: "two-opposite-both-soft",
        invariants: [29929.0, 29929.0, 6464.0, 6464.0, -40000.0, 250000.0],
        masses: [0.0, 29929.0, 0.0, 6464.0],
        reference: [
            [
                "2.07257961826532165773595707925910221E-9",
                "-2.38738794491262306188797727655146298E-9",
            ],
            [
                "6.43883060742234429855217116189433387E-10",
                "-7.41683766260345388319068427011178834E-10",
            ],
            ["0", "0"],
        ],
    },
    Case {
        name: "two-opposite-one-soft",
        invariants: [29929.0, 50000.0, 70000.0, 6464.0, -40000.0, 250000.0],
        masses: [0.0, 29929.0, 0.0, 6464.0],
        reference: [
            [
                "-3.79345784572392332540016918657976688E-10",
                "-2.65604626962067516425687433266802989E-9",
            ],
            [
                "3.21941530371117214927608558094716693E-10",
                "-3.70841883130172694159534213505589372E-10",
            ],
            ["0", "0"],
        ],
    },
    Case {
        name: "two-adjacent-collinear",
        invariants: [0.0, 29929.0, 15625.0, 8315.0, -50000.0, 200000.0],
        masses: [0.0, 0.0, 29929.0, 8315.0],
        reference: [
            [
                "3.58133283511240650923103960055080454E-10",
                "-3.16658128438326972875841919222979207E-10",
            ],
            [
                "-1.90635644320774957798385212659172227E-12",
                "-2.05048862310757747605510114244932087E-10",
            ],
            ["-6.52690800242530633160760939347547896E-11", "0"],
        ],
    },
];

/// The exact rationals of the tth case (labelling 0123), as pyAmpliCol sent them.
const TTH_EXACT: [(i64, i64); 6] = [
    (15625, 1),
    (29929, 1),
    (29929, 1),
    (1299773111931, 1562500),
    (1756052, 5),
    (976955119487, 3125000),
];

/// The 24 labellings: (a,b,c,d) takes p1 = P(a,b), p2 = P(b,c), p3 = P(c,d),
/// p4 = P(d,a), s12 = P(a,c), s23 = P(b,d) and masses (m_a, m_b, m_c, m_d).
fn labellings() -> Vec<[usize; 4]> {
    let mut result = Vec::new();
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
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

fn relabel<V: Clone>(invariants: &[V; 6], masses: &[V; 4], [a, b, c, d]: [usize; 4]) -> [V; 10] {
    let slot = |i: usize, k: usize| match (i.min(k), i.max(k)) {
        (0, 1) => 0,
        (1, 2) => 1,
        (2, 3) => 2,
        (0, 3) => 3,
        (0, 2) => 4,
        (1, 3) => 5,
        _ => unreachable!(),
    };
    let pair = |i: usize, k: usize| invariants[slot(i, k)].clone();
    [
        pair(a, b),
        pair(b, c),
        pair(c, d),
        pair(d, a),
        pair(a, c),
        pair(b, d),
        masses[a].clone(),
        masses[b].clone(),
        masses[c].clone(),
        masses[d].clone(),
    ]
}

fn reference(case: &Case) -> [[Float; 2]; 3] {
    case.reference
        .map(|parts| parts.map(|s| Float::parse(s, Some(BITS)).unwrap()))
}

/// Relative error of every coefficient. A zero reference coefficient (no
/// double pole) is compared absolutely, against the largest coefficient.
fn check_values(
    values: &[[Float; 2]; 3],
    expected: &[[Float; 2]; 3],
    tolerance: f64,
    label: &str,
    failures: &mut Vec<String>,
) {
    let norm = |z: &[Float; 2]| (z[0].clone() * &z[0] + z[1].clone() * &z[1]).sqrt();
    let scale = expected
        .iter()
        .map(norm)
        .fold(Float::new(BITS), |a, b| if b > a { b } else { a });
    for (coefficient, (value, expected)) in values.iter().zip(expected).enumerate() {
        let size = norm(expected);
        let size = if size.is_zero() { scale.clone() } else { size };
        let error = norm(&[
            value[0].clone() - &expected[0],
            value[1].clone() - &expected[1],
        ]);
        if !error.is_finite() || error > size.clone() * Float::with_val(BITS, tolerance) {
            failures.push(format!(
                "{label}, coefficient {coefficient}: {:e} {:e} (relative error {:e})",
                value[0],
                value[1],
                error / size
            ));
        }
    }
}

fn binary64_input(case: &Case, order: [usize; 4]) -> Vec<Complex<f64>> {
    relabel(&case.invariants, &case.masses, order)
        .into_iter()
        .chain([MU_SQUARED])
        .map(|v| Complex::new(v, 0.0))
        .collect()
}

fn widen(values: &[Complex<f64>]) -> [[Float; 2]; 3] {
    core::array::from_fn(|k| {
        [
            Float::with_val(BITS, values[k].re),
            Float::with_val(BITS, values[k].im),
        ]
    })
}

/// The bare public master `oneloopmaster::D0(tag, ...)`, without a function
/// map: Symbolica evaluates it through the master's numeric callbacks.
fn master_callbacks<T: NativeFloat>(prototype: &T, bits: u32) -> ExpressionEvaluator<Complex<T>>
where
    Complex<T>: EvaluationDomain,
{
    let parameters = (0..11)
        .map(|i| symbolica::symbol!(format!("d0_infrared_relabelling::bare{i}")).to_atom())
        .collect::<Vec<_>>();
    let master = oneloop::D0();
    let expressions = [0, -1, -2].map(|tag| {
        let mut arguments = vec![Atom::num(tag)];
        arguments.extend_from_slice(&parameters);
        master.call(&arguments)
    });
    Atom::evaluator_multiple(&expressions, &parameters)
        .direct_translation(true)
        .build()
        .unwrap()
        .map_coeff_with_prec(
            &|value| {
                Complex::new(
                    prototype.from_rational(&value.re),
                    prototype.from_rational(&value.im),
                )
            },
            bits,
        )
}

#[test]
fn binary64_routes_agree_in_every_labelling() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let mut handwritten = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
            let mut generated = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
            generated.set_handwritten(false);
            let mut callbacks = master_callbacks(&0.0_f64, 53);
            let mut failures = Vec::new();
            for case in &CASES {
                let expected = reference(case);
                for order in labellings() {
                    let input = binary64_input(case, order);
                    for route in [
                        "handwritten",
                        "generated",
                        "expression",
                        "symjit",
                        "callbacks",
                    ] {
                        let mut output = [Complex::new(f64::NAN, f64::NAN); 3];
                        match route {
                            "handwritten" => handwritten.evaluate(&input, &mut output).unwrap(),
                            "generated" => generated.evaluate(&input, &mut output).unwrap(),
                            "expression" => evaluate_with_backend(
                                ScalarIntegral::D0,
                                &input,
                                &mut output,
                                EvaluationBackend::Expression,
                            )
                            .unwrap(),
                            "callbacks" => callbacks.evaluate(&input, &mut output),
                            _ => evaluate_with_backend(
                                ScalarIntegral::D0,
                                &input,
                                &mut output,
                                EvaluationBackend::SymJit,
                            )
                            .unwrap(),
                        }
                        check_values(
                            &widen(&output),
                            &expected,
                            1e-12,
                            &format!("{}, {route}, labelling {order:?}", case.name),
                            &mut failures,
                        );
                    }
                }
            }
            assert!(failures.is_empty(), "{}", failures.join("\n"));
        })
        .unwrap()
        .join()
        .unwrap();
}

/// Generated code, the transparent expressions (relabelled by the caller) and
/// the master's numeric callbacks in carrier `T`.
fn check_carrier<T: NativeFloat>(
    bits: u32,
    from_f64: impl Fn(f64) -> T,
    convert: impl Fn(&T) -> Float,
    tolerance: f64,
    failures: &mut Vec<String>,
) where
    Complex<T>: EvaluationDomain,
{
    let context = OneLoopExpressions::new();
    let variables = (0..11)
        .map(|i| symbolica::symbol!(format!("d0_infrared_relabelling::x{i}")).to_atom())
        .collect::<Vec<_>>();
    let series = context.d0(
        core::array::from_fn(|i| &variables[i]),
        core::array::from_fn(|i| &variables[i + 6]),
        &variables[10],
    );
    let prototype = T::zero_at(bits);
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
    let mut callbacks = master_callbacks(&prototype, bits);
    let mut evaluate = |input: &[Complex<T>], route: &str| {
        let mut output =
            core::array::from_fn::<_, 3, _>(|_| Complex::new(prototype.zero(), prototype.zero()));
        match route {
            "generated" => native.evaluate(input, &mut output).unwrap(),
            "callbacks" => callbacks.evaluate(input, &mut output),
            _ => {
                let mut point = input.to_vec();
                oneloop::regular_d0_labelling(&mut point);
                expression.evaluate(&point, &mut output);
            }
        }
        output.map(|z| [convert(&z.re), convert(&z.im)])
    };
    for case in &CASES {
        let expected = reference(case);
        for order in labellings() {
            let input = binary64_input(case, order)
                .into_iter()
                .map(|z| Complex::new(from_f64(z.re), from_f64(z.im)))
                .collect::<Vec<_>>();
            for route in ["expression", "generated", "callbacks"] {
                check_values(
                    &evaluate(&input, route),
                    &expected,
                    tolerance,
                    &format!("{}, bits {bits}, {route}, labelling {order:?}", case.name),
                    failures,
                );
            }
        }
    }
    // pyAmpliCol's exact rational tth inputs (nearest carrier values): no
    // labelling may return NaN or lose the pole; they differ from the
    // binary64-rounded inputs at the 2e-16 relative level.
    let exact = TTH_EXACT.map(|(n, d)| prototype.from_rational(&Rational::from((n, d))));
    let masses = CASES[0].masses.map(&from_f64);
    for order in labellings() {
        let input = relabel(&exact, &masses, order)
            .into_iter()
            .chain([from_f64(MU_SQUARED)])
            .map(|v| Complex::new(v, prototype.zero()))
            .collect::<Vec<_>>();
        for route in ["expression", "generated", "callbacks"] {
            check_values(
                &evaluate(&input, route),
                &reference(&CASES[0]),
                1e-15,
                &format!("tth exact, bits {bits}, {route}, labelling {order:?}"),
                failures,
            );
        }
    }
}

#[test]
fn extended_carriers_agree_in_every_labelling() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let mut failures = Vec::new();
            check_carrier(
                106,
                DoubleFloat::from,
                |x| {
                    let pair = x.into_inner();
                    Float::with_val(BITS, pair.hi()) + Float::with_val(BITS, pair.lo())
                },
                1e-26,
                &mut failures,
            );
            for bits in [128, 256] {
                check_carrier(
                    bits,
                    |v| Float::with_val(bits, v),
                    Clone::clone,
                    1e-29,
                    &mut failures,
                );
            }
            assert!(failures.is_empty(), "{}", failures.join("\n"));
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn interpreted_float_agrees_in_every_labelling() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let bits = 128;
            let mut evaluator = PrecisionEvaluator::with_binary_precision_and_backend(
                ScalarIntegral::D0,
                bits,
                EvaluationBackend::Expression,
            )
            .unwrap();
            let mut failures = Vec::new();
            for case in &CASES {
                let expected = reference(case);
                for order in labellings() {
                    let input = binary64_input(case, order)
                        .into_iter()
                        .map(|z| {
                            Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im))
                        })
                        .collect::<Vec<_>>();
                    let mut output = core::array::from_fn::<_, 3, _>(|_| {
                        Complex::new(Float::new(bits), Float::new(bits))
                    });
                    evaluator.evaluate(&input, &mut output).unwrap();
                    check_values(
                        &output.map(|z| [z.re, z.im]),
                        &expected,
                        1e-29,
                        &format!(
                            "{}, interpreted Float {bits}, labelling {order:?}",
                            case.name
                        ),
                        &mut failures,
                    );
                }
            }
            assert!(failures.is_empty(), "{}", failures.join("\n"));
        })
        .unwrap()
        .join()
        .unwrap();
}
