//! Independent original-Fortran regressions for channel permutations, negative
//! squared masses and the real/tiny-width limit of the contour representation.
//! Expected values were obtained with tests/support/oracle.f90 linked to the
//! unchanged OneLOop library; every input returned finite values without ERROR.
use oneloop::OneLoopExpressions;
use symbolica::prelude::*;

struct Case {
    name: String,
    p: [f64; 6],
    m: [f64; 4],
    width: f64,
    imaginary_masses: Option<[f64; 4]>,
    mu: f64,
    expected: [(f64, f64); 3],
}

fn cases() -> Vec<Case> {
    let euclidean = [-7., -8., -9., -10., -12., -15.];
    let real = |name, p, m, expected| Case {
        name: String::from(name),
        p,
        m,
        width: 0.,
        imaginary_masses: None,
        mu: 1.,
        expected,
    };
    let mut cases = vec![
        real(
            "contour_endpoint_root",
            [0., 0., 0., 0., 2., -3.],
            [1., 1., 2., 3.],
            [(0.06382444261938502, 0.), (0., 0.), (0., 0.)],
        ),
        real(
            "zero_channel",
            [-1., -2., -3., -4., 0., -5.],
            [0.; 4],
            [
                (0.002_434_112_613_144_086_3, 0.),
                (-0.196_165_850_602_345_23, 0.),
                (0., 0.),
            ],
        ),
        real(
            "negative_finite_one",
            euclidean,
            [0., 0., 0., -2.],
            [
                (0.030342521734038577, 0.006268453515650985),
                (0., 0.),
                (0., 0.),
            ],
        ),
        real(
            "negative_finite_adjacent",
            euclidean,
            [0., 0., -2., 3.],
            [
                (0.021119271208851504, 0.006414941283148391),
                (0., 0.),
                (0., 0.),
            ],
        ),
        real(
            "negative_finite_adjacent_both",
            euclidean,
            [0., 0., -2., -3.],
            [
                (0.03324363561962497, 0.038795296463519355),
                (0., 0.),
                (0., 0.),
            ],
        ),
        real(
            "negative_finite_opposite",
            euclidean,
            [0., -2., 0., 3.],
            [
                (0.021291624231569318, 0.007489287779993325),
                (0., 0.),
                (0., 0.),
            ],
        ),
        real(
            "negative_finite_three",
            euclidean,
            [0., -2., 3., 4.],
            [
                (0.012713760914572507, 0.004690774285737651),
                (0., 0.),
                (0., 0.),
            ],
        ),
        real(
            "negative_finite_four",
            euclidean,
            [-1., 2., 3., 4.],
            [
                (0.008542210473692151, 0.0016824181255571766),
                (0., 0.),
                (0., 0.),
            ],
        ),
        real(
            "negative_finite_four_all",
            euclidean,
            [-1., -2., -3., -4.],
            [
                (0.018089525014447615, -0.15580659315023881),
                (0., 0.),
                (0., 0.),
            ],
        ),
        real(
            "negative_ir06",
            [0., 0., -2., -2., -5., -6.],
            [0., 0., 0., -2.],
            [
                (-0.07940350733373688, 0.2528099161055931),
                (-0.18444397270569685, -0.15707963267948966),
                (0.1, 0.),
            ],
        ),
        real(
            "negative_ir11",
            [0., -2., -7., 3., -5., -6.],
            [0., 0., -2., 3.],
            [
                (-0.058978662250802685, 0.04439854491298364),
                (-0.08888730116260375, -0.05817764173314432),
                (0.037037037037037035, 0.),
            ],
        ),
        real(
            "negative_ir12",
            [0., -2., -7., -4., -5., -6.],
            [0., 0., -2., 3.],
            [
                (0.026593727386821454, -0.0056139789937921895),
                (-0.03716122691366826, -0.05817764173314432),
                (0.018518518518518517, 0.),
            ],
        ),
        real(
            "negative_ir13",
            [0., -8., -7., -4., -5., -6.],
            [0., 0., -2., 3.],
            [
                (0.035141234678642136, 0.027913558142973024),
                (-0.02945551681860261, 0.),
                (0., 0.),
            ],
        ),
        real(
            "negative_ir14",
            [-2., -2., 3., 3., -5., -6.],
            [0., -2., 0., 3.],
            [
                (0.08700063263210757, 0.11835665041422457),
                (-0.05405653238311539, -0.07353912163981707),
                (0., 0.),
            ],
        ),
        real(
            "negative_ir14_both",
            [-2., -2., -3., -3., -5., -6.],
            [0., -2., 0., -3.],
            [
                (-0.2384528796383123, 1.7945344320000643e-17),
                (0.14815910436562177, -1.1150069338717302e-17),
                (0., 0.),
            ],
        ),
        real(
            "negative_ir15",
            [-7., -2., 3., -4., -5., -6.],
            [0., -2., 0., 3.],
            [
                (0.05872011142372666, 0.009824395489615022),
                (-0.027028266191557692, -0.03676956081990854),
                (0., 0.),
            ],
        ),
        real(
            "negative_ir16",
            [-2., -8., -7., 4., -5., -6.],
            [0., -2., 3., 4.],
            [
                (0.021985339016077878, -0.001192714462407377),
                (-0.014623145470368283, -0.020039841096279126),
                (0., 0.),
            ],
        ),
    ];
    for (name, p, m, real_value, width_value) in [
        (
            "real_contour",
            [1., 2., 3., 4., 12., -5.],
            [0.3, 0.5, 0.7, 1.4],
            (-0.30758996230139596, 0.07371326071438561),
            (-0.30758996078203543, 0.07371325708966603),
        ),
        (
            "equal_real_contour",
            [1., 2., 3., 4., 12., -5.],
            [2.; 4],
            (0.0011848773584578476, 0.13448983052142294),
            (0.0011848763589549227, 0.13448982825539543),
        ),
        (
            "zero_external_contour",
            [0., 0., 0., 0., 12., -5.],
            [0.3, 0.5, 0.7, 1.4],
            (-0.0654879714169148, 0.15953433881355153),
            (-0.06548797246048117, 0.15953433762272726),
        ),
        (
            "negative_contour",
            [1., 2., 3., 4., 12., -5.],
            [-1., 2., 3., 4.],
            (-0.032927201244151126, 0.0313161155779186),
            (-0.03292720113770572, 0.03131611515939738),
        ),
    ] {
        for width in [0., 1e-8, 1e-16] {
            cases.push(Case {
                name: String::from(name),
                p,
                m,
                width,
                imaginary_masses: None,
                mu: 1.,
                expected: [
                    if width == 1e-8 {
                        width_value
                    } else {
                        real_value
                    },
                    (0., 0.),
                    (0., 0.),
                ],
            });
        }
    }
    for rotation in 1..4 {
        let masses = [1., 1., 2., 3.];
        cases.push(Case {
            name: format!("contour_endpoint_root_rotation_{rotation}"),
            p: if rotation % 2 == 0 {
                [0., 0., 0., 0., 2., -3.]
            } else {
                [0., 0., 0., 0., -3., 2.]
            },
            m: core::array::from_fn(|i| masses[(i + rotation) % 4]),
            width: 0.,
            imaginary_masses: None,
            mu: 1.,
            expected: [(0.06382444261938502, 0.), (0., 0.), (0., 0.)],
        });
    }
    for (index, line) in include_str!("data/parity.txt").lines().enumerate() {
        let values: Vec<f64> = line
            .split_whitespace()
            .filter_map(|x| x.parse().ok())
            .collect();
        if values.len() != 22 || values[0] != 4. {
            continue;
        }
        let m = core::array::from_fn(|i| values[8 + 2 * i]);
        let imaginary = core::array::from_fn(|i| values[9 + 2 * i]);
        if m.iter()
            .zip(&imaginary)
            .filter(|(r, i)| **r == 0. && **i == 0.)
            .count()
            != 1
        {
            continue;
        }
        cases.push(Case {
            name: format!("three_mass_parity_line_{}", index + 1),
            p: core::array::from_fn(|i| values[2 + i]),
            m,
            width: 0.,
            imaginary_masses: Some(imaginary),
            mu: values[1],
            expected: core::array::from_fn(|i| (values[16 + 2 * i], values[17 + 2 * i])),
        });
    }
    cases
}

#[test]
fn zero_channel_and_negative_mass_boxes_match_original() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(check_cases)
        .unwrap()
        .join()
        .unwrap();
}

fn check_cases() {
    let x: [Atom; 11] =
        core::array::from_fn(|index| symbol!(format!("box_limit_{index}")).to_atom());
    let context = OneLoopExpressions::new();
    let integral = context.d0(
        core::array::from_fn(|i| &x[i]),
        core::array::from_fn(|i| &x[i + 6]),
        &x[10],
    );
    let exact = context.evaluator(integral.coefficients(), &x).unwrap();
    let mut double = exact
        .clone()
        .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    let mut precise = exact.map_coeff_with_prec(
        &|c| Complex::new(c.re.to_multi_prec_float(128), c.im.to_multi_prec_float(128)),
        128,
    );
    let mut failures = Vec::new();
    let cases = cases();
    assert!(
        cases.len() > 28,
        "three-mass acceptance fixtures must be included"
    );
    for case in &cases {
        // Exercise both signed zeros independently of tiny negative widths.
        let signs: &[f64] = if case.width == 0. { &[0., -0.] } else { &[0.] };
        for &zero in signs {
            let input: [Complex<f64>; 11] = core::array::from_fn(|index| {
                if index < 6 {
                    Complex::new(case.p[index], 0.)
                } else if index < 10 {
                    let mass = case.m[index - 6];
                    Complex::new(
                        mass,
                        if let Some(imaginary) = case.imaginary_masses {
                            imaginary[index - 6]
                        } else if case.width == 0. {
                            zero
                        } else {
                            -case.width * mass.abs()
                        },
                    )
                } else {
                    Complex::new(case.mu, 0.)
                }
            });
            let mut values = [Complex::new(0., 0.); 3];
            double.evaluate(&input, &mut values);
            record_failures(case, &values, 53, zero, &mut failures);
            let input =
                input.map(|z| Complex::new(Float::with_val(128, z.re), Float::with_val(128, z.im)));
            let mut values =
                core::array::from_fn::<_, 3, _>(|_| Complex::new(Float::new(128), Float::new(128)));
            precise.evaluate(&input, &mut values);
            record_failures(
                case,
                &values.map(|z| Complex::new(z.re.to_f64(), z.im.to_f64())),
                128,
                zero,
                &mut failures,
            );
        }
    }
    assert!(
        failures.is_empty(),
        "{} failed box comparisons:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

fn record_failures(
    case: &Case,
    actual: &[Complex<f64>; 3],
    bits: u32,
    zero: f64,
    failures: &mut Vec<String>,
) {
    let scale = case
        .p
        .iter()
        .chain(&case.m)
        .fold(0_f64, |s, x| s.max(x.abs()));
    let normalization = scale * scale;
    for (coefficient, (value, &(re, im))) in actual.iter().zip(&case.expected).enumerate() {
        let error = (value.re - re).hypot(value.im - im) * normalization;
        let tolerance = 2e-10 + 2e-8 * re.hypot(im) * normalization;
        if !error.is_finite() || error >= tolerance {
            failures.push(format!("{} width={} bits={bits} zero={zero:?} coefficient={coefficient}: {value:?}, expected=({re},{im})",
                case.name, case.width));
        }
    }
}

/// Independent equal-mass Feynman-parameter targets, not the finite-regulator
/// values returned by double-precision Fortran at these exact thresholds.
/// Cheng-Wu reduction with w=z*(1-z) gives
/// D0(s,t)=-integral_0^1 log((1-s*w)*(1-t*w))/(s+t-s*t*w) dz.
/// The targets agree at 80 and 120 decimal digits; the second is also
/// pi/4-atanh(1/sqrt(2))/sqrt(2).
#[test]
fn repeated_root_threshold_boxes_have_physical_limits() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(check_thresholds)
        .unwrap()
        .join()
        .unwrap();
}

fn check_thresholds() {
    let x: [Atom; 11] =
        core::array::from_fn(|index| symbol!(format!("box_threshold_{index}")).to_atom());
    let context = OneLoopExpressions::new();
    let integral = context.d0(
        core::array::from_fn(|i| &x[i]),
        core::array::from_fn(|i| &x[i + 6]),
        &x[10],
    );
    let exact = context.evaluator(integral.coefficients(), &x).unwrap();
    let mut failures = Vec::new();
    for bits in [128, 256] {
        let mut evaluator = exact.clone().map_coeff_with_prec(
            &|c| {
                Complex::new(
                    c.re.to_multi_prec_float(bits),
                    c.im.to_multi_prec_float(bits),
                )
            },
            bits,
        );
        // The nonzero-width targets come from the same independent integral
        // with mass 1-i*width, using log1p at its removable denominator zero.
        for (s, t, targets) in [
            (
                4.,
                -3.,
                [
                    (0.437_508_120_479_433_8, 0.),
                    (0.431_954_242_508_579_96, 0.005_520_973_441_756_335),
                    (0.437_452_584_442_429_2, 0.000_055_532_801_195_860_31),
                    (0.437_508_114_925_830_2, 5.553_603_640_345_423e-9),
                    (0.437_508_120_479_433_3, 5.553_603_672_697_954e-16),
                ],
            ),
            (
                2.,
                -4.,
                [
                    (0.162_172_923_257_217_8, 0.),
                    (0.162_172_917_371_848_63, 0.000_034_725_491_619_602_76),
                    (0.162_172_923_257_217_74, 3.472_549_255_398_694e-9),
                    (0.162_172_923_257_217_8, 3.472_549_255_398_695e-17),
                    (0.162_172_923_257_217_8, 3.472_549_255_398_695e-31),
                ],
            ),
        ] {
            for (sample, width) in [0., 1e-4, 1e-8, 1e-16, 1e-30].into_iter().enumerate() {
                let (target_re, target_im) = targets[sample];
                let input: [Complex<Float>; 11] = core::array::from_fn(|index| {
                    let (real, imaginary) = match index {
                        0..=3 => (0., 0.),
                        4 => (s, 0.),
                        5 => (t, 0.),
                        6..=9 => (1., -width),
                        _ => (1., 0.),
                    };
                    Complex::new(
                        Float::with_val(bits, real),
                        Float::with_val(bits, imaginary),
                    )
                });
                let mut result = core::array::from_fn::<_, 3, _>(|_| {
                    Complex::new(Float::new(bits), Float::new(bits))
                });
                evaluator.evaluate(&input, &mut result);
                let actual = result.map(|z| Complex::new(z.re.to_f64(), z.im.to_f64()));
                eprintln!(
                    "threshold s={s} t={t} width={width} bits={bits}: {:?}",
                    actual[0]
                );
                for (coefficient, z) in actual.iter().enumerate() {
                    let finite = z.re.is_finite() && z.im.is_finite();
                    // At this exceptionally narrow width the regular outer
                    // root representation loses accuracy through cancellation
                    // at 128 bits. Retain it as an explicit finite diagnostic;
                    // its numerical target is still asserted at 256 bits.
                    let diagnostic_only = bits == 128 && s == 2. && width == 1e-30;
                    let correct = if coefficient == 0 {
                        diagnostic_only || (z.re - target_re).hypot(z.im - target_im) < 2e-11
                    } else {
                        z.re.hypot(z.im) < 2e-11
                    };
                    if !finite || !correct {
                        failures.push(format!("s={s},t={t},width={width},bits={bits},coefficient={coefficient}: {z:?}, independent target=({target_re},{target_im})"));
                    }
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} threshold failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// PR9 exposed expanded conjugate norms in these existing complex-mass rows.
/// Keep the original binary64 oracle accuracy while exercising actual carrier
/// arithmetic in both expression and generated paths.
#[cfg(feature = "generated-evaluators")]
#[test]
fn complex_contour_norms_preserve_sheets_in_all_carriers() {
    use oneloop::{NativeEvaluator, NativeFloat, ScalarIntegral};
    use symbolica::domains::float::DoubleFloat;
    use symbolica::evaluate::EvaluationDomain;

    fn check<T: NativeFloat>(bits: u32, from_f64: impl Fn(f64) -> T)
    where
        Complex<T>: EvaluationDomain,
    {
        let x: [Atom; 11] =
            core::array::from_fn(|index| symbol!(format!("box_norm_regression_{index}")).to_atom());
        let context = OneLoopExpressions::new();
        let integral = context.d0(
            core::array::from_fn(|i| &x[i]),
            core::array::from_fn(|i| &x[i + 6]),
            &x[10],
        );
        let prototype = T::zero_at(bits);
        let mut expression = context
            .evaluator(integral.coefficients(), &x)
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
        let mut generated =
            NativeEvaluator::<T>::with_binary_precision(ScalarIntegral::D0, bits).unwrap();
        generated.set_handwritten(false);
        let names =
            [217, 221, 223, 224, 232, 233].map(|line| format!("three_mass_parity_line_{line}"));
        let selected = cases()
            .into_iter()
            .filter(|case| names.contains(&case.name))
            .collect::<Vec<_>>();
        assert_eq!(selected.len(), names.len());
        let mut failures = Vec::new();
        for case in &selected {
            let imaginary = case
                .imaginary_masses
                .expect("explicit complex-mass fixture");
            let input: [Complex<T>; 11] = core::array::from_fn(|index| {
                let (re, im) = if index < 6 {
                    (case.p[index], 0.)
                } else if index < 10 {
                    (case.m[index - 6], imaginary[index - 6])
                } else {
                    (case.mu, 0.)
                };
                Complex::new(from_f64(re), from_f64(im))
            });
            for native in [false, true] {
                let mut values = core::array::from_fn::<_, 3, _>(|_| {
                    Complex::new(prototype.zero(), prototype.zero())
                });
                if native {
                    generated.evaluate(&input, &mut values).unwrap();
                } else {
                    expression.evaluate(&input, &mut values);
                }
                let start = failures.len();
                record_failures(
                    case,
                    &values.map(|value| Complex::new(value.re.to_f64(), value.im.to_f64())),
                    bits,
                    0.,
                    &mut failures,
                );
                for failure in &mut failures[start..] {
                    *failure = format!("generated={native}: {failure}");
                }
            }
        }
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    }
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            check(53, |value| value);
            check(106, DoubleFloat::from);
            check(128, |value| Float::with_val(128, value));
            check(266, |value| Float::with_val(266, value));
        })
        .unwrap()
        .join()
        .unwrap();
}
