use super::*;
use symbolica::{
    domains::{float::Complex, rational::Rational},
    parse,
};

#[test]
fn removable_log_limit_is_finite() {
    let one = Atom::num(1);
    assert_eq!(log_over_one_minus(&one), Atom::num(-1));
}

fn evaluate_series<const N: usize>(
    expressions: &[Atom; 3],
    parameters: &[Atom],
    input: &[Complex<f64>; N],
) -> [Complex<f64>; 3] {
    let views = expressions.iter().map(Atom::as_view).collect::<Vec<_>>();
    let mut map = FunctionMap::new();
    sheet_exact::register(&mut map);
    let evaluator = Atom::evaluator_multiple(&views, parameters)
        .function_map(map.into())
        .build()
        .unwrap();
    let mut evaluator = evaluator
        .map_coeff(&|value: &Complex<Rational>| Complex::new(value.re.to_f64(), value.im.to_f64()));
    let mut output = [Complex::new(0.0, 0.0); 3];
    evaluator.evaluate(input, &mut output);
    output
}

#[test]
fn a0_is_a_transparent_native_expression() {
    let m = parse!("m");
    let mu2 = parse!("mu2");
    let coefficients = a0(&m, &mu2).into_coefficients();

    let expected = if_nonzero(&m, &m * (1 - physical_log(&(&m / &mu2))));
    assert_eq!(coefficients[0], expected);
    assert_eq!(coefficients[1], m);
    assert_eq!(coefficients[2], parse!("0"));
}

#[test]
fn a0_zero_mass_is_scaleless() {
    let coefficients = a0(&parse!("0"), &parse!("mu2")).into_coefficients();
    assert_eq!(coefficients, [parse!("0"), parse!("0"), parse!("0")]);
}

#[test]
fn bubble_axis_lips_survive_arbitrary_precision_arithmetic() {
    use symbolica::prelude::{Float, RealLike};
    let p = parse!("axis_p");
    let m0 = parse!("axis_m0");
    let m1 = parse!("axis_m1");
    let linear = &m1 - &m0 - &p;
    let delta = (&linear * &linear - Atom::num(4) * &p * &m0).sqrt();
    let root = (-linear - &delta) / (Atom::num(2) * &p);
    let ratio = Atom::num(1) - Atom::num(1) / &root;
    let off_axis = off_real_axis(&ratio);
    let expressions = [
        delta,
        root,
        ratio.clone(),
        off_axis,
        upper_log(&ratio),
        physical_log(&ratio),
    ];
    let views = expressions.each_ref().map(Atom::as_view);
    let e = Atom::evaluator_multiple(&views, &[p, m0, m1])
        .build()
        .unwrap();
    for bits in [256, 384] {
        let mut evaluator = e.clone().map_coeff_with_prec(
            &|c| {
                Complex::new(
                    c.re.to_multi_prec_float(bits),
                    c.im.to_multi_prec_float(bits),
                )
            },
            bits,
        );
        let input = [7, 1, 2].map(|v| Complex::new(Float::with_val(bits, v), Float::new(bits)));
        let mut out =
            core::array::from_fn::<_, 6, _>(|_| Complex::new(Float::new(bits), Float::new(bits)));
        evaluator.evaluate(&input, &mut out);
        eprintln!("branch diagnostic {bits}bits: {out:?}");
        assert_eq!(out[3].re.to_f64(), 0.);
        assert_eq!(out[3].im.to_f64(), 0.);
        assert!((out[4].im.to_f64() - std::f64::consts::PI).abs() < 1e-12);
        assert!((out[5].im.to_f64() + std::f64::consts::PI).abs() < 1e-12);
    }
}

#[test]
fn tensor_tadpoles_use_exact_rationals() {
    let coefficients = an(4, &parse!("m"), &parse!("mu2")).unwrap();
    assert_eq!(coefficients.len(), 3);
    let m = parse!("m");
    let logarithm = physical_log(&parse!("m/mu2"));
    assert_eq!(
        coefficients[1].coefficients()[0],
        if_nonzero(&m, &m * &m * (parse!("3/2") - &logarithm) / 4),
    );
    assert_eq!(
        coefficients[2].coefficients()[0],
        if_nonzero(&m, &m * &m * &m * (parse!("11/6") - logarithm) / 24),
    );
    assert_eq!(coefficients[1].coefficients()[1], parse!("m^2/4"));
    assert_eq!(coefficients[2].coefficients()[1], parse!("m^3/24"));
}

#[test]
fn tensor_rank_contract_matches_oneloop() {
    assert_eq!(an(1, &parse!("m"), &parse!("mu2")).unwrap().len(), 1);
    assert_eq!(an(3, &parse!("m"), &parse!("mu2")).unwrap().len(), 2);
    assert_eq!(
        an(5, &parse!("m"), &parse!("mu2")).unwrap_err(),
        UnsupportedRank(5)
    );
}

#[test]
fn a0_matches_the_existing_rust_reference() {
    let mass = parse!("m");
    let mu_squared = parse!("mu2");
    let expressions = a0(&mass, &mu_squared).into_coefficients();
    let views = expressions.iter().map(Atom::as_view).collect::<Vec<_>>();
    let parameters = [mass, mu_squared];
    let evaluator = Atom::evaluator_multiple(&views, &parameters)
        .build()
        .unwrap();
    let mut evaluator = evaluator
        .map_coeff(&|value: &Complex<Rational>| Complex::new(value.re.to_f64(), value.im.to_f64()));
    let mut output = [Complex::new(0.0, 0.0); 3];
    evaluator.evaluate(
        &[Complex::new(100.0, -1.4), Complex::new(1.0, 0.0)],
        &mut output,
    );

    let expected = Complex::new(-360.507_218_918_917_4, 6.447_283_991_027_792);
    let finite_difference = output[0] - expected;
    assert!(finite_difference.re.hypot(finite_difference.im) < 2.0e-11);
    let pole_difference = output[1] - Complex::new(100.0, -1.4);
    assert!(pole_difference.re.hypot(pole_difference.im) < 2.0e-13);
    assert_eq!(output[2], Complex::new(0.0, 0.0));

    evaluator.evaluate(
        &[
            Complex::new(1.0, 0.0),
            Complex::new(91.1876_f64.powi(2), 0.0),
        ],
        &mut output,
    );
    assert!((output[0].re - 10.025_837_845_847_876).abs() < 2.0e-13);
    assert_eq!(output[0].im, 0.0);
}

#[test]
fn b0_is_a_transparent_native_expression() {
    let coefficients =
        b0(&parse!("p2"), &parse!("m0"), &parse!("m1"), &parse!("mu2")).into_coefficients();
    let printed = coefficients[0].to_string();

    assert!(printed.contains("if("));
    assert!(printed.contains("log("));
    assert!(printed.contains("conj("));
    assert!(!printed.contains("oneloop"));
    assert_eq!(coefficients[2], parse!("0"));
}

#[test]
fn b0_matches_existing_reference_and_elementary_limits() {
    let p2 = parse!("p2");
    let m0 = parse!("m0");
    let m1 = parse!("m1");
    let mu2 = parse!("mu2");
    let expressions = b0(&p2, &m0, &m1, &mu2).into_coefficients();
    let views = expressions.iter().map(Atom::as_view).collect::<Vec<_>>();
    let parameters = [p2, m0, m1, mu2];
    let evaluator = Atom::evaluator_multiple(&views, &parameters)
        .build()
        .unwrap();
    let mut evaluator = evaluator
        .map_coeff(&|value: &Complex<Rational>| Complex::new(value.re.to_f64(), value.im.to_f64()));
    let mut output = [Complex::new(0.0, 0.0); 3];

    evaluator.evaluate(
        &[
            Complex::new(3.0, 0.0),
            Complex::new(0.7, -0.03),
            Complex::new(1.4, -0.08),
            Complex::new(1.0, 0.0),
        ],
        &mut output,
    );
    let expected = Complex::new(0.720_288_906_550_704_1, 0.118_845_872_934_281_17);
    let difference = output[0] - expected;
    assert!(difference.re.hypot(difference.im) < 2.0e-11);
    assert_eq!(output[1], Complex::new(1.0, 0.0));

    evaluator.evaluate(
        &[
            Complex::new(-1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
        ],
        &mut output,
    );
    assert!((output[0].re - 2.0).abs() < 2.0e-13);
    assert!(output[0].im.abs() < 2.0e-13);

    evaluator.evaluate(
        &[
            Complex::new(0.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(1.0, 0.0),
        ],
        &mut output,
    );
    assert!((output[0].re + 2.0_f64.ln()).abs() < 2.0e-13);
    assert!(output[0].im.abs() < 2.0e-13);

    let additional_cases = [
        (
            [
                Complex::new(3.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.901_387_711_331_890_2, core::f64::consts::PI),
        ),
        (
            [
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(5.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-1.220_298_400_350_203_8, 0.0),
        ),
        (
            [
                Complex::new(-1.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.090_457_495_115_561_67, 0.0),
        ),
        (
            [
                Complex::new(2.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(1.306_852_819_440_054_6, 0.0),
        ),
    ];
    for (input, expected) in additional_cases {
        evaluator.evaluate(&input, &mut output);
        let difference = output[0] - expected;
        assert!(difference.re.hypot(difference.im) < 2.0e-11);
    }
}

#[test]
fn db0_is_a_transparent_native_expression() {
    let coefficients =
        db0(&parse!("p2"), &parse!("m0"), &parse!("m1"), &parse!("mu2")).into_coefficients();
    let printed = coefficients[0].to_string();

    assert!(printed.contains("if("));
    assert!(printed.contains("log("));
    assert!(!printed.contains("oneloop::internal"));
    assert_eq!(coefficients[2], parse!("0"));
}

#[test]
fn db0_matches_existing_reference_branches() {
    let p2 = parse!("p2");
    let m0 = parse!("m0");
    let m1 = parse!("m1");
    let mu2 = parse!("mu2");
    let expressions = db0(&p2, &m0, &m1, &mu2).into_coefficients();
    let views = expressions.iter().map(Atom::as_view).collect::<Vec<_>>();
    let parameters = [p2, m0, m1, mu2];
    let evaluator = Atom::evaluator_multiple(&views, &parameters)
        .build()
        .unwrap();
    let mut evaluator = evaluator
        .map_coeff(&|value: &Complex<Rational>| Complex::new(value.re.to_f64(), value.im.to_f64()));
    let mut output = [Complex::new(0.0, 0.0); 3];
    let cases = [
        (
            [
                Complex::new(3.0, 0.0),
                Complex::new(0.7, -0.03),
                Complex::new(1.4, -0.08),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.430_820_441_739_602_66, 0.065_385_279_130_906_14),
            Complex::new(0.0, 0.0),
        ),
        (
            [
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(5.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.049_521_951_157_720_36, 0.0),
            Complex::new(0.0, 0.0),
        ),
        (
            [
                Complex::new(-1.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.189_069_783_783_671_26, 0.0),
            Complex::new(0.0, 0.0),
        ),
        (
            [
                Complex::new(3.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-1.0 / 3.0, 0.0),
            Complex::new(0.0, 0.0),
        ),
        (
            [
                Complex::new(2.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.326_713_204_860_013_65, 0.0),
            Complex::new(-0.25, 0.0),
        ),
    ];

    for (input, expected_finite, expected_pole) in cases {
        evaluator.evaluate(&input, &mut output);
        let finite_difference = output[0] - expected_finite;
        let pole_difference = output[1] - expected_pole;
        assert!(finite_difference.re.hypot(finite_difference.im) < 2.0e-10);
        assert!(pole_difference.re.hypot(pole_difference.im) < 2.0e-13);
        assert_eq!(output[2], Complex::new(0.0, 0.0));
    }
}

#[test]
fn scaleless_db0_is_not_silently_zero() {
    let parameters = [
        parse!("scaleless_db0_p"),
        parse!("scaleless_db0_m0"),
        parse!("scaleless_db0_m1"),
        parse!("scaleless_db0_mu"),
    ];
    let generic = db0(
        &parameters[0],
        &parameters[1],
        &parameters[2],
        &parameters[3],
    );
    let output = evaluate_series(
        generic.coefficients(),
        &parameters,
        &[
            Complex::new(0., 0.),
            Complex::new(0., 0.),
            Complex::new(0., 0.),
            Complex::new(1., 0.),
        ],
    );
    assert!(!output[0].re.is_finite() || !output[0].im.is_finite());
    let zero = Atom::num(0);
    let direct = db0(&zero, &zero, &zero, &Atom::num(1));
    assert_eq!(
        direct.coefficients()[0],
        Atom::num(Coefficient::Indeterminate)
    );
}

#[test]
fn b1_is_a_transparent_native_expression() {
    let coefficients =
        b1(&parse!("p2"), &parse!("m0"), &parse!("m1"), &parse!("mu2")).into_coefficients();
    let printed = coefficients[0].to_string();

    assert!(printed.contains("if("));
    assert!(printed.contains("log("));
    assert!(!printed.contains("oneloop"));
    assert_eq!(coefficients[2], parse!("0"));
}

#[test]
fn b1_matches_existing_reference_branches() {
    let p2 = parse!("p2");
    let m0 = parse!("m0");
    let m1 = parse!("m1");
    let mu2 = parse!("mu2");
    let expressions = b1(&p2, &m0, &m1, &mu2).into_coefficients();
    let views = expressions.iter().map(Atom::as_view).collect::<Vec<_>>();
    let parameters = [p2, m0, m1, mu2];
    let evaluator = Atom::evaluator_multiple(&views, &parameters)
        .build()
        .unwrap();
    let mut evaluator = evaluator
        .map_coeff(&|value: &Complex<Rational>| Complex::new(value.re.to_f64(), value.im.to_f64()));
    let mut output = [Complex::new(0.0, 0.0); 3];
    let cases = [
        (
            [
                Complex::new(3.0, 0.0),
                Complex::new(0.7, -0.03),
                Complex::new(1.4, -0.08),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.271_938_401_813_376, -0.057_835_382_756_500_02),
        ),
        (
            [
                Complex::new(3.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.450_693_855_665_945_1, -core::f64::consts::FRAC_PI_2),
        ),
        (
            [
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(5.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.684_432_126_911_682_5, 0.0),
        ),
        (
            [
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.403_426_409_720_027_36, 0.0),
        ),
        (
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.096_573_590_279_972_64, 0.0),
        ),
        (
            [
                Complex::new(2.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.153_426_409_720_027_36, 0.0),
        ),
    ];

    for (input, expected_finite) in cases {
        evaluator.evaluate(&input, &mut output);
        let finite_difference = output[0] - expected_finite;
        assert!(finite_difference.re.hypot(finite_difference.im) < 2.0e-10);
        assert_eq!(output[1], Complex::new(-0.5, 0.0));
        assert_eq!(output[2], Complex::new(0.0, 0.0));
    }

    evaluator.evaluate(
        &[
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
        ],
        &mut output,
    );
    assert_eq!(output, [Complex::new(0.0, 0.0); 3]);
}

#[test]
fn b00_matches_existing_reference_branches() {
    let p2 = parse!("p2");
    let m0 = parse!("m0");
    let m1 = parse!("m1");
    let mu2 = parse!("mu2");
    let expressions = b00(&p2, &m0, &m1, &mu2).into_coefficients();
    let views = expressions.iter().map(Atom::as_view).collect::<Vec<_>>();
    let parameters = [p2, m0, m1, mu2];
    let evaluator = Atom::evaluator_multiple(&views, &parameters)
        .build()
        .unwrap();
    let mut evaluator = evaluator
        .map_coeff(&|value: &Complex<Rational>| Complex::new(value.re.to_f64(), value.im.to_f64()));
    let mut output = [Complex::new(0.0, 0.0); 3];
    let cases = [
        (
            [
                Complex::new(3.0, 0.0),
                Complex::new(0.7, -0.03),
                Complex::new(1.4, -0.08),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.404_032_001_204_310_1, -0.017_748_356_492_285_27),
            Complex::new(0.275, -0.0275),
        ),
        (
            [
                Complex::new(3.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.392_013_594_499_639_24, -core::f64::consts::FRAC_PI_4),
            Complex::new(-0.25, 0.0),
        ),
        (
            [
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(5.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.496_946_590_717_727_23, 0.0),
            Complex::new(1.75, 0.0),
        ),
        (
            [
                Complex::new(2.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.324_506_495_368_907_13, 0.0),
            Complex::new(1.0 / 3.0, 0.0),
        ),
    ];

    for (input, expected_finite, expected_pole) in cases {
        evaluator.evaluate(&input, &mut output);
        let finite_difference = output[0] - expected_finite;
        let pole_difference = output[1] - expected_pole;
        assert!(finite_difference.re.hypot(finite_difference.im) < 2.0e-10);
        assert!(pole_difference.re.hypot(pole_difference.im) < 2.0e-13);
        assert_eq!(output[2], Complex::new(0.0, 0.0));
    }
}

#[test]
fn b11_matches_existing_reference_branches() {
    let p2 = parse!("p2");
    let m0 = parse!("m0");
    let m1 = parse!("m1");
    let mu2 = parse!("mu2");
    let expressions = b11(&p2, &m0, &m1, &mu2).into_coefficients();
    let views = expressions.iter().map(Atom::as_view).collect::<Vec<_>>();
    let parameters = [p2, m0, m1, mu2];
    let evaluator = Atom::evaluator_multiple(&views, &parameters)
        .build()
        .unwrap();
    let mut evaluator = evaluator
        .map_coeff(&|value: &Complex<Rational>| Complex::new(value.re.to_f64(), value.im.to_f64()));
    let mut output = [Complex::new(0.0, 0.0); 3];
    let cases = [
        (
            [
                Complex::new(3.0, 0.0),
                Complex::new(0.7, -0.03),
                Complex::new(1.4, -0.08),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.124_287_648_889_346_21, 0.034_846_047_862_857_39),
        ),
        (
            [
                Complex::new(3.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.356_018_125_999_518_93, core::f64::consts::PI / 3.0),
        ),
        (
            [
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(5.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.478_829_006_058_197_2, 0.0),
        ),
        (
            [
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.380_062_050_924_462_66, 0.0),
        ),
        (
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.119_937_949_075_537_34, 0.0),
        ),
        (
            [
                Complex::new(2.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.008_826_837_964_426_194, 0.0),
        ),
    ];

    for (input, expected_finite) in cases {
        evaluator.evaluate(&input, &mut output);
        let finite_difference = output[0] - expected_finite;
        assert!(finite_difference.re.hypot(finite_difference.im) < 2.0e-10);
        assert!((output[1].re - 1.0 / 3.0).abs() < 2.0e-13);
        assert!(output[1].im.abs() < 2.0e-13);
        assert_eq!(output[2], Complex::new(0.0, 0.0));
    }
}

#[test]
fn b001_matches_existing_reference_branches() {
    let p2 = parse!("p2");
    let m0 = parse!("m0");
    let m1 = parse!("m1");
    let mu2 = parse!("mu2");
    let expressions = b001(&p2, &m0, &m1, &mu2).into_coefficients();
    let views = expressions.iter().map(Atom::as_view).collect::<Vec<_>>();
    let parameters = [p2, m0, m1, mu2];
    let evaluator = Atom::evaluator_multiple(&views, &parameters)
        .build()
        .unwrap();
    let mut evaluator = evaluator
        .map_coeff(&|value: &Complex<Rational>| Complex::new(value.re.to_f64(), value.im.to_f64()));
    let mut output = [Complex::new(0.0, 0.0); 3];
    let cases = [
        (
            [
                Complex::new(3.0, 0.0),
                Complex::new(0.7, -0.03),
                Complex::new(1.4, -0.08),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.212_244_325_353_773_4, 0.007_259_890_043_764_92),
            Complex::new(-1.0 / 6.0, 0.015_833_333_333_333_33),
        ),
        (
            [
                Complex::new(3.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.196_006_797_249_819_6, core::f64::consts::PI / 8.0),
            Complex::new(0.125, 0.0),
        ),
        (
            [
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(5.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.402_675_635_998_978_04, 0.0),
            Complex::new(-1.0, 0.0),
        ),
        (
            [
                Complex::new(2.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.201_713_204_860_013_68, 0.0),
            Complex::new(-0.25, 0.0),
        ),
    ];

    for (input, expected_finite, expected_pole) in cases {
        evaluator.evaluate(&input, &mut output);
        let finite_difference = output[0] - expected_finite;
        let pole_difference = output[1] - expected_pole;
        assert!(finite_difference.re.hypot(finite_difference.im) < 2.0e-10);
        assert!(pole_difference.re.hypot(pole_difference.im) < 2.0e-13);
        assert_eq!(output[2], Complex::new(0.0, 0.0));
    }
}

#[test]
fn b111_matches_existing_reference_branches() {
    let p2 = parse!("p2");
    let m0 = parse!("m0");
    let m1 = parse!("m1");
    let mu2 = parse!("mu2");
    let expressions = b111(&p2, &m0, &m1, &mu2).into_coefficients();
    let views = expressions.iter().map(Atom::as_view).collect::<Vec<_>>();
    let parameters = [p2, m0, m1, mu2];
    let evaluator = Atom::evaluator_multiple(&views, &parameters)
        .build()
        .unwrap();
    let mut evaluator = evaluator
        .map_coeff(&|value: &Complex<Rational>| Complex::new(value.re.to_f64(), value.im.to_f64()));
    let mut output = [Complex::new(0.0, 0.0); 3];
    let cases = [
        (
            [
                Complex::new(3.0, 0.0),
                Complex::new(0.7, -0.03),
                Complex::new(1.4, -0.08),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.061_060_888_224_836_24, -0.023_726_854_411_009_246),
        ),
        (
            [
                Complex::new(3.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(-0.308_680_261_166_305_8, -core::f64::consts::FRAC_PI_4),
        ),
        (
            [
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(5.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.368_684_627_151_776_97, 0.0),
        ),
        (
            [
                Complex::new(2.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            Complex::new(0.048_286_795_139_986_32, 0.0),
        ),
    ];

    for (input, expected_finite) in cases {
        evaluator.evaluate(&input, &mut output);
        let finite_difference = output[0] - expected_finite;
        assert!(finite_difference.re.hypot(finite_difference.im) < 2.0e-10);
        let pole_difference = output[1] - Complex::new(-0.25, 0.0);
        assert!(pole_difference.re.hypot(pole_difference.im) < 2.0e-13);
        assert_eq!(output[2], Complex::new(0.0, 0.0));
    }
}

#[test]
fn rank_four_bubbles_match_existing_reference() {
    let parameters = [parse!("p2"), parse!("m0"), parse!("m1"), parse!("mu2")];
    let complex_input = [
        Complex::new(3.0, 0.0),
        Complex::new(0.7, -0.03),
        Complex::new(1.4, -0.08),
        Complex::new(1.0, 0.0),
    ];
    let zero_input = [
        Complex::new(0.0, 0.0),
        Complex::new(2.0, 0.0),
        Complex::new(5.0, 0.0),
        Complex::new(1.0, 0.0),
    ];
    type BubbleTensor = fn(&Atom, &Atom, &Atom, &Atom) -> LaurentSeries;
    let functions: [BubbleTensor; 3] = [b0000, b0011, b1111];
    let expected_complex = [
        [
            Complex::new(0.084_177_235_220_498_43, -0.011_335_877_599_270_464),
            Complex::new(0.048_762_5, -0.008_291_666_666_666_666),
            Complex::new(0.0, 0.0),
        ],
        [
            Complex::new(0.148_472_977_849_870_76, -0.003_433_632_771_996_77),
            Complex::new(0.129_166_666_666_666_65, -0.011_25),
            Complex::new(0.0, 0.0),
        ],
        [
            Complex::new(0.029_939_862_301_517_683, 0.017_531_313_325_083_266),
            Complex::new(0.2, 0.0),
            Complex::new(0.0, 0.0),
        ],
    ];
    let expected_zero = [
        [
            Complex::new(0.262_019_977_641_903, 0.0),
            Complex::new(1.625, 0.0),
            Complex::new(0.0, 0.0),
        ],
        [
            Complex::new(-0.323_522_613_452_528_86, 0.0),
            Complex::new(0.708_333_333_333_333_3, 0.0),
            Complex::new(0.0, 0.0),
        ],
        // B1111 = ∫x⁴ [1/eps - log(2+3x)] dx. The former reference
        // +0.2998... carried OneLOop's sign error in its rank-four mass swap.
        [
            Complex::new(-0.299_847_502_997_086_6, 0.0),
            Complex::new(0.2, 0.0),
            Complex::new(0.0, 0.0),
        ],
    ];

    for ((function, expected_complex), expected_zero) in functions
        .into_iter()
        .zip(expected_complex)
        .zip(expected_zero)
    {
        let expressions = function(
            &parameters[0],
            &parameters[1],
            &parameters[2],
            &parameters[3],
        )
        .into_coefficients();
        for (input, expected) in [
            (&complex_input, expected_complex),
            (&zero_input, expected_zero),
        ] {
            let output = evaluate_series(&expressions, &parameters, input);
            for (actual, expected) in output.into_iter().zip(expected) {
                let difference = actual - expected;
                assert!(
                    difference.re.hypot(difference.im) < 2.0e-10,
                    "actual={actual:?}, expected={expected:?}"
                );
            }
        }
    }
}

#[test]
fn c0_matches_reference_sectors() {
    large_stack(c0_reference_sectors);
}

pub(super) fn large_stack(f: fn()) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(f)
        .unwrap()
        .join()
        .unwrap();
}

#[test]
#[ignore = "diagnoses construction cost of individual native definition groups"]
fn native_definition_construction_stages() {
    large_stack(|| {
        let mut map = FunctionMap::new();
        eprintln!("constructing sheet definitions");
        sheet_exact::register(&mut map);
        eprintln!("constructing triangle HV definitions");
        triangle_hv::register(&mut map);
        eprintln!("constructing complex box definitions");
        box_complex::register(&mut map);
        eprintln!("constructing triangle sector definitions");
        register_triangles(&mut map);
        eprintln!("constructing box sector definitions");
        register_boxes(&mut map);
        eprintln!("constructing public master definitions");
        register_masters(&mut map);
        eprintln!("all native definitions constructed");
    });
}

#[test]
fn triangle_single_width_regression() {
    large_stack(|| {
        let args = [
            "single_width_p1",
            "single_width_p2",
            "single_width_p3",
            "single_width_m",
        ]
        .map(|s| symbolica::symbol!(s).to_atom());
        let expression = triangle_finite_one_mass(&args[0], &args[1], &args[2], &args[3]);
        let mut map = FunctionMap::new();
        sheet_exact::register(&mut map);
        crate::triangle::register_helpers(&mut map);
        let mut evaluator = expression
            .evaluator(&args)
            .function_map(map.into())
            .direct_translation(true)
            .build()
            .unwrap()
            .map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64()));
        let out = evaluator.evaluate_single(&[
            Complex::new(-2733.0531146121098, 0.),
            Complex::new(-985.3438480599564, 0.),
            Complex::new(-3751.7953341772745, 0.),
            Complex::new(1770.7324644808666, -37.0257759674679),
        ]);
        assert!(
            (out.re + 0.0005352086907030319).hypot(out.im + 0.000004086322878714895) < 1e-10,
            "{out:?}"
        );
    });
}

fn c0_reference_sectors() {
    let parameters = [
        parse!("p1"),
        parse!("p2"),
        parse!("p3"),
        parse!("m1"),
        parse!("m2"),
        parse!("m3"),
        parse!("mu2"),
    ];
    let mapped = c0(
        [&parameters[0], &parameters[1], &parameters[2]],
        [&parameters[3], &parameters[4], &parameters[5]],
        &parameters[6],
    );
    let evaluator = mapped.evaluator(&parameters).unwrap();
    let mut evaluator = evaluator
        .map_coeff(&|value: &Complex<Rational>| Complex::new(value.re.to_f64(), value.im.to_f64()));
    let z = Complex::new(0.0, 0.0);
    let mut output = [z; 3];
    let cases = [
        (
            [0.0, 0.0, -3.0],
            [z, z, z],
            [
                Complex::new(-0.201_158_160_135_430_36, 0.0),
                Complex::new(0.366_204_096_222_703_3, 0.0),
                Complex::new(-1.0 / 3.0, 0.0),
            ],
        ),
        (
            [0.0, -2.0, -3.0],
            [z, z, z],
            [
                Complex::new(-0.363_247_973_447_190_3, 0.0),
                Complex::new(0.405_465_108_108_164_4, 0.0),
                z,
            ],
        ),
        (
            [-1.0, -2.0, -3.0],
            [z, z, z],
            [Complex::new(-1.254_138_667_411_437_2, 0.0), z, z],
        ),
        (
            [-1.0, -2.0, -3.0],
            [z, z, Complex::new(1.4, -0.08)],
            [
                Complex::new(-0.768_387_288_344_432_2, -0.013_967_008_239_948_784),
                z,
                z,
            ],
        ),
        (
            [-1.0, -2.0, -3.0],
            [z, Complex::new(0.7, -0.03), Complex::new(1.4, -0.08)],
            [
                Complex::new(-0.501_878_434_080_612_1, -0.013_567_392_613_106_676),
                z,
                z,
            ],
        ),
        (
            [-1.0, -2.0, -3.0],
            [
                Complex::new(0.5, -0.02),
                Complex::new(0.7, -0.03),
                Complex::new(1.4, -0.08),
            ],
            [
                Complex::new(-0.387_383_148_490_874_8, -0.012_089_012_569_460_895),
                z,
                z,
            ],
        ),
    ];

    for (momenta, masses, expected) in cases {
        let input = [
            Complex::new(momenta[0], 0.0),
            Complex::new(momenta[1], 0.0),
            Complex::new(momenta[2], 0.0),
            masses[0],
            masses[1],
            masses[2],
            Complex::new(1.0, 0.0),
        ];
        evaluator.evaluate(&input, &mut output);
        for (actual, expected) in output.into_iter().zip(expected) {
            let difference = actual - expected;
            assert!(
                difference.re.hypot(difference.im) < 2.0e-9,
                "momenta={momenta:?}, actual={actual:?}, expected={expected:?}"
            );
        }
    }
}

#[test]
fn massless_d0_infrared_sectors_match_reference() {
    let parameters = [
        parse!("p1"),
        parse!("p2"),
        parse!("p3"),
        parse!("p4"),
        parse!("p12"),
        parse!("p23"),
        parse!("mu2"),
    ];
    let p = [
        &parameters[0],
        &parameters[1],
        &parameters[2],
        &parameters[3],
        &parameters[4],
        &parameters[5],
    ];
    let cases = [
        (
            [0.0, 0.0, 0.0, 0.0, -5.0, -6.0],
            0,
            [
                -0.136_738_438_720_416_47,
                -0.226_746_492_110_810_33,
                0.133_333_333_333_333_33,
            ],
        ),
        (
            [0.0, 0.0, 0.0, -4.0, -5.0, -6.0],
            1,
            [
                -0.019_955_432_024_197_674,
                -0.134_326_868_036_151,
                0.066_666_666_666_666_67,
            ],
        ),
        (
            [0.0, 0.0, -3.0, -4.0, -5.0, -6.0],
            2,
            [
                0.037_500_182_171_644_69,
                -0.090_268_340_036_740_32,
                0.033_333_333_333_333_33,
            ],
        ),
        (
            [0.0, -2.0, 0.0, -4.0, -5.0, -6.0],
            5,
            [0.044_778_998_436_892_37, -0.120_159_621_816_574_5, 0.0],
        ),
        (
            [0.0, -2.0, -3.0, -4.0, -5.0, -6.0],
            3,
            [0.102_271_049_257_966_12, -0.060_079_810_908_287_25, 0.0],
        ),
    ];
    for (input, case, expected) in cases {
        let expressions = box_massless_low_case(p, case, &parameters[6]).into_coefficients();
        let views = expressions.iter().map(Atom::as_view).collect::<Vec<_>>();
        // The infrared formulas continue their logarithms on sheets.
        let mut map = FunctionMap::new();
        sheet_exact::register(&mut map);
        let evaluator = Atom::evaluator_multiple(&views, &parameters)
            .function_map(map.into())
            .build()
            .unwrap();
        let mut evaluator = evaluator.map_coeff(&|value: &Complex<Rational>| {
            Complex::new(value.re.to_f64(), value.im.to_f64())
        });
        let input = input
            .map(|value| Complex::new(value, 0.0))
            .into_iter()
            .chain([Complex::new(1.0, 0.0)])
            .collect::<Vec<_>>();
        let mut output = [Complex::new(0.0, 0.0); 3];
        evaluator.evaluate(&input, &mut output);
        for (actual, expected) in output.into_iter().zip(expected) {
            assert!(
                (actual.re - expected).abs() < 2.0e-9,
                "case={case}, actual={actual:?}"
            );
            assert!(actual.im.abs() < 2.0e-9);
        }
    }
}

#[test]
fn massive_d0_infrared_formulas_match_reference() {
    large_stack(massive_d0_infrared_reference_cases);
}

fn massive_d0_infrared_reference_cases() {
    let n = |value: i64| Atom::num(value);
    let mu = n(1);
    let m2 = n(2);
    let m3 = n(3);
    let m4 = n(4);
    let cases = [
        (
            box_one_mass_06(&n(-5), &n(-6), &m2, &mu),
            [0.016_077_113_897_630_62, -0.126_879_345_380_845_67, 0.05],
        ),
        (
            box_one_mass_07(&n(-3), &n(-5), &n(-6), &m2, &mu),
            [
                -0.024_114_641_085_640_903,
                -0.095_307_737_326_992_48,
                0.0375,
            ],
        ),
        (
            box_one_mass_08(&n(-3), &n(-4), &n(-5), &n(-6), &m2, &mu),
            [0.007_892_358_021_349_452, -0.059_178_090_353_290_43, 0.025],
        ),
        (
            box_one_mass_09(&n(-1), &n(-4), &n(-5), &n(-6), &m2, &mu),
            [0.018_400_361_421_865_807, -0.083_557_646_595_849_09, 0.0125],
        ),
        (
            box_one_mass_10(&n(-1), &n(-4), &n(-3), &n(-5), &n(-6), &m2, &mu),
            [0.039_897_243_462_062_45, -0.059_412_615_476_566_73, 0.0],
        ),
        (
            box_two_adjacent_11(&n(-3), &n(-5), &n(-6), &m2, &m3, &mu),
            [
                -0.012_934_410_692_339_639,
                -0.051_543_730_028_214_37,
                0.015_873_015_873_015_872,
            ],
        ),
        (
            box_two_adjacent_12(&n(-3), &n(-4), &n(-5), &n(-6), &m2, &m3, &mu),
            [
                0.016_303_862_354_627_337,
                -0.029_375_412_492_956_3,
                0.007_936_507_936_507_936,
            ],
        ),
        (
            box_two_adjacent_13(&n(-2), &n(-3), &n(-4), &n(-5), &n(-6), &m2, &m3, &mu),
            [0.015_350_152_013_953_197, -0.023_169_434_749_037_963, 0.0],
        ),
        (
            box_two_opposite_14(&n(-5), &n(-6), &m2, &m3, &mu),
            [0.094_666_988_330_808_07, -0.058_819_906_999_478_17, 0.0],
        ),
        (
            box_two_opposite_15(&n(-1), &m3, &n(-5), &n(-6), &m2, &m3, &mu),
            [0.121_951_458_233_307_58, -0.029_409_953_499_739_086, 0.0],
        ),
        (
            box_three_masses_16(&n(-2), &n(-3), &n(-5), &n(-6), &m2, &m3, &m4, &mu),
            [0.012_416_923_322_605_182, -0.016_350_857_574_315_11, 0.0],
        ),
    ];

    for (series, expected) in cases {
        let expressions = series.into_coefficients();
        let output = evaluate_series(&expressions, &[], &[]);
        for (actual, expected) in output.into_iter().zip(expected) {
            let difference = actual - Complex::new(expected, 0.0);
            assert!(
                difference.re.hypot(difference.im) < 3.0e-11,
                "actual={actual:?}, expected={expected:?}"
            );
        }
    }
}

#[test]
fn equal_mass_box_vacuum_limit() {
    large_stack(|| {
        let zero = Atom::num(0);
        let mass = parse!("vacuum_box_mass");
        let mu = parse!("vacuum_box_scale");
        let expression = d0([&zero; 6], [&mass; 4], &mu);
        let mut evaluator = expression
            .evaluator(&[mass, mu])
            .unwrap()
            .map_coeff(&|x: &Complex<Rational>| Complex::new(x.re.to_f64(), x.im.to_f64()));
        for (mass, scale) in [
            (Complex::new(1., 0.), 1.),
            (Complex::new(2., 0.), 7.),
            (Complex::new(0.5, -0.1), 3.),
        ] {
            let mut output = [Complex::new(0., 0.); 3];
            evaluator.evaluate(&[mass, Complex::new(scale, 0.)], &mut output);
            let expected = Complex::new(1. / 6., 0.) / (mass * mass);
            let error = output[0] - expected;
            assert!(error.re.hypot(error.im) < 2e-13);
            assert_eq!(output[1], Complex::new(0., 0.));
            assert_eq!(output[2], Complex::new(0., 0.));
        }
    });
}

#[test]
fn reciprocal_root_real_axis_preserves_its_sheet() {
    let q = parse!("real_axis_root_q");
    let (root, difference) = crate::triangle::r_function_body(&q);
    for (input, expected) in [
        (-7. / 6_f64.sqrt(), Complex::new(-6_f64.sqrt(), 0.)),
        (3., Complex::new((3. + 5_f64.sqrt()) / 2., 0.)),
        (1., Complex::new(0.5, -3_f64.sqrt() / 2.)),
        (-1., Complex::new(-0.5, -3_f64.sqrt() / 2.)),
    ] {
        let output = evaluate_series(
            &[root.clone(), difference.clone(), Atom::num(0)],
            std::slice::from_ref(&q),
            &[Complex::new(input, 0.)],
        );
        let error = output[0] - expected;
        assert!(error.re.hypot(error.im) < 2e-14);
        if input.abs() >= 2. {
            assert_eq!(output[0].im, 0.);
        } else {
            assert!(output[0].im < 0.);
        }
    }
}

#[test]
fn mapped_d0_routes_reference_sectors() {
    large_stack(d0_reference_sectors);
}

fn d0_reference_sectors() {
    let parameters = [
        parse!("p1"),
        parse!("p2"),
        parse!("p3"),
        parse!("p4"),
        parse!("p12"),
        parse!("p23"),
        parse!("m1"),
        parse!("m2"),
        parse!("m3"),
        parse!("m4"),
        parse!("mu2"),
    ];
    let mapped = d0_mapped(
        [
            &parameters[0],
            &parameters[1],
            &parameters[2],
            &parameters[3],
            &parameters[4],
            &parameters[5],
        ],
        [
            &parameters[6],
            &parameters[7],
            &parameters[8],
            &parameters[9],
        ],
        &parameters[10],
    );
    let evaluator = mapped.evaluator(&parameters).unwrap();
    let mut evaluator = evaluator
        .map_coeff(&|value: &Complex<Rational>| Complex::new(value.re.to_f64(), value.im.to_f64()));
    let z = Complex::new(0.0, 0.0);
    let r = |value: f64| Complex::new(value, 0.0);
    let cases = [
        (
            [0., 0., 2., 2., -5., -6., 0., 0., 0., 2., 1.],
            [0.016_077_113_897_630_62, -0.126_879_345_380_845_67, 0.05],
        ),
        (
            [0., 0., -3., 2., -5., -6., 0., 0., 0., 2., 1.],
            [
                -0.024_114_641_085_640_903,
                -0.095_307_737_326_992_48,
                0.0375,
            ],
        ),
        (
            [0., 0., -3., -4., -5., -6., 0., 0., 0., 2., 1.],
            [0.007_892_358_021_349_452, -0.059_178_090_353_290_43, 0.025],
        ),
        (
            [-1., 0., 2., -4., -5., -6., 0., 0., 0., 2., 1.],
            [0.018_400_361_421_865_807, -0.083_557_646_595_849_09, 0.0125],
        ),
        (
            [-1., 0., -3., -4., -5., -6., 0., 0., 0., 2., 1.],
            [0.039_897_243_462_062_45, -0.059_412_615_476_566_73, 0.],
        ),
        (
            [0., 2., -3., 3., -5., -6., 0., 0., 2., 3., 1.],
            [
                -0.012_934_410_692_339_639,
                -0.051_543_730_028_214_37,
                0.015_873_015_873_015_872,
            ],
        ),
        (
            [0., 2., -3., -4., -5., -6., 0., 0., 2., 3., 1.],
            [
                0.016_303_862_354_627_337,
                -0.029_375_412_492_956_3,
                0.007_936_507_936_507_936,
            ],
        ),
        (
            [0., -2., -3., -4., -5., -6., 0., 0., 2., 3., 1.],
            [0.015_350_152_013_953_197, -0.023_169_434_749_037_963, 0.],
        ),
        (
            [2., 2., 3., 3., -5., -6., 0., 2., 0., 3., 1.],
            [0.094_666_988_330_808_07, -0.058_819_906_999_478_17, 0.],
        ),
        (
            [-1., 2., 3., 3., -5., -6., 0., 2., 0., 3., 1.],
            [0.121_951_458_233_307_58, -0.029_409_953_499_739_086, 0.],
        ),
        (
            [2., -2., -3., 4., -5., -6., 0., 2., 3., 4., 1.],
            [0.012_416_923_322_605_182, -0.016_350_857_574_315_11, 0.],
        ),
        (
            [-1., -2., -3., -4., -5., -6., 0., 0., 0., 0., 1.],
            [0.230_720_567_898_758, 0., 0.],
        ),
        (
            [-1., -2., -3., -4., -5., -6., 0., 0., 0., 1.4, 1.],
            [0.159_028_182_108_436_3, 0., 0.],
        ),
        (
            [-1., -2., -3., -4., -5., -6., 0., 0., 0.7, 1.4, 1.],
            [0.116_778_930_088_213_91, 0., 0.],
        ),
        (
            [-1., -2., -3., -4., -5., -6., 0., 0.5, 0., 1.4, 1.],
            [0.107_217_099_967_768_14, 0., 0.],
        ),
        (
            [-1., -2., -3., -4., -5., -6., 0.3, 0.5, 0.7, 1.4, 1.],
            [0.064_264_380_690_262_37, 0., 0.],
        ),
    ];
    let mut output = [z; 3];
    for (input, expected) in cases {
        let input = input.map(r);
        evaluator.evaluate(&input, &mut output);
        for (actual, expected) in output.into_iter().zip(expected) {
            let difference = actual - r(expected);
            assert!(
                difference.re.hypot(difference.im) < 5.0e-9,
                "actual={actual:?}, expected={expected:?}"
            );
        }
    }
}
