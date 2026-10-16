//! Exercise registration independently of the prepared integral caches.
#[path = "../src/definitions.rs"]
mod definitions;

use definitions::FunctionMap;
use symbolica::{
    evaluate::{FunctionRegistrationOptions, InliningPolicy},
    prelude::*,
};

#[test]
fn interpreted_log_retains_accuracy_with_a_weak_tiny_component() {
    let x = symbol!("mixed_precision_log_argument").to_atom();
    let map = FunctionMap::new();
    let (expressions, parameters) =
        map.normalize_inputs(&[x.log(), Symbol::ABS.call((&x,))], &[x], false);
    let mut evaluator = Atom::evaluator_multiple(&expressions, &parameters)
        .function_map(map.as_symbolica().clone())
        .build()
        .unwrap()
        .map_coeff_with_prec(
            &|c| Complex::new(c.re.to_multi_prec_float(128), c.im.to_multi_prec_float(128)),
            128,
        );
    let mut tiny = Float::with_val(128, 0.5).pow(120);
    tiny.set_prec(2);
    for sign in [-1, 1] {
        let large = Float::with_val(128, 2 * sign);
        for z in [
            Complex::new(large.clone(), tiny.clone()),
            Complex::new(tiny.clone(), large),
        ] {
            let mut output =
                core::array::from_fn::<_, 2, _>(|_| Complex::new(Float::new(128), Float::new(128)));
            evaluator.evaluate(std::slice::from_ref(&z), &mut output);
            let [value, norm] = output;
            assert!(norm.re.prec() >= 120, "{norm:?}");
            assert!((norm.re - Float::with_val(128, 2)).norm().to_f64() < 1e-35);
            let expected = Float::with_val(128, 2).log();
            assert!(value.re.prec() >= 120, "{value:?}");
            assert!((value.re - expected).norm().to_f64() < 1e-35);
            if z.re > Float::new(128) && z.re.norm() > z.im.norm() {
                assert!(
                    value.im.prec() <= 4,
                    "the uncertain phase must stay uncertain"
                );
            }
        }
    }
}

#[test]
fn replacing_a_helper_keeps_tagged_callees_and_the_original_map() {
    let x = symbol!("replaced_helper_x");
    let f = symbol!("replaced_helper_f");
    let g = symbol!("replaced_helper_g");
    let mut original = FunctionMap::new();
    let options = FunctionRegistrationOptions::new().inlining(InliningPolicy::Never);
    original
        .add_tagged_function_with_options(
            f,
            vec![Atom::num(2)],
            vec![x],
            x.to_atom() + Symbol::PI,
            options.clone(),
        )
        .unwrap();
    original
        .add_function_with_options(g, vec![x], f.call((2, x)), options)
        .unwrap();
    let replaced = original.with_replaced_body(g, 2 * f.call((2, x))).unwrap();
    assert!(original.with_replaced_body(f, Atom::num(0)).is_err());
    for (map, factor) in [(original, 1.), (replaced, 2.)] {
        let mut evaluator = g
            .call((x,))
            .evaluator(&[x.to_atom()])
            .function_map(map.into())
            .build()
            .unwrap()
            .map_coeff(&|c| c.re.to_f64());
        assert_eq!(
            evaluator.evaluate_single(&[1.]),
            factor * (1. + std::f64::consts::PI)
        );
    }
}

#[test]
fn interpreted_log_retains_the_precise_small_real_part_near_one() {
    let x = symbol!("mixed_precision_near_unit_argument").to_atom();
    let map = FunctionMap::new();
    let (expressions, parameters) = map.normalize_inputs(&[x.log()], &[x], false);
    let exact = Atom::evaluator_multiple(&expressions, &parameters)
        .function_map(map.as_symbolica().clone())
        .build()
        .unwrap();
    for bits in [128, 512, 3456] {
        let mut evaluator = exact.clone().map_coeff_with_prec(
            &|c| {
                Complex::new(
                    c.re.to_multi_prec_float(bits),
                    c.im.to_multi_prec_float(bits),
                )
            },
            bits,
        );
        // The axis value is exact despite its low working precision. The
        // accurately known tiny component determines the small real logarithm.
        let tiny = Float::with_val(bits, 0.5).pow(bits as u64);
        let z = Complex::new(Float::with_val(2, 1), tiny.clone());
        let value = evaluator.evaluate_single(&[z]);
        let expected = tiny.clone() * tiny / Float::with_val(bits, 2);
        assert!(!value.re.is_fully_zero(), "{value:?}");
        assert!(value.re.prec() >= bits - 8, "{value:?}");
        let error = (value.re / expected - Float::with_val(bits, 1)).norm();
        assert!(error < Float::with_val(bits, 0.5).pow((bits - 8) as u64));
    }
}

#[test]
fn explicit_constants_preserve_nested_tagged_calls_and_inspection() {
    let x = symbol!("constant_forward_x");
    let y = symbol!("constant_forward_y");
    let f = symbol!("constant_forward_f");
    let g = symbol!("constant_forward_g");
    let options = FunctionRegistrationOptions::new().inlining(InliningPolicy::Never);
    let body = parse!("pi + polylog(2, constant_forward_x)");
    let mut map = FunctionMap::new();
    map.add_tagged_function_with_options(
        f,
        vec![Atom::num(2)],
        vec![x],
        body.clone(),
        options.clone(),
    )
    .unwrap();
    map.add_function_with_options(
        g,
        vec![y, x],
        f.call((2, y)) + f.call((2, y.to_atom() / 2)),
        options,
    )
    .unwrap();
    let inspected = f.call((2, y));
    let (tags, args, original) = map.get_definition(inspected.as_view()).unwrap();
    assert_eq!(tags, 1);
    assert_eq!(args.len(), 1);
    assert_eq!(original, &body);
    assert!(map.get_definition(f.call((3, y)).as_view()).is_none());

    let expression = g.call((y, 7));
    let parameters = [y.to_atom()];
    let (expressions, normalized) = map.normalize_inputs(&[expression], &parameters, false);
    let exact = Atom::evaluator_multiple(&expressions, &normalized)
        .function_map(map.into())
        .build()
        .unwrap();
    // Export used to panic after a lifted pi changed external-function indices.
    assert_eq!(exact.export_instructions().output_count, 1);
    let mut actual = exact.map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    let bytes = bincode::encode_to_vec(&actual, bincode::config::standard()).unwrap();
    let (mut restored, consumed): (
        symbolica::evaluate::ExpressionEvaluator<Complex<f64>>,
        usize,
    ) = bincode::decode_from_slice(&bytes, bincode::config::standard()).unwrap();
    assert_eq!(consumed, bytes.len());
    assert_eq!(restored.get_input_len(), 1);
    assert_eq!(restored.get_output_len(), 1);
    let mut compiled = restored
        .jit_compile(
            symbolica::evaluate::JITCompilationSettings::new()
                .optimization_level(2)
                .direct_translation(false)
                .with_option("fastmath", "false")
                .with_option("fast_complex", "true")
                .with_option("use_threads", "false")
                .with_option("use_simd", "false")
                .with_option("simd_branch", "false")
                .with_option("enable_simd512", "false"),
        )
        .unwrap();
    let mut expected =
        parse!("2*pi + polylog(2,constant_forward_y) + polylog(2,constant_forward_y/2)")
            .evaluator(&parameters)
            .build()
            .unwrap()
            .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    for input in [Complex::new(0.25, 0.), Complex::new(-0.5, 0.125)] {
        let expected = expected.evaluate_single(&[input]);
        let mut output = [Complex::new(0., 0.)];
        compiled.evaluate(&[input], &mut output);
        for actual in [
            actual.evaluate_single(&[input]),
            restored.evaluate_single(&[input]),
            output[0],
        ] {
            let delta = actual - expected;
            assert!(delta.re.hypot(delta.im) < 1e-13);
        }
    }
}

#[test]
fn functions_capture_global_inputs_and_shadow_only_their_own_parameters() {
    let x = symbol!("forward_scope_x");
    let f = symbol!("forward_scope_f");
    let g = symbol!("forward_scope_g");
    let mut map = FunctionMap::new();
    map.add_function(f, Vec::<Symbol>::new(), x.to_atom() + Symbol::PI)
        .unwrap();
    // Since Symbolica 4cd26c33, ordinary functions see global inputs and
    // their own formals. Caller-local capture now requires add_aliases.
    map.add_function(g, vec![x], f.call(()) + x).unwrap();
    let expressions = [g.call((2,)), g.call((3,))];
    let mut evaluator = Atom::evaluator_multiple(&expressions, &[x.to_atom()])
        .function_map(map.into())
        .build()
        .unwrap()
        .map_coeff(&|c| c.re.to_f64());
    let mut output = [0.; 2];
    for global in [5., -7.] {
        evaluator.evaluate(&[global], &mut output);
        assert_eq!(
            output,
            [
                global + std::f64::consts::PI + 2.,
                global + std::f64::consts::PI + 3.
            ]
        );
    }
}

#[test]
fn jit_lowering_preserves_complex_branches_in_nested_and_restored_evaluators() {
    let x = symbol!("jit_lowering_x");
    let root = symbol!("jit_lowering_root");
    let z = x.to_atom();
    let mut map = FunctionMap::new();
    map.add_function_with_options(
        root,
        vec![x],
        z.sqrt() + Symbol::PI,
        FunctionRegistrationOptions::new().inlining(InliningPolicy::Never),
    )
    .unwrap();
    let expressions = [
        root.call((x,)),
        z.pow(Atom::num((-1, 2))),
        z.log(),
        z.pow(-2),
        Symbol::IF.call((&z, 7, 11)),
        Symbol::IF.call((&z - z.conj(), 7, 11)),
        &z * z.conj(),
    ];
    let (plain, plain_parameters) =
        map.normalize_inputs(&expressions, std::slice::from_ref(&z), false);
    let mut reference = Atom::evaluator_multiple(&plain, &plain_parameters)
        .function_map(map.as_symbolica().clone())
        .direct_translation(true)
        .horner_iterations(0)
        .build()
        .unwrap()
        .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    let (lowered, parameters) = map.normalize_inputs(&expressions, &[z], true);
    let source = Atom::evaluator_multiple(&lowered, &parameters)
        .function_map(map.as_jit_symbolica().clone())
        .direct_translation(true)
        .horner_iterations(0)
        .build()
        .unwrap()
        .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    let bytes = bincode::encode_to_vec(&source, bincode::config::standard()).unwrap();
    let (restored, consumed): (
        symbolica::evaluate::ExpressionEvaluator<Complex<f64>>,
        usize,
    ) = bincode::decode_from_slice(&bytes, bincode::config::standard()).unwrap();
    assert_eq!(consumed, bytes.len());
    let inputs = [
        (4., 0.),
        (4., -0.),
        (-4., 0.),
        (-4., -0.),
        (0., 4.),
        (-0., 4.),
        (0., -4.),
        (-0., -4.),
        (1e-20, 1.),
        (1e-20, -1.),
        (0.7, -0.03),
        (1. / 7.31, -1. / 13.71),
    ]
    .map(|(re, im)| Complex::new(re, im));
    for source in [source, restored] {
        let mut compiled = source
            .jit_compile(
                symbolica::evaluate::JITCompilationSettings::new()
                    .optimization_level(2)
                    .direct_translation(false)
                    .with_option("fastmath", "false")
                    .with_option("fast_complex", "true")
                    .with_option("use_threads", "false")
                    .with_option("use_simd", "false"),
            )
            .unwrap();
        let mut batch = vec![Complex::new(0., 0.); inputs.len() * expressions.len()];
        compiled.batch_evaluate(&inputs, &mut batch, inputs.len());
        for (input, batch) in inputs.iter().zip(batch.chunks_exact(expressions.len())) {
            let mut expected = vec![Complex::new(0., 0.); expressions.len()];
            let mut scalar = expected.clone();
            reference.evaluate(&[*input], &mut expected);
            compiled.evaluate(&[*input], &mut scalar);
            for actual in [scalar.as_slice(), batch] {
                assert_eq!(
                    actual.last().unwrap().im,
                    0.,
                    "a conjugate product must be exactly real"
                );
                for (actual, expected) in actual.iter().zip(&expected) {
                    let delta = *actual - *expected;
                    assert!(
                        delta.re.hypot(delta.im) < 1e-13,
                        "{input:?}: {actual:?} != {expected:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn helpers_preserve_the_original_floating_point_sum_order() {
    let x = symbol!("ordered_argument");
    let (f, g, h) = symbol!("ordered_first", "ordered_second", "ordered_third");
    let mut map = FunctionMap::new();
    for (symbol, factor) in [
        (h, 1),
        (g, -10_000_000_000_000_000i64),
        (f, 10_000_000_000_000_000i64),
    ] {
        map.add_function_with_options(
            symbol,
            vec![x],
            Atom::num(factor) * x,
            FunctionRegistrationOptions::new().inlining(InliningPolicy::Never),
        )
        .unwrap();
    }
    let sum = f.call((x,)) + g.call((x,)) + h.call((x,));
    let (expressions, parameters) = map.normalize_inputs(&[sum], &[x.to_atom()], false);
    let mut evaluator = Atom::evaluator_multiple(&expressions, &parameters)
        .function_map(map.as_symbolica().clone())
        .direct_translation(true)
        .horner_iterations(0)
        .build()
        .unwrap()
        .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    assert_eq!(
        evaluator.evaluate_single(&[Complex::new(1., 0.)]),
        Complex::new(1., 0.)
    );
}

#[test]
fn jit_square_root_callback_preserves_extreme_finite_scales() {
    let z = symbol!("jit_sqrt_scale_argument").to_atom();
    let map = FunctionMap::new();
    let (expressions, parameters) = map.normalize_inputs(&[z.sqrt()], &[z], true);
    let exact = Atom::evaluator_multiple(&expressions, &parameters)
        .function_map(map.as_jit_symbolica().clone())
        .build()
        .unwrap();
    let mut evaluator = exact
        .jit_compile::<Complex<f64>>(
            symbolica::evaluate::JITCompilationSettings::new()
                .optimization_level(2)
                .direct_translation(false)
                .with_option("fastmath", "false")
                .with_option("fast_complex", "true")
                .with_option("use_simd", "false")
                .with_option("use_threads", "false"),
        )
        .unwrap();
    for (input, expected) in [
        (Complex::new(1e300, 0.), Complex::new(1e150, 0.)),
        (Complex::new(1e-300, 0.), Complex::new(1e-150, 0.)),
        (Complex::new(0., -2e300), Complex::new(1e150, -1e150)),
        (Complex::new(0., -2e-300), Complex::new(1e-150, -1e-150)),
    ] {
        let mut output = [Complex::new(f64::NAN, f64::NAN)];
        evaluator.evaluate(&[input], &mut output);
        for (actual, expected) in [(output[0].re, expected.re), (output[0].im, expected.im)] {
            assert!(actual.is_finite(), "sqrt({input:?}) = {:?}", output[0]);
            if expected == 0. {
                assert_eq!(actual, 0.);
            } else {
                assert!((actual / expected - 1.).abs() < 1e-14);
            }
        }
    }
}
