//! Region selection preserves parametric values while removing decided IFs.
use oneloop::{ScalarIntegral, select_branch, select_branch_with};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::Coefficient,
    domains::float::{Complex, Float},
    id::{Pattern, ReplaceWith, Replacement},
    parse, symbol,
    transformer::Transformer,
};

fn on_stack(test: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(test)
        .unwrap()
        .join()
        .unwrap();
}

fn has_if(expression: &Atom) -> bool {
    let mut pending = vec![expression.as_view()];
    while let Some(view) = pending.pop() {
        match view {
            AtomView::Fun(fun) => {
                if fun.get_symbol() == Symbol::IF {
                    return true;
                }
                pending.extend(fun.iter());
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                pending.extend([base, exponent]);
            }
            AtomView::Mul(product) => pending.extend(product.iter()),
            AtomView::Add(sum) => pending.extend(sum.iter()),
            _ => {}
        }
    }
    false
}

#[test]
fn replacements_affect_only_conditions_including_nested_conditions() {
    on_stack(|| {
        let x = symbol!("branch_scope_x").to_atom();
        let y = symbol!("branch_scope_y").to_atom();
        let z = symbol!("branch_scope_z").to_atom();
        let nested = Symbol::IF.call((&y, &x + &y, &z));
        let condition = Symbol::IF.call((&x, &y, &z));
        let expression = Symbol::IF.call((condition, nested, &x - &z));
        let rules = [
            Replacement::new(x.clone(), 2),
            Replacement::new(y.clone(), 3),
        ];
        assert_eq!(select_branch(expression.as_view(), &rules), &x + &y);

        let unresolved = Symbol::IF.call((&x + &z, Symbol::IF.call((&y, &x, &z)), &y));
        assert_eq!(
            select_branch(&unresolved, &rules),
            Symbol::IF.call((&x + &z, &x, &y)),
            "an unresolved condition must not retain its substituted sample"
        );
        assert_eq!(select_branch(&(&x + &y), &rules), &x + &y);
    });
}

#[test]
fn visits_all_expression_positions_and_skips_inactive_arms() {
    on_stack(|| {
        let x = symbol!("branch_positions_x").to_atom();
        let y = symbol!("branch_positions_y").to_atom();
        let f = symbol!("branch_positions_f");
        let branch = Symbol::IF.call((&x, &y, 17));
        let expression = f.call((&branch,)) + (&branch + 1).pow(branch.clone()) * &branch;
        assert_eq!(
            select_branch(expression, &[Replacement::new(x.clone(), 1)]),
            f.call((&y,)) + (&y + 1).pow(y.clone()) * &y
        );
        let dead = Symbol::IF.call((&y, &y, 0));
        let expression = Symbol::IF.call((&x, &x, dead));
        let selected = select_branch_with(&expression, |condition| {
            assert_eq!(
                condition,
                x.as_view(),
                "inactive conditions must not be visited"
            );
            Ok::<_, ()>(Atom::num(1))
        })
        .unwrap();
        assert_eq!(selected, x);
        assert_eq!(
            select_branch_with(expression, |_| Err::<Atom, _>("replacement error")),
            Err("replacement error")
        );
    });
}

#[test]
fn native_replacement_patterns_and_rule_order_are_honored() {
    on_stack(|| {
        let x = symbol!("branch_pattern_x").to_atom();
        let y = symbol!("branch_pattern_y").to_atom();
        let f = symbol!("branch_pattern_f");
        let wildcard = symbol!("branch_pattern_a_").to_atom();
        let expression = Symbol::IF.call((f.call((&x + 1,)), &x, &y));
        let rules = [
            Replacement::new(f.call((&wildcard,)), 0),
            Replacement::new(f.call((&wildcard,)), 1),
        ];
        assert_eq!(select_branch(expression, &rules), y);
    });
}

#[test]
fn numeric_roots_complex_values_and_high_precision_conditions_use_native_evaluation() {
    on_stack(|| {
        let x = symbol!("branch_numeric_x").to_atom();
        let y = symbol!("branch_numeric_y").to_atom();
        for (condition, sample, expected) in [
            (x.clone(), Atom::num((1, 2)), x.clone()),
            (x.clone(), Atom::num(Float::with_val(256, 3.23)), x.clone()),
            (x.clone(), Atom::num(0), y.clone()),
            (x.clone(), parse!("1e-1000"), x.clone()),
            (x.clone(), parse!("1+1i"), x.clone()),
            (x.clone().sqrt() - 2, Atom::num(4), y.clone()),
            (Symbol::ABS.call((&x,)) - &x, Atom::num(4), y.clone()),
            (Symbol::ABS.call((&x,)) - &x, Atom::num(-4), x.clone()),
            (Symbol::CONJ.call((&x,)) - &x, parse!("1-2i"), x.clone()),
        ] {
            let expression = Symbol::IF.call((condition, &x, &y));
            assert_eq!(
                select_branch(expression, &[Replacement::new(x.clone(), sample)]),
                expected
            );
        }
    });
}

#[test]
fn all_deeply_nested_ifs_are_selected_without_a_fixed_depth_limit() {
    on_stack(|| {
        let condition = symbol!("branch_deep_condition").to_atom();
        let x = symbol!("branch_deep_x").to_atom();
        let mut expression = x.clone();
        for _ in 0..4096 {
            expression = Symbol::IF.call((&condition, expression, 0));
        }
        assert_eq!(
            select_branch(expression, &[Replacement::new(condition, 1)]),
            x
        );
    });
}

#[test]
fn deeply_nested_conditions_and_unknown_or_nonfinite_samples_are_preserved_safely() {
    on_stack(|| {
        let x = symbol!("branch_condition_depth_x").to_atom();
        let y = symbol!("branch_condition_depth_y").to_atom();
        let mut condition = y.clone();
        for _ in 0..1024 {
            condition = Symbol::IF.call((&x, condition, 0));
        }
        let expression = Symbol::IF.call((condition, &x + &y, &y));
        assert_eq!(
            select_branch(
                expression,
                &[
                    Replacement::new(x.clone(), 1),
                    Replacement::new(y.clone(), 2)
                ]
            ),
            &x + &y
        );
        let expression = Symbol::IF.call((&x, &x + &y, &y));
        assert_eq!(select_branch(&expression, &[]), expression);
        let safe_sample =
            Symbol::IF.call((Symbol::PI, 1, Atom::num(Coefficient::positive_infinity())));
        assert_eq!(
            select_branch(&expression, &[Replacement::new(x.clone(), safe_sample)]),
            &x + &y
        );
        assert_eq!(
            select_branch(
                &expression,
                &[Replacement::new(x.clone(), Atom::num(0).log())]
            ),
            expression
        );
        let expression = Symbol::IF.call((x.clone().sqrt() - 1, &x, &y));
        let nearly_one = parse!(&format!("1+1/10^{}", 400));
        assert_eq!(
            select_branch(&expression, &[Replacement::new(x.clone(), nearly_one)]),
            x
        );
        let nearly_one = Float::parse(&format!("1.{}1", "0".repeat(399)), Some(1600)).unwrap();
        assert_eq!(
            select_branch(
                expression,
                &[Replacement::new(x.clone(), Atom::num(nearly_one))]
            ),
            x
        );
        let expression = Symbol::IF.call(((&x + 1).sqrt() - 1, &x, &y));
        let tiny = Float::parse("1e-1000", Some(256)).unwrap();
        assert_eq!(
            select_branch(expression, &[Replacement::new(x.clone(), Atom::num(tiny))]),
            x
        );
    });
}

#[test]
fn float_rhs_callbacks_and_transformers_are_protected_before_condition_normalization() {
    on_stack(|| {
        let x = symbol!("branch_callback_x").to_atom();
        let y = symbol!("branch_callback_y").to_atom();
        let expression = Symbol::IF.call(((&x + 1).sqrt() - 1, &x, &y));
        let map = ReplaceWith::Map(Box::new(|_| {
            Atom::num(Float::parse("1e-1000", Some(256)).unwrap())
        }));
        assert_eq!(
            select_branch(&expression, &[Replacement::new(x.clone(), map)]),
            x
        );
        let transform = Pattern::Transformer(Box::new((
            Some(Pattern::Literal(Atom::num(0))),
            vec![Transformer::Map(Box::new(|_, _, out| {
                *out = Atom::num(Float::parse("1e-1000", Some(256)).unwrap());
                Ok(())
            }))],
        )));
        assert_eq!(
            select_branch(expression, &[Replacement::new(x.clone(), transform)]),
            x
        );
    });
}

#[test]
fn equal_mass_b0_becomes_if_free_and_remains_parametric_for_float_integer_and_fraction_rules() {
    on_stack(|| {
        let s = symbol!("branch_b0_psq"; Real).to_atom();
        let mass = symbol!("branch_b0_mass_squared"; Positive).to_atom();
        let expressions = ScalarIntegral::B0
            .get_expression(&[s.clone(), mass.clone(), mass.clone(), Atom::num(1)])
            .unwrap();
        assert!(has_if(&expressions.coefficients()[0]));
        let discriminant = (&s * &s - Atom::num(4) * &s * &mass).sqrt();
        let root_plus = (&s + &discriminant) / (&s * 2);
        let root_minus = (&s - &discriminant) / (&s * 2);
        let root_term = |root: &Atom| {
            Atom::num(-1) + (Atom::num(1) - root) * (Atom::num(1) - Atom::num(1) / root).log()
        };
        let below_threshold = -mass.log() - root_term(&root_plus) - root_term(&root_minus);
        for (sample, sample_float) in [
            (Atom::num(Float::with_val(128, 3.23)), 3.23_f64),
            (Atom::num(3), 3.0),
            (Atom::num((1, 2)), 0.5),
            (Atom::num(-2), -2.0),
            (Atom::num(5), 5.0),
            (Atom::num(0), 0.0),
            (Atom::num(4), 4.0),
        ] {
            let rules = [
                Replacement::new(s.clone(), sample),
                Replacement::new(mass.clone(), 1),
            ];
            let selected = expressions
                .coefficients()
                .each_ref()
                .map(|coefficient| select_branch(coefficient, &rules));
            assert!(
                selected.iter().all(|coefficient| !has_if(coefficient)),
                "sample {sample_float}: {selected:?}"
            );
            assert_eq!(selected[1], Atom::num(1));
            assert_eq!(selected[2], Atom::num(0));
            if sample_float > 0.0 && sample_float < 4.0 {
                // The selected formula is the parametric root/log form. Its roots
                // are evaluated without cancellation, so compare the closed forms
                // as functions at parametric points below threshold.
                for (momentum_value, mass_value) in [(sample_float, 1.0), (1.3, 0.9), (0.4, 2.0)] {
                    let rules = [
                        Replacement::new(
                            s.clone(),
                            Atom::num(Float::with_val(128, momentum_value)),
                        ),
                        Replacement::new(mass.clone(), Atom::num(Float::with_val(128, mass_value))),
                    ];
                    let value = |expression: &Atom| -> Complex<f64> {
                        expression
                            .replace_multiple(&rules)
                            .evaluate::<Atom, _>(&std::collections::HashMap::new())
                            .unwrap()
                    };
                    let (actual, expected) = (value(&selected[0]), value(&below_threshold));
                    assert!(
                        (actual.re - expected.re).hypot(actual.im - expected.im) < 1e-12,
                        "the selected parametric root/log formula at ({momentum_value}, {mass_value}): {actual:?} vs {expected:?}"
                    );
                }
            }
            let variables = selected[0].get_all_symbols(false);
            assert!(
                variables.contains(&mass.get_symbol().unwrap()),
                "mass must remain parametric"
            );
            if sample_float != 0.0 {
                assert!(
                    variables.contains(&s.get_symbol().unwrap()),
                    "momentum must remain parametric"
                );
            }
            let nearby = if sample_float == 0.0 {
                (0.0, 1.25)
            } else if sample_float == 4.0 {
                (5.0, 1.25)
            } else {
                (sample_float + sample_float.signum() * 0.17, 1.0)
            };
            for (momentum_value, mass_value) in [(sample_float, 1.0), nearby] {
                let mut reference = [Complex::new(0.0, 0.0); 3];
                oneloop::evaluate(
                    ScalarIntegral::B0,
                    &[
                        Complex::new(momentum_value, 0.0),
                        Complex::new(mass_value, 0.0),
                        Complex::new(mass_value, 0.0),
                        Complex::new(1.0, 0.0),
                    ],
                    &mut reference,
                )
                .unwrap();
                let check_rules = [
                    Replacement::new(s.clone(), Atom::num(Float::with_val(128, momentum_value))),
                    Replacement::new(mass.clone(), Atom::num(Float::with_val(128, mass_value))),
                ];
                let actual: Complex<f64> = selected[0]
                    .replace_multiple(&check_rules)
                    .evaluate::<Atom, _>(&std::collections::HashMap::new())
                    .unwrap();
                assert!(
                    (actual.re - reference[0].re).abs() < 1e-12,
                    "sample {sample_float}, check ({momentum_value}, {mass_value}): {actual:?} vs {reference:?}"
                );
                assert!(
                    (actual.im - reference[0].im).abs() < 1e-12,
                    "sample {sample_float}, check ({momentum_value}, {mass_value}): {actual:?} vs {reference:?}"
                );
            }
        }
    });
}
