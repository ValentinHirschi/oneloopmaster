//! Exact but compact sheet expressions used by the mapped box sectors.
use crate::definitions::FunctionMap;

use core::ops::{Div, Mul};

use symbolica::atom::{Atom, AtomCore, Symbol};

fn choose(c: &Atom, yes: Atom, no: Atom) -> Atom {
    Symbol::IF.call((c, yes, no))
}

pub(crate) fn i() -> Atom {
    Atom::num(-1).sqrt()
}

pub(crate) fn re(z: &Atom) -> Atom {
    symbolica::symbol!("__olo_real_part").call((z,))
}

pub(crate) fn im(z: &Atom) -> Atom {
    symbolica::symbol!("__olo_imaginary_part").call((z,))
}

/// Lower-lip square root with a bound argument, shared by C0/D0 formulas.
pub(crate) fn sqrt_lower(z: &Atom) -> Atom {
    symbolica::symbol!("__olo_sqrt_lower").call((z,))
}

fn zero(x: &Atom) -> Atom {
    choose(x, Atom::num(0), Atom::num(1))
}

pub(crate) fn sign_nonnegative(x: &Atom) -> Atom {
    // Evaluate x once: expanding a composite x in x-|x| can leave a rounding
    // residual, which native if interprets as nonzero.
    symbolica::symbol!("__olo_sign_nonnegative").call((x,))
}

pub(crate) fn negative(x: &Atom) -> Atom {
    (Atom::num(1) - sign_nonnegative(x)) / 2
}

fn xor(a: &Atom, b: &Atom) -> Atom {
    a + b - Atom::num(2) * a * b
}

fn physical(value: &Atom, odd: &Atom) -> Atom {
    (Atom::num(1) - Atom::num(2) * odd) * value
}

fn log_body(value: &Atom, phase: &Atom, odd: &Atom) -> Atom {
    let i_pi = i() * Symbol::PI.to_atom();
    let imaginary = im(value);
    let real_axis = zero(&imaginary);
    let real = value.log() + &i_pi * phase;
    let adjusted_phase = phase + odd * sign_nonnegative(&imaginary);
    let complex = physical(value, odd).log() + i_pi * adjusted_phase;
    choose(&real_axis, real, complex)
}

fn dilog_body(value: &Atom, phase: &Atom, odd: &Atom, at_one: Atom, at_zero: Atom) -> Atom {
    let i_pi = i() * Symbol::PI.to_atom();
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    let even = phase - odd;
    let argument = physical(value, odd);
    let nonunit = &argument - 1;
    let inverted = negative(&(Atom::num(1) - value * value.conj()));
    let normalized = choose(&inverted, Atom::num(1) / &argument, argument);
    let logarithm = choose(&inverted, -value.log(), value.log());
    let even = choose(&inverted, -&even, even);
    let odd = choose(&inverted, -odd, odd.clone());
    let sheet_logarithm = &logarithm + &i_pi * (&even + &odd);
    let log_one_minus = (Atom::num(1) - &normalized).log();
    let below_half = negative(&(Atom::num(2) * re(&normalized) - 1));
    let below =
        &pi_squared / 6 - crate::dilog_atom(&normalized) - &log_one_minus * &sheet_logarithm;
    let above = crate::dilog_atom(Atom::num(1) - &normalized) - &log_one_minus * &i_pi * even;
    let base = choose(&below_half, below, above);
    let ordinary = choose(
        &inverted,
        -&base - &sheet_logarithm * &sheet_logarithm / 2,
        base,
    );
    choose(value, choose(&nonunit, ordinary, at_one), at_zero)
}

fn log_symbol() -> Symbol {
    symbolica::symbol!("__olo_sheet_log")
}

fn dilog_symbol() -> Symbol {
    symbolica::symbol!("__olo_sheet_dilog")
}

fn log_ratio_symbol() -> Symbol {
    symbolica::symbol!("__olo_sheet_log_over_one_minus")
}

fn difference_symbol() -> Symbol {
    symbolica::symbol!("__olo_sheet_dilog_difference")
}

fn sheet_parameters() -> [Symbol; 3] {
    [
        symbolica::symbol!("__olo_sheet_value"),
        symbolica::symbol!("__olo_sheet_phase"),
        symbolica::symbol!("__olo_sheet_odd"),
    ]
}

fn log_ratio_body(arguments: &[Atom; 3], neighborhood: Atom) -> Atom {
    let difference = Atom::num(1) - physical(&arguments[0], &arguments[2]);
    let logarithm = log_symbol().call(arguments.as_slice());
    let quotient = &logarithm / &difference;
    // A nonzero sheet logarithm at one is a pole, not a removable limit.
    let at_one = choose(&logarithm, quotient.clone(), neighborhood);
    choose(&difference, quotient, at_one)
}

/// Generate the complete local polynomial before Dualizer selects a zero primal.
/// The requested order is the largest total degree of the actual Dual shape.
pub(crate) fn with_neighborhoods(map: &FunctionMap, order: usize) -> Result<FunctionMap, String> {
    let variable = symbolica::symbol!("__olo_log_ratio_neighborhood");
    let x = variable.to_atom();
    let local = (x.log() / (Atom::num(1) - &x))
        .series(
            variable,
            1,
            symbolica::poly::series::SeriesDepth::absolute(order),
        )
        .map_err(|error| error.to_string())?
        .to_atom();
    let arguments = sheet_parameters().map(Symbol::to_atom);
    let local = local
        .replace(variable)
        .with(physical(&arguments[0], &arguments[2]));
    let mut map = map.with_replaced_body(log_ratio_symbol(), log_ratio_body(&arguments, local))?;
    let arguments = product_parameters().map(Symbol::to_atom);
    let parts = product_parts(&arguments, true, false);
    for (component, body) in parts.into_iter().enumerate() {
        map = map.with_replaced_body(product_symbol(true, component), body)?;
    }
    with_dilog_neighborhoods(map, order)
}

fn difference_parameters() -> [Symbol; 6] {
    [
        "__olo_dd_v1",
        "__olo_dd_p1",
        "__olo_dd_o1",
        "__olo_dd_v2",
        "__olo_dd_p2",
        "__olo_dd_o2",
    ]
    .map(|s| symbolica::symbol!(s))
}

fn difference_sheets(arguments: &[Atom; 6]) -> [SheetAtom; 2] {
    core::array::from_fn(|i| SheetAtom {
        value: arguments[3 * i].clone(),
        phase: arguments[3 * i + 1].clone(),
        odd: arguments[3 * i + 2].clone(),
    })
}

fn difference_body(first: &SheetAtom, second: &SheetAtom, coincident: Atom) -> Atom {
    let difference = first.value() - second.value();
    choose(
        &difference,
        (first.dilog() - second.dilog()) / &difference,
        coincident,
    )
}

fn with_dilog_neighborhoods(mut map: FunctionMap, order: usize) -> Result<FunctionMap, String> {
    let depth = order
        .checked_add(1)
        .ok_or("dilogarithm neighborhood order overflow")?;
    let delta = symbolica::symbol!("__olo_dilog_offset");
    let center = symbolica::symbol!("__olo_dilog_center");
    let other = symbolica::symbol!("__olo_dilog_other_offset");
    let h = delta.to_atom();
    let y = center.to_atom();
    // Series and cancellation are performed before specializing either argument.
    // In particular, the joint unit limit needs the native polylog-at-zero series.
    let local_numerator = crate::dilog_atom(-&h)
        .series(delta, 0, depth)
        .map_err(|error| error.to_string())?
        .to_atom();
    let unit_difference = ((&local_numerator - local_numerator.replace(delta).with(other))
        / (&h - other.to_atom()))
    .cancel();
    // Keep Li2's argument as the symbolic center until log(1-y) has been
    // replaced by the continued sheet logarithm below.
    let numerator = crate::dilog_atom(&y + &h)
        .series(delta, 0, depth)
        .map_err(|error| error.to_string())?
        .to_atom();
    let general = (-(numerator - crate::dilog_atom(&y)) / &h)
        .cancel()
        .series(delta, 0, order)
        .map_err(|error| error.to_string())?;
    let constant = general
        .coefficient(0.into())
        .ok_or("missing dilogarithm constant")?;
    let correction = general.to_atom() - constant;

    let arguments = sheet_parameters().map(Symbol::to_atom);
    let value = physical(&arguments[0], &arguments[2]);
    let logarithm = log_symbol().call(arguments.as_slice());
    let unit = choose(
        &logarithm,
        &logarithm / (Atom::num(1) - &value),
        local_numerator.replace(delta).with(&value - 1),
    );
    // Li2(1-z) is finite at z=0, but its derivatives are nonanalytic there.
    // Keep the native primitive so its singular derivatives remain visible.
    map = map.with_replaced_body(
        dilog_symbol(),
        dilog_body(
            &arguments[0],
            &arguments[1],
            &arguments[2],
            unit,
            crate::dilog_atom(Atom::num(1) - &value),
        ),
    )?;

    let arguments = difference_parameters().map(Symbol::to_atom);
    let [first, second] = difference_sheets(&arguments);
    let difference = first.value() - second.value();
    let undefined = (first.dilog() - second.dilog()) / &difference;
    // The positive offset powers vanish exactly at coincidence, preserving
    // the existing scalar primal while keeping both arguments' derivatives.
    let regular = second.log_over_one_minus()
        + correction
            .replace((Atom::num(1) - &y).log())
            .with(second.log())
            .replace(center)
            .with(Atom::num(1) - second.value())
            .replace(delta)
            .with(-&difference);
    let unit = unit_difference
        .replace(delta)
        .with(first.value() - 1)
        .replace(other)
        .with(second.value() - 1);
    let coincident = choose(
        &(first.phase.clone() - &second.phase),
        undefined.clone(),
        choose(
            &(first.odd.clone() - &second.odd),
            undefined.clone(),
            choose(
                &second.value(),
                choose(
                    &(second.value() - 1),
                    regular,
                    choose(&second.log(), undefined.clone(), unit),
                ),
                undefined,
            ),
        ),
    );
    map.with_replaced_body(
        difference_symbol(),
        difference_body(&first, &second, coincident),
    )
}

fn product_parameters() -> [Symbol; 6] {
    [
        "__olo_q_v1",
        "__olo_q_p1",
        "__olo_q_o1",
        "__olo_q_v2",
        "__olo_q_p2",
        "__olo_q_o2",
    ]
    .map(|s| symbolica::symbol!(s))
}

fn product_parts(arguments: &[Atom; 6], divide: bool, scalar_zero_ratio: bool) -> [Atom; 3] {
    let [left, left_phase, left_odd, right, right_phase, right_odd] = arguments;
    let raw = if divide {
        let quotient = left / right;
        // The quotient is linear in its numerator. Its exact local expression
        // preserves a primal value of one at equal arguments without erasing
        // independent jets or introducing reciprocal-multiplication roundoff.
        let centered = quotient
            .series(product_parameters()[0], right, 1)
            .expect("a linear quotient has an exact local expansion")
            .to_atom();
        let zero_ratio = if scalar_zero_ratio {
            // Preserve only the historical scalar 0/0 convention.
            Atom::num(1)
        } else {
            // Without a specified common vanishing factor, a 0/0 jet has no
            // unique continuation. The derivative constructor leaves it nonfinite.
            quotient.clone()
        };
        choose(
            &(left - right),
            quotient,
            choose(right, centered, zero_ratio),
        )
    } else {
        left * right
    };
    let n = negative(&re(&raw));
    let phase = if divide {
        left_phase - right_phase - &n * sign_nonnegative(&im(right))
    } else {
        left_phase + right_phase + &n * sign_nonnegative(&im(right))
    };
    [
        (1 - Atom::num(2) * &n) * raw,
        phase,
        xor(&xor(left_odd, right_odd), &n),
    ]
}

pub(crate) fn register(function_map: &mut FunctionMap) {
    register_impl(function_map, false);
}

pub(crate) fn register_for_inspection(function_map: &mut FunctionMap) {
    register_impl(function_map, true);
}

fn register_impl(function_map: &mut FunctionMap, grouped_predicates: bool) {
    let sign_arg = symbolica::symbol!("__olo_sign_argument");
    let x = sign_arg.to_atom();
    let component = if grouped_predicates {
        super::grouped(&x)
    } else {
        x.clone()
    };
    function_map
        .add_function_with_options(
            symbolica::symbol!("__olo_sqrt_lower"),
            vec![sign_arg],
            super::physical_sqrt(&x),
            super::native_function_options(),
        )
        .unwrap();
    function_map
        .add_function(
            symbolica::symbol!("__olo_real_part"),
            vec![sign_arg],
            (&component + component.conj()) / 2,
        )
        .unwrap();
    function_map
        .add_function(
            symbolica::symbol!("__olo_imaginary_part"),
            vec![sign_arg],
            (&component - component.conj()) / (Atom::num(2) * i()),
        )
        .unwrap();
    // In a FunctionMap, x is a bound parameter evaluated once. In an inlined
    // expression, preserve that arithmetic grouping so x-|x| cannot acquire a
    // nonzero rounding residue from differently flattened copies of x.
    let sign_value = if grouped_predicates {
        super::grouped(&x)
    } else {
        x.clone()
    };
    function_map
        .add_function(
            symbolica::symbol!("__olo_sign_nonnegative"),
            vec![sign_arg],
            choose(
                &(&sign_value - Symbol::ABS.call((&sign_value,))),
                Atom::num(-1),
                Atom::num(1),
            ),
        )
        .unwrap();
    let args = product_parameters();
    let x = args.map(Symbol::to_atom);
    for divide in [false, true] {
        let parts = product_parts(&x, divide, true);
        for (k, body) in parts.into_iter().enumerate() {
            function_map
                .add_function_with_options(
                    product_symbol(divide, k),
                    args.to_vec(),
                    body,
                    super::native_function_options(),
                )
                .unwrap();
        }
    }
    let variables = sheet_parameters();
    let atoms = variables.map(Symbol::to_atom);
    function_map
        .add_function_with_options(
            log_symbol(),
            variables.to_vec(),
            log_body(&atoms[0], &atoms[1], &atoms[2]),
            super::native_function_options(),
        )
        .expect("unique sheet logarithm");
    function_map
        .add_function_with_options(
            log_ratio_symbol(),
            variables.to_vec(),
            log_ratio_body(&atoms, Atom::num(-1)),
            super::native_function_options(),
        )
        .expect("unique sheet logarithm ratio");
    function_map
        .add_function_with_options(
            dilog_symbol(),
            variables.to_vec(),
            dilog_body(
                &atoms[0],
                &atoms[1],
                &atoms[2],
                Atom::num(0),
                Symbol::PI.to_atom().pow(2) / 6,
            ),
            super::native_function_options(),
        )
        .expect("unique sheet dilogarithm");
    let variables = difference_parameters();
    let [first, second] = difference_sheets(&variables.map(Symbol::to_atom));
    function_map
        .add_function_with_options(
            difference_symbol(),
            variables.to_vec(),
            difference_body(&first, &second, first.log_over_one_minus()),
            super::native_function_options(),
        )
        .expect("unique sheet dilogarithm difference");
}

#[derive(Clone)]
pub(crate) struct SheetAtom {
    value: Atom,
    phase: Atom,
    odd: Atom,
}

impl SheetAtom {
    pub(crate) fn lower(value: Atom) -> Self {
        Self::from_value(value, false)
    }

    pub(crate) fn upper(value: Atom) -> Self {
        Self::from_value(value, true)
    }

    fn from_value(value: Atom, upper: bool) -> Self {
        Self::with_sign(value, if upper { Atom::num(1) } else { Atom::num(-1) })
    }

    pub(crate) fn with_real_axis_sign(value: Atom, upper: bool) -> Self {
        Self::from_value(value, upper)
    }

    pub(crate) fn scale_positive(&self, scale: &Atom) -> Self {
        Self {
            value: &self.value * scale,
            phase: self.phase.clone(),
            odd: self.odd.clone(),
        }
    }

    pub(crate) fn with_sign(value: Atom, zero_sign: Atom) -> Self {
        let n = negative(&re(&value));
        let imaginary = im(&value);
        let s = choose(&imaginary, sign_nonnegative(&imaginary), zero_sign);
        Self {
            value: (Atom::num(1) - Atom::num(2) * &n) * value,
            phase: &n * s,
            odd: n,
        }
    }

    pub(crate) fn value(&self) -> Atom {
        physical(&self.value, &self.odd)
    }

    pub(crate) fn log(&self) -> Atom {
        log_symbol().call((&self.value, &self.phase, &self.odd))
    }

    pub(crate) fn log_over_one_minus(&self) -> Atom {
        log_ratio_symbol().call((&self.value, &self.phase, &self.odd))
    }

    pub(crate) fn dilog(&self) -> Atom {
        dilog_symbol().call((&self.value, &self.phase, &self.odd))
    }
}

impl Mul for &SheetAtom {
    type Output = SheetAtom;
    fn mul(self, rhs: Self) -> Self::Output {
        product(self, rhs, false)
    }
}

impl Div for &SheetAtom {
    type Output = SheetAtom;
    fn div(self, rhs: Self) -> Self::Output {
        product(self, rhs, true)
    }
}

fn product_symbol(divide: bool, component: usize) -> Symbol {
    symbolica::symbol!(format!(
        "__olo_sheet_{}_{}",
        if divide { "quotient" } else { "product" },
        component
    ))
}

fn product(left: &SheetAtom, right: &SheetAtom, divide: bool) -> SheetAtom {
    let args = [
        left.value.clone(),
        left.phase.clone(),
        left.odd.clone(),
        right.value.clone(),
        right.phase.clone(),
        right.odd.clone(),
    ];
    SheetAtom {
        value: product_symbol(divide, 0).call(args.as_slice()),
        phase: product_symbol(divide, 1).call(args.as_slice()),
        odd: product_symbol(divide, 2).call(args.as_slice()),
    }
}

pub(crate) fn divided_difference(first: &SheetAtom, second: &SheetAtom) -> Atom {
    difference_symbol().call((
        &first.value,
        &first.phase,
        &first.odd,
        &second.value,
        &second.phase,
        &second.odd,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::prelude::{Complex, Float, RealLike};
    #[test]
    fn sheet_product_and_quotient_reference_values() {
        let args = ["sheet_pair_a", "sheet_pair_b", "sheet_pair_s"]
            .map(|s| symbolica::symbol!(s).to_atom());
        let x = SheetAtom::with_sign(args[0].clone(), args[2].clone());
        let y = SheetAtom::with_sign(args[1].clone(), -&args[2]);
        let mut map = FunctionMap::new();
        register(&mut map);
        let exact = [&x * &y, &x / &y].map(|q| {
            let values = [q.value.clone(), q.phase.clone(), q.log(), q.dilog()];
            Atom::evaluator_multiple(&values.iter().map(Atom::as_view).collect::<Vec<_>>(), &args)
                .function_map(map.clone().into())
                .direct_translation(true)
                .build()
                .unwrap()
        });
        let mut evaluators = exact
            .clone()
            .map(|e| e.map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64())));
        let mut precise = exact.map(|e| {
            e.map_coeff_with_prec(
                &|v| Complex::new(v.re.to_multi_prec_float(128), v.im.to_multi_prec_float(128)),
                128,
            )
        });
        let mut failures = vec![];
        let mut operation_counts = [0; 2];
        for (line, text) in include_str!("../tests/data/sheets.txt").lines().enumerate() {
            if text.starts_with('#') {
                continue;
            }
            let v = text
                .split_whitespace()
                .map(|s| s.parse::<f64>().unwrap())
                .collect::<Vec<_>>();
            assert_eq!(v.len(), 13, "invalid sheet fixture line {line}");
            assert!(v.iter().all(|x| x.is_finite()));
            assert!(
                v[5] == 0.0 || v[5] == 1.0,
                "invalid operation at line {line}"
            );
            operation_counts[v[5] as usize] += 1;
            let input = [
                Complex::new(v[0], v[1]),
                Complex::new(v[2], v[3]),
                Complex::new(v[4], 0.),
            ];
            let mut out = [Complex::new(0., 0.); 4];
            evaluators[v[5] as usize].evaluate(&input, &mut out);
            let precise_input =
                input.map(|z| Complex::new(Float::with_val(128, z.re), Float::with_val(128, z.im)));
            let mut precise_out =
                core::array::from_fn::<_, 4, _>(|_| Complex::new(Float::new(128), Float::new(128)));
            precise[v[5] as usize].evaluate(&precise_input, &mut precise_out);
            for (k, (a, b)) in out.iter().zip(precise_out).enumerate() {
                let error = (a.re - b.re.to_f64()).hypot(a.im - b.im.to_f64());
                assert!(
                    error.is_finite() && error < 2e-12,
                    "precision mismatch line={line} component={k}: {a:?} vs {b:?}"
                );
            }
            let expected = [
                Complex::new(v[6], v[7]),
                Complex::new(v[8], 0.),
                Complex::new(v[9], v[10]),
                Complex::new(v[11], v[12]),
            ];
            for (k, (a, b)) in out.into_iter().zip(expected).enumerate() {
                // The auxiliary Rust dilog oracle uses a truncated 24-term
                // series; its worst unit-circle error is about 8e-8. Raw values,
                // integer sheets and logs stay strict, as does the independent
                // Fortran integral suite. Native Li2 is separately checked at 128 bits above.
                let tolerance = if k == 3 { 2e-7 } else { 1e-11 };
                if !(a.re - b.re).hypot(a.im - b.im).is_finite()
                    || (a.re - b.re).hypot(a.im - b.im) > tolerance
                {
                    failures.push(format!(
                        "line={line} component={k} actual={a:?} expected={b:?}"
                    ));
                }
            }
        }
        assert_eq!(operation_counts, [162, 162], "incomplete sheet fixtures");
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    }
    #[test]
    fn continued_dilog_reference_values() {
        let x = symbolica::parse!("sheet_test_x");
        let phase = symbolica::parse!("sheet_test_phase");
        let odd = symbolica::parse!("sheet_test_odd");
        let mut map = FunctionMap::new();
        register(&mut map);
        let mut eval = dilog_body(
            &x,
            &phase,
            &odd,
            Atom::num(0),
            Symbol::PI.to_atom().pow(2) / 6,
        )
        .evaluator(&[x, phase, odd])
        .function_map(map.into())
        .build()
        .unwrap()
        .map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64()));
        for (r, im, p, er, ei) in [
            (
                1.4637302044942475,
                0.8724394147326465,
                0,
                -0.523198871257367,
                -0.689201652707516,
            ),
            (0.2, 0.3, -3_i32, 0.0018756431260822382, 2.316496618720058),
            (0.2, 0.3, -1, 1.5411219798572924, 0.9804786253650036),
            (0.2, 0.3, 1, 3.080368316588503, -0.3555393679900509),
            (2., 3., -1, -0.07881407482088121, 3.5866974706558468),
            (2., 0., -1, 2.3201804233130985, 3.4513922952232026),
            (0.2, 0., 2, 1.0747946000082484, 1.4020522830093165),
        ] {
            let out = eval.evaluate_single(&[
                Complex::new(r, im),
                Complex::new(p as f64, 0.),
                Complex::new(p.rem_euclid(2) as f64, 0.),
            ]);
            assert!(
                (out.re - er).hypot(out.im - ei) < 1e-12,
                "r={r} im={im} phase={p} actual={out:?} expected={er}+{ei}i"
            );
        }
    }
}
