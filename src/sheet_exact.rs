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
    // Bind the real argument once so the sign predicate uses the same value
    // and precision in its cancelling pair.
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

/// The continued `Li2(1-x)` on the sheet `(value, phase, odd)`.
///
/// With `complement`, the given `w = 1 - x` decides the inversion and, above
/// half, replaces `1 - x` in `log(1-x)` and `Li2(1-x)` (`-w/x` in the
/// inverted form): next to `x = 1` on a sheet with windings the dilogarithm
/// carries `i pi n log(1-x)`, and `1 - x` formed from a rounded `x` loses
/// `|1 - x|` relative digits ([complement_difference]).
fn dilog_body(
    value: &Atom,
    phase: &Atom,
    odd: &Atom,
    complement: Option<&Atom>,
    at_one: Atom,
    at_zero: Atom,
) -> Atom {
    let i_pi = i() * Symbol::PI.to_atom();
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    let even = phase - odd;
    let argument = physical(value, odd);
    let (nonunit, inverted) = match complement {
        None => (
            &argument - 1,
            negative(&re(&(Atom::num(1) - value * value.conj()))),
        ),
        // 1 - |1 - w|²: the side of the unit circle, which sheets with two or
        // more windings jump across, decided by w itself next to x = 1.
        Some(w) => (
            w.clone(),
            negative(&(Atom::num(2) * re(w) - re(&(w * w.conj())))),
        ),
    };
    let given = complement.map(|w| choose(&inverted, -w / &argument, w.clone()));
    let normalized = choose(&inverted, Atom::num(1) / &argument, argument);
    let below_half = negative(&(Atom::num(2) * re(&normalized) - 1));
    // Below half `1 - n` keeps its relative precision; far from one the given
    // complement (about -x) does not, so it is taken only above half.
    let one_minus = match given {
        None => Atom::num(1) - &normalized,
        Some(given) => choose(&below_half, Atom::num(1) - &normalized, given),
    };
    let logarithm = choose(&inverted, -value.log(), value.log());
    let even = choose(&inverted, -&even, even);
    let odd = choose(&inverted, -odd, odd.clone());
    let sheet_logarithm = &logarithm + &i_pi * (&even + &odd);
    let log_one_minus = one_minus.log();
    let below =
        &pi_squared / 6 - crate::dilog_atom(&normalized) - &log_one_minus * &sheet_logarithm;
    let above = crate::dilog_atom(&one_minus) - &log_one_minus * &i_pi * even;
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

fn log_given_symbol() -> Symbol {
    symbolica::symbol!("__olo_sheet_log_given")
}

fn difference_symbol() -> Symbol {
    symbolica::symbol!("__olo_sheet_dilog_difference")
}

fn complement_dilog_symbol() -> Symbol {
    symbolica::symbol!("__olo_sheet_dilog_complement")
}

fn complement_difference_symbol() -> Symbol {
    symbolica::symbol!("__olo_sheet_dilog_difference_complement")
}

fn complement_parameters() -> [Symbol; 2] {
    ["__olo_dd_w1", "__olo_dd_w2"].map(|s| symbolica::symbol!(s))
}

fn sheet_parameters() -> [Symbol; 3] {
    [
        symbolica::symbol!("__olo_sheet_value"),
        symbolica::symbol!("__olo_sheet_phase"),
        symbolica::symbol!("__olo_sheet_odd"),
    ]
}

/// `near` takes the logarithm of an unreflected value within 1/2 of one
/// through [log_one_plus_body] (as `log(v) + i pi phase`): the native and
/// SymJIT `log|v|` is only accurate in absolute terms, and `log(q)/(1-q)` of a
/// quotient of nearly coincident roots lost their relative gap in digits.
/// The Dual neighborhoods keep the plain logarithm, whose derivatives are
/// holomorphic.
fn log_ratio_body(arguments: &[Atom; 3], neighborhood: Atom, near: bool) -> Atom {
    let difference = Atom::num(1) - physical(&arguments[0], &arguments[2]);
    let mut logarithm = log_symbol().call(arguments.as_slice());
    if near {
        let [value, phase, odd] = arguments;
        let offset = value - Atom::num(1);
        let close = negative(&(Atom::num(4) * re(&(&offset * offset.conj())) - 1));
        let accurate =
            log_one_plus_symbol().call((value, &offset)) + i() * Symbol::PI.to_atom() * phase;
        logarithm = choose(odd, logarithm.clone(), choose(&close, accurate, logarithm));
    }
    let quotient = &logarithm / &difference;
    // A nonzero sheet logarithm at one is a pole, not a removable limit.
    let at_one = choose(&logarithm, quotient.clone(), neighborhood);
    choose(&difference, quotient, at_one)
}

/// The sheet logarithm `log(v) + i pi phase` given the offset `w = x - 1` of
/// the physical value, which the caller forms without rounding `x` (the
/// contour box's quotients of a pole and a root next to it,
/// [SheetAtom::log_given]): `near` takes an unreflected value within 1/2 of
/// one from `w` through [log_one_plus_body], as [log_ratio_body] does from
/// `v - 1`. The Dual neighborhoods keep the plain logarithm.
fn log_given_body(arguments: &[Atom; 3], w: &Atom, near: bool) -> Atom {
    let logarithm = log_symbol().call(arguments.as_slice());
    if !near {
        return logarithm;
    }
    let [value, phase, odd] = arguments;
    let close = negative(&(Atom::num(4) * re(&(w * w.conj())) - 1));
    let accurate = log_one_plus_symbol().call((value, w)) + i() * Symbol::PI.to_atom() * phase;
    choose(odd, logarithm.clone(), choose(&close, accurate, logarithm))
}

fn log_given_parameters() -> [Symbol; 4] {
    let [value, phase, odd] = sheet_parameters();
    [value, phase, odd, symbolica::symbol!("__olo_sheet_offset")]
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
    let mut map =
        map.with_replaced_body(log_ratio_symbol(), log_ratio_body(&arguments, local, false))?;
    let offset = log_given_parameters()[3].to_atom();
    map = map.with_replaced_body(log_given_symbol(), log_given_body(&arguments, &offset, false))?;
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

/// `near` selects [near_difference] where the direct difference cancels. The
/// Dual neighborhoods keep the direct difference: [log_one_plus_body] forms
/// real parts, which carry no holomorphic derivative.
///
/// With `complements` (`1 - x` of both arguments, [complement_difference]) the
/// dilogarithms take them, and the gap is `w2 - w1` when both arguments are
/// within 1/2 of one. Farther out `x1 - x2` loses nothing, and the
/// complements, formed by the caller from other expressions than the values,
/// need not round as the values whose dilogarithms the numerator holds.
fn difference_body(
    first: &SheetAtom,
    second: &SheetAtom,
    complements: Option<[&Atom; 2]>,
    coincident: Atom,
    near: bool,
) -> Atom {
    let difference = first.value() - second.value();
    let ([upper, lower], gap) = match complements {
        None => ([first.dilog(), second.dilog()], difference.clone()),
        Some([w1, w2]) => {
            let norm = |z: &Atom| re(&(z * z.conj()));
            let near = |w: &Atom| negative(&(Atom::num(4) * norm(w) - 1));
            let near_one = near(w1) * near(w2);
            (
                [first.dilog_complement(w1), second.dilog_complement(w2)],
                choose(&near_one, w2 - w1, difference.clone()),
            )
        }
    };
    let noncoincident = if near {
        // Symbolica's evaluator calls a registered function at every
        // occurrence of the call: the dilogarithms are passed to the
        // selection once, as its arguments.
        select_symbol().call((
            &upper,
            &lower,
            &gap,
            &first.value,
            &first.phase,
            &first.odd,
            &second.value,
            &second.phase,
            &second.odd,
        ))
    } else {
        (&upper - &lower) / &gap
    };
    choose(&difference, noncoincident, coincident)
}

fn select_symbol() -> Symbol {
    symbolica::symbol!("__olo_sheet_dilog_difference_select")
}

fn select_parameters() -> [Symbol; 3] {
    ["__olo_dds_upper", "__olo_dds_lower", "__olo_dds_gap"].map(|s| symbolica::symbol!(s))
}

/// The divided difference `(upper - lower)/gap` of the dilogarithms of
/// `first` and `second` ([difference_body]), or [near_difference] where it
/// cancels by more than 16.
///
/// The test takes `upper` and the numerator more than once. Symbolica's
/// evaluator calls a registered function anew at every occurrence of its
/// call, and with the dilogarithm calls inside the test the interpreted D0 and
/// C0 took twice as long at 700 bits; as arguments of this function they are
/// evaluated once.
fn select_body(first: &SheetAtom, second: &SheetAtom, upper: &Atom, lower: &Atom, gap: &Atom) -> Atom {
    let numerator = upper - lower;
    let direct = &numerator / gap;
    // The values cancel by more than 16: the direct difference loses more
    // than four bits. Only then are the other forms considered.
    let norm = |z: &Atom| re(&(z * z.conj()));
    let cancels = negative(&(Atom::num(256) * norm(&numerator) - norm(upper)));
    choose(&cancels, near_difference(first, second, direct.clone()), direct)
}

fn log_one_plus_symbol() -> Symbol {
    symbolica::symbol!("__olo_log_one_plus")
}

/// The principal `log(q)` of `q = 1 + w`, given both, with its relative
/// precision for small `|w|`.
///
/// `log(1 + w)` would take `log|1 + w|` of a rounded modulus, which next to
/// one has an absolute error of a few units: the native and SymJIT complex
/// logarithms form `|z|` within a few ulps, and their `log(q)/(1 - q)` of a
/// quotient `q = 1 + w` lost `|w|` relative digits. Here, for `|w| < 1/2`,
/// `log|1 + w| = log1p(2 Re w + |w|²)/2` with Kahan's
/// `log1p(t) = log(u) t/(u-1)`, `u = 1 + t` (a real logarithm), and the
/// argument is the angle of `1 + w`. Farther from one the logarithm of `q` is
/// direct (`q` itself is formed by the caller without rounding `1 + w`).
fn log_one_plus_body(q: &Atom, w: &Atom) -> Atom {
    let norm = re(&(w * w.conj()));
    let t = Atom::num(2) * re(w) + &norm;
    let u = super::grouped(&(Atom::num(1) + &t));
    let offset = &u - Atom::num(1);
    let log1p = choose(&offset, u.log() * &t / &offset, t.clone());
    // Im(1 + w) = Im w exactly, so the angle keeps its relative precision too.
    let angle = im(&(Atom::num(1) + w).log());
    let near = negative(&(Atom::num(4) * norm - 1));
    choose(&near, log1p / 2 + i() * angle, q.log())
}

fn series_symbol(degree: u32) -> Symbol {
    if degree == SERIES_DEGREES[0] {
        symbolica::symbol!("__olo_sheet_dilog_series_short")
    } else {
        symbolica::symbol!("__olo_sheet_dilog_series")
    }
}

/// Highest powers of the series variables kept in [series_body]: the short
/// expansion suffices in binary64 for `|s| < 1/4`, the long one for `|s| < 1.24`.
const SERIES_DEGREES: [u32; 2] = [10, 22];

/// `(f(s1) - f(s2))/(s1 - s2)` for `f(s) = Li2(1 - e^-s) = sum_n B_n s^(n+1)/(n+1)!`,
/// through `s^degree` (OneLOop's `sumterms_c`). `h_n = (s1^(n+1) -
/// s2^(n+1))/(s1 - s2)` is built from `h_n = s1 h_(n-1) + s2^n`, which has no
/// cancellation when `s1` and `s2` nearly coincide.
fn series_body(first: &Atom, second: &Atom, degree: u32) -> Atom {
    // a_n = B_n/n! from sum_(j<=n) a_j/(n+1-j)! = 0 (n >= 1), exactly.
    let factorial = |n: u32| Atom::num(symbolica::domains::integer::Integer::factorial(n));
    let mut bernoulli = vec![Atom::num(1)];
    for n in 1..=degree {
        let sum = (0..n).fold(Atom::num(0), |sum, j| {
            sum + &bernoulli[j as usize] / factorial(n + 1 - j)
        });
        bernoulli.push(-sum);
    }
    let mut power = Atom::num(1);
    let mut complete = Atom::num(1);
    let mut series = Atom::num(1);
    for n in 1..=degree {
        power = &power * second;
        complete = first * &complete + &power;
        if !bernoulli[n as usize].is_zero() {
            series += &bernoulli[n as usize] / Atom::num(n + 1) * &complete;
        }
    }
    series
}

/// The divided difference where `(Li2(1-x1) - Li2(1-x2))/(x1 - x2)` cancels.
///
/// The direct difference loses the digits that the two values share: about
/// `|x|/|x1 - x2|` for nearly coincident arguments (nearly coincident
/// triangle roots lost three digits in binary64), and `pi²/6` divided by
/// `|x1 - x2| |log x|` when both arguments are small, which includes a zero
/// partner. As OneLOop's `dilog2_c`, the representation of [dilog_body]
/// (reflected below `Re = 1/2`, inverted outside the unit circle) is divided
/// term by term instead: the polylogarithm through the Bernoulli expansion in
/// `s = -log(n)` or `-log(1-n)` of the normalized `n` ([series_body]), the
/// logarithms through `log(q)/(1-q)` of their quotients, formed from the
/// quotients' offsets from one ([log_one_plus_body]). No term cancels, and
/// `pi²/6` drops out.
///
/// [difference_body] calls it only where the two dilogarithms cancel by more
/// than 16 (four bits). The expansion is then used when both arguments share
/// their phase and parity; and their side of the unit circle, unless the
/// phase is -1, 0 or 1, where the continued `Li2(1-x)` is continuous across it
/// and the representation inside holds on both sides; and when the first
/// omitted term, bounded by `8 (max|s|²/39)^12 > 8 (max|s|/(2 pi))^24`, does
/// not change `1` at the working precision (the expansion through `s^10` is
/// used when `8 (max|s|²/39)^6` does not). In binary64 that holds for
/// `|s| < 1.24`, which covers nearly coincident arguments anywhere, but at 256
/// bits only for `|s| < 3e-3`: higher precisions mostly keep the direct
/// difference, whose loss stays visible in the tracked precision. Double-double
/// arithmetic represents `1 + bound` exactly and keeps it too. A zero partner
/// of an argument below 1/2 in modulus (`Li2(1) = pi²/6` on every sheet) has
/// the closed form `-Li2(x)/x - L(x) log(1-x)/x`.
fn near_difference(first: &SheetAtom, second: &SheetAtom, direct: Atom) -> Atom {
    let one = Atom::num(1);
    let i_pi = i() * Symbol::PI.to_atom();
    let norm = |z: &Atom| re(&(z * z.conj()));
    let less = |a: Atom, b: &Atom| negative(&(a - b));
    let log_one_plus = |q: &Atom, w: &Atom| log_one_plus_symbol().call((q, w));
    // log(q)/(1-q) of a quotient q = 1 + w on its principal sheet.
    let ratio_log = |q: &Atom, w: &Atom| choose(w, log_one_plus(q, w) / -w, Atom::num(-1));
    // The sheet logarithm L(x) = log(v) + i pi phase.
    let logarithm = |s: &SheetAtom| s.value.log() + &i_pi * &s.phase;
    let [x1, x2] = [first.value(), second.value()];
    let [m1, m2] = [norm(&x1), norm(&x2)];
    let small = |m: &Atom| less(Atom::num(4) * m, &one);
    let zero_partner = |x: &Atom, s: &SheetAtom| {
        -(crate::dilog_atom(x) / x) - log_one_plus(&(&one - x), &-x) / x * logarithm(s)
    };

    // dilog_body's inversion test, on the same arithmetic.
    let inverted = |s: &SheetAtom| negative(&re(&(&one - &s.value * s.value.conj())));
    let [inverted_1, inverted_2] = [inverted(first), inverted(second)];
    // Across the unit circle the continued Li2(1-x) is continuous for the
    // phases -1, 0, 1, and the representation inside holds on both sides;
    // higher windings jump there (OneLOop's dilog2_c adds the jump).
    let straddle = &inverted_1 - &inverted_2;
    let phase = &second.phase;
    let unwound = choose(&(phase * (phase * phase - 1)), Atom::num(0), Atom::num(1));
    let inversion = choose(&straddle, Atom::num(0), inverted_2.clone());
    let normalized = |x: &Atom| choose(&inversion, &one / x, x.clone());
    let [n1, n2] = [normalized(&x1), normalized(&x2)];
    // 1 - n, formed without cancellation next to n = 1.
    let complement = |x: &Atom| choose(&inversion, (x - &one) / x, &one - x);
    let [c1, c2] = [complement(&x1), complement(&x2)];
    // n1 - n2, without forming the rounded reciprocals' difference.
    let gap_n = choose(&inversion, (&x2 - &x1) * &n1 * &n2, &x1 - &x2);
    // 1/n2, without another reciprocal outside the unit circle.
    let inverse_n2 = choose(&inversion, x2.clone(), &one / &x2);
    let below = negative(&(Atom::num(2) * re(&n2) - 1));
    let sign = &one - Atom::num(2) * &inversion;
    // -log(1-n); above 1/2 the parity is even and -log(n) = -sign log(v).
    let variable = |n: &Atom, c: &Atom, s: &SheetAtom| {
        choose(&below, -log_one_plus(c, &-n), -(&sign * s.value.log()))
    };
    let [s1, s2] = [variable(&n1, &c1, first), variable(&n2, &c2, second)];
    let [r1, r2] = [norm(&s1), norm(&s2)];
    let largest = choose(&less(r1.clone(), &r2), r2, r1);
    // The first omitted terms of the short and long expansions.
    let truncated = |power: i64| {
        super::grouped(&(&one + Atom::num(8) * (&largest / 39).pow(power))) - &one
    };
    let series = choose(
        &truncated(6),
        series_symbol(SERIES_DEGREES[1]).call((&s1, &s2)),
        series_symbol(SERIES_DEGREES[0]).call((&s1, &s2)),
    );
    let truncated = truncated(12);
    // log(q)/(1-q) for q = n1/n2 and (1-n1)/(1-n2), from their offsets from one.
    let rho = ratio_log(&(&n1 * &inverse_n2), &(&gap_n * &inverse_n2)) * &inverse_n2;
    let sigma = ratio_log(&(&c1 / &c2), &(-&gap_n / &c2)) / &c2;
    let [l1, l2] = [logarithm(first), logarithm(second)];
    let even = &sign * (&second.phase - &second.odd);
    // pi²/6 - Li2(n) - log(1-n) Lambda(n), with Lambda = sign L.
    let reflected = &sigma * (&series - &sign * &l2) - &s1 * &rho;
    // Li2(1-n) - i pi even log(1-n).
    let unreflected = &series * &rho - choose(&even, &i_pi * &even * &sigma, Atom::num(0));
    let normalized_difference = choose(&below, reflected, unreflected);
    // Outside: Li2(1-x) = -B(1/x) - L²/2 with B the normalized form above.
    let offset = (&first.value - &second.value) / &second.value;
    let outer = &n1 * &n2 * &normalized_difference
        + (&l1 + &l2) / 2 * ratio_log(&(&first.value / &second.value), &offset) * &n2;
    let expansion = choose(&inversion, outer, normalized_difference);

    // The form to use: 0 the direct difference, 1 the expansion, 2 and 3 the
    // zero-partner limits of x1 and x2. The direct difference appears once.
    let [direct_kind, expansion_kind] = [Atom::num(0), Atom::num(1)];
    let candidate = choose(
        &(&first.phase - &second.phase),
        direct_kind.clone(),
        choose(
            &(&first.odd - &second.odd),
            direct_kind.clone(),
            choose(
                &straddle,
                choose(
                    &unwound,
                    choose(&truncated, direct_kind.clone(), expansion_kind.clone()),
                    direct_kind.clone(),
                ),
                choose(&truncated, direct_kind.clone(), expansion_kind),
            ),
        ),
    );
    let kind = choose(
        &second.value,
        choose(
            &first.value,
            candidate,
            choose(&small(&m2), Atom::num(3), direct_kind.clone()),
        ),
        choose(&small(&m1), Atom::num(2), direct_kind),
    );
    choose(
        &kind,
        choose(
            &(&kind - 1),
            choose(&(&kind - 2), zero_partner(&x2, second), zero_partner(&x1, first)),
            expansion,
        ),
        direct,
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
            None,
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
    map = map.with_replaced_body(
        difference_symbol(),
        difference_body(&first, &second, None, coincident, false),
    )?;
    // The complements only restore digits of 1 - x: the neighborhoods take
    // the plain difference, whose derivatives they define.
    map.with_replaced_body(
        complement_difference_symbol(),
        divided_difference(&first, &second),
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
    // This helper is used on real components or inside real-axis branches.
    // Bind x once, and test x+|x|: x-|x| multiplies |x| by a prepared -1
    // coefficient, which can round a stronger tracked Float input and leave
    // a nonzero residue even for positive x. The sum cancels a negative x
    // against its own norm without that multiplication. Zero keeps sign +1.
    let sign_value = if grouped_predicates {
        super::grouped(&x)
    } else {
        x.clone()
    };
    function_map
        .add_function_with_options(
            symbolica::symbol!("__olo_sign_nonnegative"),
            vec![sign_arg],
            choose(
                &(&sign_value + Symbol::ABS.call((&sign_value,))),
                Atom::num(1),
                choose(&sign_value, Atom::num(-1), Atom::num(1)),
            ),
            super::native_function_options(),
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
            log_ratio_body(&atoms, Atom::num(-1), std::env::var("ZZ_GEN_LOG").is_err()),
            super::native_function_options(),
        )
        .expect("unique sheet logarithm ratio");
    let given = log_given_parameters();
    function_map
        .add_function_with_options(
            log_given_symbol(),
            given.to_vec(),
            log_given_body(&atoms, &given[3].to_atom(), true),
            super::native_function_options(),
        )
        .expect("unique sheet logarithm of a given offset");
    function_map
        .add_function_with_options(
            dilog_symbol(),
            variables.to_vec(),
            dilog_body(
                &atoms[0],
                &atoms[1],
                &atoms[2],
                None,
                Atom::num(0),
                Symbol::PI.to_atom().pow(2) / 6,
            ),
            super::native_function_options(),
        )
        .expect("unique sheet dilogarithm");
    let complement = symbolica::symbol!("__olo_sheet_complement");
    function_map
        .add_function_with_options(
            complement_dilog_symbol(),
            [variables.as_slice(), &[complement]].concat(),
            dilog_body(
                &atoms[0],
                &atoms[1],
                &atoms[2],
                Some(&complement.to_atom()),
                Atom::num(0),
                Symbol::PI.to_atom().pow(2) / 6,
            ),
            super::native_function_options(),
        )
        .expect("unique sheet dilogarithm of a complement");
    let [q, w] = [
        symbolica::symbol!("__olo_log_one_plus_value"),
        symbolica::symbol!("__olo_log_one_plus_offset"),
    ];
    function_map
        .add_function_with_options(
            log_one_plus_symbol(),
            vec![q, w],
            log_one_plus_body(&q.to_atom(), &w.to_atom()),
            super::native_function_options(),
        )
        .expect("unique log(1 + w)");
    let [a, b] = [
        symbolica::symbol!("__olo_dilog_series_a"),
        symbolica::symbol!("__olo_dilog_series_b"),
    ];
    for degree in SERIES_DEGREES {
        function_map
            .add_function_with_options(
                series_symbol(degree),
                vec![a, b],
                series_body(&a.to_atom(), &b.to_atom(), degree),
                super::native_function_options(),
            )
            .expect("unique dilogarithm series");
    }
    let variables = difference_parameters();
    let [first, second] = difference_sheets(&variables.map(Symbol::to_atom));
    // At coincidence the plain log(x)/(1-x), the primal of the Dual
    // neighborhoods' limit (`with_dilog_neighborhoods`).
    let coincident = log_ratio_body(
        &[first.value.clone(), first.phase.clone(), first.odd.clone()],
        Atom::num(-1),
        false,
    );
    function_map
        .add_function_with_options(
            difference_symbol(),
            variables.to_vec(),
            difference_body(
                &first,
                &second,
                None,
                coincident.clone(),
                std::env::var("ZZ_GEN_DD").is_err(),
            ),
            super::native_function_options(),
        )
        .expect("unique sheet dilogarithm difference");
    let selection = select_parameters();
    let [upper, lower, gap] = selection.map(Symbol::to_atom);
    function_map
        .add_function_with_options(
            select_symbol(),
            [selection.as_slice(), variables.as_slice()].concat(),
            select_body(&first, &second, &upper, &lower, &gap),
            super::native_function_options(),
        )
        .expect("unique sheet dilogarithm difference selection");
    let complements = complement_parameters();
    let [w1, w2] = complements.map(Symbol::to_atom);
    function_map
        .add_function_with_options(
            complement_difference_symbol(),
            [variables.as_slice(), complements.as_slice()].concat(),
            difference_body(&first, &second, Some([&w1, &w2]), coincident, true),
            super::native_function_options(),
        )
        .expect("unique sheet dilogarithm difference of complements");
}

#[derive(Clone)]
pub(crate) struct SheetAtom {
    value: Atom,
    phase: Atom,
    odd: Atom,
}

impl SheetAtom {
    #[cfg(test)]
    pub(crate) fn raw(value: Atom, phase: Atom, odd: Atom) -> Self {
        Self { value, phase, odd }
    }

    pub(crate) fn lower(value: Atom) -> Self {
        Self::from_value(value, false)
    }

    pub(crate) fn upper(value: Atom) -> Self {
        Self::from_value(value, true)
    }

    fn from_value(value: Atom, upper: bool) -> Self {
        Self::with_sign(value, if upper { Atom::num(1) } else { Atom::num(-1) })
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

    /// [Self::log] given `w = x - 1`, formed by the caller without rounding
    /// the physical value `x` ([log_given_body]).
    pub(crate) fn log_given(&self, w: &Atom) -> Atom {
        log_given_symbol().call((&self.value, &self.phase, &self.odd, w))
    }

    /// This sheet without its winding, and the winding `n` (in units of
    /// `i pi`): `log(self) = log(unwound) + i pi n`. Only an unreflected value
    /// (`odd = 0`, the case of a quotient near one) is unwound.
    pub(crate) fn split_winding(&self) -> (SheetAtom, Atom) {
        let winding = (1 - &self.odd) * &self.phase;
        (
            SheetAtom {
                value: self.value.clone(),
                phase: &self.phase - &winding,
                odd: self.odd.clone(),
            },
            winding,
        )
    }

    /// This value `turns` sheets further: `log(self) + 2 pi i turns`.
    pub(crate) fn wound(&self, turns: &Atom) -> Self {
        Self {
            value: self.value.clone(),
            phase: &self.phase + Atom::num(2) * turns,
            odd: self.odd.clone(),
        }
    }

    pub(crate) fn dilog(&self) -> Atom {
        dilog_symbol().call((&self.value, &self.phase, &self.odd))
    }

    /// [SheetAtom::dilog] given the complement `1 - x` of the value.
    pub(crate) fn dilog_complement(&self, complement: &Atom) -> Atom {
        complement_dilog_symbol().call((&self.value, &self.phase, &self.odd, complement))
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

/// [divided_difference] given the complements `w1 = 1 - x1` and
/// `w2 = 1 - x2` of the two arguments, for arguments that can lie next to one.
///
/// On a sheet with `n` windings the dilogarithm carries `-i pi n log(1-x)`,
/// which next to `x = 1` is the logarithm of a small quantity: of an invariant
/// between two massless lines, for the box arguments of a collinear limit.
/// Taken from a rounded `x`, `1 - x` loses `|1 - x|` relative digits and the
/// dilogarithm `eps/|1 - x|` absolute ones (a one-mass box with p1² = -1e-10
/// kept 2e-5 in binary64). The caller forms the complements from the
/// kinematics; the dilogarithms take them ([dilog_body]), and the gap is
/// `w2 - w1` where it is the smaller form. Cancelling differences take
/// [near_difference] from the values, as [divided_difference] does.
pub(crate) fn complement_difference(
    first: &SheetAtom,
    first_complement: &Atom,
    second: &SheetAtom,
    second_complement: &Atom,
) -> Atom {
    complement_difference_symbol().call((
        &first.value,
        &first.phase,
        &first.odd,
        &second.value,
        &second.phase,
        &second.odd,
        first_complement,
        second_complement,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::prelude::{Complex, Float, RealLike};

    /// [near_difference] against the direct difference at 400 bits, on pairs
    /// of one sheet (phases -3..3) where the direct difference cancels: both
    /// arguments small, nearly coincident (relative gaps 1e-9 to 0.3, also
    /// across the unit circle for phases -1, 0, 1), and a zero partner; and
    /// both large, where the direct difference is used. Binary64 must be
    /// accurate to about a hundred units (the direct difference of gaps above
    /// 1/32 on winding sheets loses up to seven bits), and 256 bits too (the direct
    /// difference, whose cancellation is at most 2^50 here, or the expansion
    /// where its bound allows it). The direct difference alone fails this by up
    /// to eight orders of magnitude.
    #[test]
    fn near_divided_differences_match_high_precision() {
        let s = ["dd_v1", "dd_p1", "dd_o1", "dd_v2", "dd_p2", "dd_o2"]
            .map(|name| symbolica::symbol!(name).to_atom());
        let first = SheetAtom::raw(s[0].clone(), s[1].clone(), s[2].clone());
        let second = SheetAtom::raw(s[3].clone(), s[4].clone(), s[5].clone());
        let expressions = [
            divided_difference(&first, &second),
            (first.dilog() - second.dilog()) / (first.value() - second.value()),
        ];
        let mut map = FunctionMap::new();
        register(&mut map);
        let exact = Atom::evaluator_multiple(
            &expressions.iter().map(Atom::as_view).collect::<Vec<_>>(),
            &s,
        )
        .function_map(map.into())
        .build()
        .unwrap();
        let mut binary64 = exact
            .clone()
            .map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64()));
        let precise = |bits: u32| {
            exact.clone().map_coeff_with_prec(
                &|v| {
                    Complex::new(
                        v.re.to_multi_prec_float(bits),
                        v.im.to_multi_prec_float(bits),
                    )
                },
                bits,
            )
        };
        let (mut p256, mut p400) = (precise(256), precise(400));
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut uniform = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let polar = |modulus: f64, uniform: &mut dyn FnMut() -> f64| {
            let angle = (uniform() - 0.5) * std::f64::consts::PI;
            Complex::new(modulus * angle.cos(), modulus * angle.sin())
        };
        let (mut worst, mut worst_256) = (0f64, 0f64);
        let mut failures = Vec::new();
        for k in 0..4000 {
            let phase = (uniform() * 7.0) as i64 - 3;
            let (v1, v2) = match k % 4 {
                0 => (
                    polar(10f64.powf(-6.0 + 5.09 * uniform()), &mut uniform),
                    polar(10f64.powf(-6.0 + 5.09 * uniform()), &mut uniform),
                ),
                1 => (
                    polar(10f64.powf(0.35 + 4.5 * uniform()), &mut uniform),
                    polar(10f64.powf(0.35 + 4.5 * uniform()), &mut uniform),
                ),
                2 => {
                    let mut modulus = 10f64.powf(-1.3 + 2.6 * uniform());
                    if (modulus - 1.0).abs() < 0.02 {
                        modulus = 1.5;
                    }
                    let v2 = polar(modulus, &mut uniform);
                    let gap = 10f64.powf(-9.0 + 8.5 * uniform());
                    let angle = 2.0 * std::f64::consts::PI * uniform();
                    let factor = Complex::new(1.0 + gap * angle.cos(), gap * angle.sin());
                    let v1 = Complex::new(
                        v2.re * factor.re - v2.im * factor.im,
                        v2.re * factor.im + v2.im * factor.re,
                    );
                    let crosses = (v1.norm_squared() - 1.0) * (v2.norm_squared() - 1.0) <= 0.0;
                    if v1.re < 0.0 || (crosses && phase.abs() > 1) {
                        continue;
                    }
                    (v1, v2)
                }
                _ => (
                    polar(10f64.powf(-6.0 + 5.65 * uniform()), &mut uniform),
                    Complex::new(0.0, 0.0),
                ),
            };
            let odd = phase.rem_euclid(2) as f64;
            let input = [
                v1,
                Complex::new(phase as f64, 0.0),
                Complex::new(odd, 0.0),
                v2,
                Complex::new(phase as f64, 0.0),
                Complex::new(odd, 0.0),
            ];
            let mut out = [Complex::new(0.0, 0.0); 2];
            binary64.evaluate(&input, &mut out);
            let exact_input = |bits: u32| {
                input.map(|z| Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im)))
            };
            let fresh = |bits: u32| {
                core::array::from_fn::<_, 2, _>(|_| Complex::new(Float::new(bits), Float::new(bits)))
            };
            let (mut o256, mut o400) = (fresh(256), fresh(400));
            p256.evaluate(&exact_input(256), &mut o256);
            p400.evaluate(&exact_input(400), &mut o400);
            let reference = &o400[1];
            let scale = reference.re.to_f64().hypot(reference.im.to_f64());
            let error = (out[0].re - reference.re.to_f64()).hypot(out[0].im - reference.im.to_f64()) / scale;
            let error_256 = (o256[0].re.clone() - &reference.re)
                .to_f64()
                .hypot((o256[0].im.clone() - &reference.im).to_f64())
                / scale;
            worst = worst.max(error);
            worst_256 = worst_256.max(error_256);
            if !(error <= 3e-14) || !(error_256 <= 1e-60) {
                failures.push(format!("{input:?}: binary64 {error:e}, 256 bits {error_256:e}"));
            }
        }
        eprintln!("worst binary64 {worst:e}, worst 256 bits {worst_256:e}");
        assert!(failures.is_empty(), "{}", failures[..failures.len().min(10)].join("\n"));
    }

    /// [complement_difference] next to one: `x1 = 1 - w1` with `|w1|` from
    /// 1e-14 to 1e-2 on sheets with windings (phases -4..4, even), `x2 = 1 - w2`
    /// anywhere else, both rounded to binary64 and passed with their exact
    /// complements. Binary64 must be within 1e-14 of the plain difference of
    /// `1 - w1` and `1 - w2` at 400 bits (it is within 1e-15; the plain
    /// difference of the rounded arguments was off by up to 3e-4).
    #[test]
    fn complement_differences_next_to_one() {
        let s = ["cd_v1", "cd_p1", "cd_o1", "cd_v2", "cd_p2", "cd_o2", "cd_w1", "cd_w2"]
            .map(|name| symbolica::symbol!(name).to_atom());
        let first = SheetAtom::raw(s[0].clone(), s[1].clone(), s[2].clone());
        let second = SheetAtom::raw(s[3].clone(), s[4].clone(), s[5].clone());
        let expressions = [
            complement_difference(&first, &s[6], &second, &s[7]),
            divided_difference(&first, &second),
        ];
        let mut map = FunctionMap::new();
        register(&mut map);
        let exact = Atom::evaluator_multiple(
            &expressions.iter().map(Atom::as_view).collect::<Vec<_>>(),
            &s,
        )
        .function_map(map.into())
        .build()
        .unwrap();
        let mut binary64 = exact
            .clone()
            .map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64()));
        let mut p400 = exact.map_coeff_with_prec(
            &|v| Complex::new(v.re.to_multi_prec_float(400), v.im.to_multi_prec_float(400)),
            400,
        );
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut uniform = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let one = Complex::new(1.0, 0.0);
        // (stored value, parity) of x on a sheet with the given phase parity.
        let stored = |x: Complex<f64>| if x.re < 0.0 { (-x, 1.0) } else { (x, 0.0) };
        let (mut worst, mut worst_plain) = (0f64, 0f64);
        let mut failures = Vec::new();
        for _ in 0..2000 {
            let phase = 2.0 * ((uniform() * 5.0) as i64 - 2) as f64;
            let angle = 2.0 * std::f64::consts::PI * uniform();
            let modulus = 10f64.powf(-14.0 + 12.0 * uniform());
            let w1 = Complex::new(modulus * angle.cos(), modulus * angle.sin());
            let angle = 2.0 * std::f64::consts::PI * uniform();
            let modulus = 10f64.powf(-0.5 + 1.2 * uniform());
            let w2 = Complex::new(modulus * angle.cos(), modulus * angle.sin());
            let (v1, o1) = stored(one - w1);
            let (v2, o2) = stored(one - w2);
            let p2 = phase + o2;
            let input = [
                v1,
                Complex::new(phase, 0.0),
                Complex::new(o1, 0.0),
                v2,
                Complex::new(p2, 0.0),
                Complex::new(o2, 0.0),
                w1,
                w2,
            ];
            let mut out = [Complex::new(0.0, 0.0); 2];
            binary64.evaluate(&input, &mut out);
            // The reference takes the exact arguments 1 - w.
            let float = |z: Complex<f64>| {
                Complex::new(Float::with_val(400, z.re), Float::with_val(400, z.im))
            };
            let exact_value = |w: Complex<f64>, odd: f64| {
                let x = Complex::new(Float::with_val(400, 1.0), Float::with_val(400, 0.0))
                    - float(w);
                if odd == 0.0 { x } else { -x }
            };
            let high = [
                exact_value(w1, o1),
                float(input[1]),
                float(input[2]),
                exact_value(w2, o2),
                float(input[4]),
                float(input[5]),
                float(w1),
                float(w2),
            ];
            let mut o400 = core::array::from_fn::<_, 2, _>(|_| {
                Complex::new(Float::new(400), Float::new(400))
            });
            p400.evaluate(&high, &mut o400);
            let reference = Complex::new(o400[1].re.to_f64(), o400[1].im.to_f64());
            let error = |value: Complex<f64>| {
                (value.re - reference.re).hypot(value.im - reference.im)
                    / reference.re.hypot(reference.im)
            };
            worst = worst.max(error(out[0]));
            worst_plain = worst_plain.max(error(out[1]));
            if !(error(out[0]) <= 1e-14) {
                failures.push(format!("{input:?}: {:e}", error(out[0])));
            }
        }
        eprintln!("worst binary64 {worst:e} (plain difference {worst_plain:e})");
        assert!(failures.is_empty(), "{}", failures[..failures.len().min(10)].join("\n"));
    }

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
            None,
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
