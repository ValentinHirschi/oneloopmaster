//! Triangle sectors and their native expressions.
use super::*;
use sheet_exact::sqrt_lower as physical_sqrt;

pub(crate) fn negative_quadratic_roots(a: &Atom, b: &Atom, c: &Atom) -> (Atom, Atom) {
    let discriminant = physical_sqrt(&(b * b - Atom::num(4) * a * c));
    (
        (b + &discriminant) / (Atom::num(2) * a),
        (b - discriminant) / (Atom::num(2) * a),
    )
}

/// Labelled roots `(-b + d)/(2a)` and `(-b - d)/(2a)` of `a x² + b x + c`,
/// where `d² = b² - 4ac`.
///
/// The root of larger magnitude is evaluated directly and the other one from
/// Vieta's product `c/a`, as in OneLOop's `solabc`; the direct formula loses
/// all digits of a small root in `-b ∓ d`. The labels still follow the sign
/// in front of `d`, so callers keep their analytic lips. Requires `a c ≠ 0`.
pub(crate) fn labelled_quadratic_roots(a: &Atom, b: &Atom, c: &Atom, d: &Atom) -> (Atom, Atom) {
    let plus = d - b;
    let minus = -b - d;
    let size_difference = grouped(&(Symbol::ABS.call((&plus,)) - Symbol::ABS.call((&minus,))));
    let minus_larger = &size_difference - Symbol::ABS.call((&size_difference,));
    let two_a = Atom::num(2) * a;
    let two_c = Atom::num(2) * c;
    (
        if_nonzero_else(&minus_larger, &two_c / &minus, &plus / &two_a),
        if_nonzero_else(&minus_larger, &minus / &two_a, &two_c / &plus),
    )
}

/// Roots `X1, X2` of `a X² - b X + c`, the negated roots of OneLOop's
/// `solabc(a, b, c)`, in its order: `X1` has the larger magnitude and is
/// evaluated directly, `X2 = c/(a X1)` from Vieta's product. `X2` is then
/// exactly zero when `c = 0`, which the finite formulas require (they divide
/// by `X2` only in terms that vanish with it), and a small root keeps its
/// digits. The finite formulas are symmetric in the two labels otherwise.
///
/// At `c = 0`, `X1` is Vieta's sum `b/a - X2`, OneLOop's `b/a` in value:
/// the direct formula would add `b` to a square root that equals `±b` only
/// up to rounding. The difference keeps `X1`'s dependence on `c`.
///
/// `discriminant` is a square root of `b² - 4ac`; callers pass it from an
/// exact kinematic identity so that coincident roots stay exactly equal.
pub(crate) fn ordered_negative_quadratic_roots(
    a: &Atom,
    b: &Atom,
    c: &Atom,
    discriminant: &Atom,
) -> (Atom, Atom) {
    let plus = b + discriminant;
    let minus = b - discriminant;
    let size_difference = grouped(&(Symbol::ABS.call((&plus,)) - Symbol::ABS.call((&minus,))));
    let minus_larger = &size_difference - Symbol::ABS.call((&size_difference,));
    let larger = if_nonzero_else(&minus_larger, minus, plus);
    let small = if_nonzero(&larger, Atom::num(2) * c / &larger);
    let first = if_nonzero_else(c, &larger / (Atom::num(2) * a), b / a - &small);
    let second = if_nonzero_else(discriminant, small, first.clone());
    (first, second)
}

/// Sign selecting the real-axis lip of a root `x` of `a x² - b x + c`.
///
/// OneLOop shifts external invariants as `p -> p (1 + i eps)` and lets the
/// roots of its finite formulas inherit the imaginary parts. `velocity` holds
/// `D(a), D(b), D(c)`, where `dq/d eps = i D(q)`. To first order
/// `Im x = eps D(x)` with `D(x) = (D(b) x - D(a) x² - D(c)) / (2 a x - b)`.
/// An exactly real root is placed on the upper lip when `D(x)` vanishes, as
/// by OneLOop's `qonv(x, 1)`. Only real-axis roots use this sign.
pub(crate) fn root_lip(x: &Atom, a: &Atom, b: &Atom, velocity: [&Atom; 3]) -> Atom {
    root_lip_or(x, a, b, velocity, Atom::num(1))
}

/// [root_lip] with an explicit sign for a vanishing first-order shift.
pub(crate) fn root_lip_or(
    x: &Atom,
    a: &Atom,
    b: &Atom,
    velocity: [&Atom; 3],
    fallback: Atom,
) -> Atom {
    let [da, db, dc] = velocity;
    let numerator = db * x - da * x * x - dc;
    let slope = Atom::num(2) * a * x - b;
    let shift = grouped(&sheet_exact::re(&(numerator * slope.conj())));
    if_nonzero_else(&shift, sheet_exact::sign_nonnegative(&shift), fallback)
}

/// Velocities `D(r), D(d)` of `r_function(k)` for a velocity `D(k)`.
///
/// Away from `k = ±2`, `dr = r dk/d` and `dd = k dk/d` with `d = r - 1/r`.
/// At a threshold `r = ±1` the root instead moves like `sqrt(dk)`, along
/// `r - i` for both signs (OneLOop's `rfun` keeps `|r| >= 1`), which dominates
/// every first-order velocity: in units of `sqrt(eps)`, `D(r) = -sqrt|D(k)|`
/// and `D(d) = 2 D(r)`. The returned weight is zero when any of the
/// r-functions sits at such a threshold; it multiplies the first-order
/// velocities, and the regular r-function velocities are dropped as well.
pub(crate) fn r_function_velocities<const N: usize>(
    functions: [[&Atom; 4]; N],
) -> ([[Atom; 2]; N], Atom) {
    let singular = functions.map(|[_, _, d, dk]| {
        if_nonzero_else(
            d,
            Atom::num(0),
            if_nonzero_else(dk, Atom::num(1), Atom::num(0)),
        )
    });
    let any_singular = any_nonzero(singular.iter().cloned());
    let velocities = core::array::from_fn(|index| {
        let [k, r, d, dk] = functions[index];
        let step = -Symbol::ABS.call((dk,)).sqrt();
        // A vanishing D(k) gives no shift, even where d vanishes as well.
        let regular = [r * dk / d, k * dk / d].map(|v| if_nonzero(d, v));
        let threshold = [
            &singular[index] * &step,
            Atom::num(2) * &singular[index] * &step,
        ];
        let [regular_r, regular_d] = regular;
        let [threshold_r, threshold_d] = threshold;
        [
            if_nonzero_else(&any_singular, threshold_r, regular_r),
            if_nonzero_else(&any_singular, threshold_d, regular_d),
        ]
    });
    (
        velocities,
        if_nonzero_else(&any_singular, Atom::num(0), Atom::num(1)),
    )
}

/// [r_function_velocities] for a single r-function: `[D(r), D(d), weight]`.
pub(crate) fn r_function_velocity(k: &Atom, r: &Atom, d: &Atom, dk: &Atom) -> [Atom; 3] {
    let ([[dr, dd]], weight) = r_function_velocities([[k, r, d, dk]]);
    [dr, dd, weight]
}

/// [r_function] of `k = (m_i + m_j - p)/(s_i s_j)`, with `s = sqrt(m)`.
///
/// Exactly known roots are returned with their exact values, decided by exact
/// tests on the inputs rather than on the rounded `k`: near `|k| = 2` a
/// rounding error in `k` becomes a square-root sized error in `r`, can move a
/// threshold root off the real axis, and turns exact zeros of later
/// coefficients into noise. No tolerance is involved.
///
/// - At `p = 0` with unequal masses the root is `E + δ` from
///   [zero_momentum_root]: exactly `E = s_j/s_i` or `s_i/s_j` in value, with
///   the analytic dependence of the root on every input, including `p`.
/// - At `p = 0` with equal masses the root is 1, and at a (pseudo-)threshold
///   `(p - m_i - m_j)² = 4 m_i m_j`, where `k = ±2`, it is `±1`; both with
///   `r - 1/r = 0`. These are square-root branch points, where the root has
///   no finite derivatives in `p` anyway.
pub(crate) fn r_function_of(
    momentum: &Atom,
    mass_i: &Atom,
    mass_j: &Atom,
    sqrt_i: &Atom,
    sqrt_j: &Atom,
) -> (Atom, Atom) {
    let k = (mass_i + mass_j - momentum) / (sqrt_i * sqrt_j);
    let (root, difference) = r_function(&k);
    let unequal = mass_i - mass_j;
    let [nearer_ratio, exact, delta] = zero_momentum_root(momentum, mass_i, mass_j, sqrt_i, sqrt_j);
    let ratio = sqrt_j / sqrt_i;
    let inverse = sqrt_i / sqrt_j;
    let exact_difference = if_nonzero_else(&nearer_ratio, &ratio - &inverse, inverse - ratio);
    let shifted = &exact + &delta;
    // r - 1/r = (E - 1/E) + δ (1 + 1/(E (E + δ))), exactly E - 1/E at p = 0.
    let shifted_difference =
        exact_difference + &delta * (Atom::num(1) + Atom::num(1) / (&exact * &shifted));
    // Equal masses test p - 4 m directly: the Källén form would round 3 m.
    let off_threshold = if_nonzero_else(
        &unequal,
        (momentum - mass_i - mass_j).pow(2) - Atom::num(4) * mass_i * mass_j,
        momentum - Atom::num(4) * mass_i,
    );
    let threshold_root = sheet_exact::sign_nonnegative(&sheet_exact::re(&k));
    (
        if_nonzero_else(
            momentum,
            if_nonzero_else(&off_threshold, root, threshold_root),
            if_nonzero_else(&unequal, shifted, Atom::num(1)),
        ),
        if_nonzero_else(
            momentum,
            if_nonzero_else(&off_threshold, difference, Atom::num(0)),
            if_nonzero_else(&unequal, shifted_difference, Atom::num(0)),
        ),
    )
}

/// `s_i/s_j - 1/r` for the root `r` of [r_function_of].
///
/// It vanishes identically at `p = 0` when `r = s_j/s_i`, and is then formed
/// as `δ/(E (E + δ))` with the [zero_momentum_root] data, which is exactly
/// zero there while keeping the analytic dependence on every input. The
/// direct difference would be rounding noise, and a later quadratic root
/// that must vanish (OneLOop's `c = 0` in `trif3`) would become noise as well.
pub(crate) fn r_function_inverse_gap(
    momentum: &Atom,
    mass_i: &Atom,
    mass_j: &Atom,
    sqrt_i: &Atom,
    sqrt_j: &Atom,
    root: &Atom,
) -> Atom {
    let regular = sqrt_i / sqrt_j - Atom::num(1) / root;
    let [nearer_ratio, exact, delta] = zero_momentum_root(momentum, mass_i, mass_j, sqrt_i, sqrt_j);
    let vanishing = &delta / (&exact * (&exact + &delta));
    if_nonzero_else(
        momentum,
        regular.clone(),
        if_nonzero_else(
            &(mass_i - mass_j),
            if_nonzero_else(&nearer_ratio, vanishing, regular),
            Atom::num(0),
        ),
    )
}

/// `[selects s_j/s_i, E, δ]` for the r-function root near `p = 0` with
/// unequal masses `m_i != m_j`.
///
/// At `p = 0`, `k = R + 1/R` with `R = s_j/s_i`, so the root is `E = R` or
/// `1/R`, the one nearest to the computed root (the first entry is 1 when it
/// is `R`). Writing the root as `E + δ`, the quadratic `s_i s_j r² - (m_i +
/// m_j - p) r + s_i s_j = 0` becomes `s_i s_j δ² + B δ + p E = 0` with
/// `B = m_E - m_o + p`, where `m_E = s_i s_j E` is the mass in the numerator of
/// `E` and `m_o` the other one. Its small root
/// `δ = -2 p E / (B (1 + sqrt(1 - 4 p m_E / B²)))` is exactly zero at `p = 0`
/// (`B = m_E - m_o != 0`) and analytic around it.
fn zero_momentum_root(
    momentum: &Atom,
    mass_i: &Atom,
    mass_j: &Atom,
    sqrt_i: &Atom,
    sqrt_j: &Atom,
) -> [Atom; 3] {
    let ratio = sqrt_j / sqrt_i;
    let inverse = sqrt_i / sqrt_j;
    let (root_zero, _) = r_function(&((mass_i + mass_j) / (sqrt_i * sqrt_j)));
    let nearer_ratio = 1 - sheet_exact::negative(
        &(Symbol::ABS.call((&(&root_zero - &inverse),))
            - Symbol::ABS.call((&(&root_zero - &ratio),))),
    );
    let small_root = |exact: &Atom, numerator_mass: &Atom, other_mass: &Atom| {
        let linear = numerator_mass - other_mass + momentum;
        let radicand =
            Atom::num(1) - Atom::num(4) * momentum * numerator_mass / (&linear * &linear);
        Atom::num(-2) * momentum * exact / (&linear * (Atom::num(1) + radicand.sqrt()))
    };
    let delta = if_nonzero_else(
        &nearer_ratio,
        small_root(&ratio, mass_j, mass_i),
        small_root(&inverse, mass_i, mass_j),
    );
    let exact = if_nonzero_else(&nearer_ratio, ratio, inverse);
    [nearer_ratio, exact, delta]
}

pub(crate) fn r_function(value: &Atom) -> (Atom, Atom) {
    // OneLOop selects the root outside the unit circle for non-real q and
    // real |q|>2, but the lower-half-plane root on the real interval [-2,2].
    // Choosing the reciprocal preserves q=r+1/r but changes continued Li2 terms.
    let difference = (value - 2).sqrt() * (value + 2).sqrt();
    let outside = (value + &difference) / 2;
    let discriminant: Atom = value * value - 4;
    // Keep the real-axis formulas on their exact axes: multiplying two
    // numerically evaluated negative-real square roots can leak a tiny
    // imaginary component and incorrectly select a continuation sheet.
    let inside_difference = -sheet_exact::i() * (-&discriminant).sqrt();
    let outside_difference = sheet_exact::sign_nonnegative(value) * discriminant.sqrt();
    let select = |off_axis: Atom, inside: Atom, outside: Atom| {
        if_nonzero_else(
            &sheet_exact::im(value),
            off_axis,
            if_nonzero_else(&sheet_exact::negative(&discriminant), inside, outside),
        )
    };
    // r - 1/r is exactly the square root added to q in each branch. Taking it
    // directly keeps d = 0 exactly at a threshold and its sign reliable next
    // to one, where forming r - 1/r cancels and can even flip the sign.
    let root = select(
        outside,
        (value + &inside_difference) / 2,
        (value + &outside_difference) / 2,
    );
    (
        root,
        select(difference, inside_difference, outside_difference),
    )
}

pub(crate) fn triangle_ir_one(mass: &Atom, mu_squared: &Atom) -> LaurentSeries {
    let factor = Atom::num(1) / (Atom::num(2) * mass);
    LaurentSeries::new(
        &factor * (2 + physical_log(&(mass / mu_squared))),
        -factor,
        Atom::num(0),
    )
}

pub(crate) fn triangle_ir_two(momentum: &Atom, mass: &Atom, mu_squared: &Atom) -> LaurentSeries {
    let difference = mass - momentum;
    let q_mass = SheetAtom::lower(mass.clone());
    let q_difference = SheetAtom::lower(difference.clone());
    let log_mass = q_mass.scale_positive(&(Atom::num(1) / mu_squared)).log();
    let ratio = &q_mass / &q_difference;
    let logarithm = ratio.log();
    let double_pole = Atom::num(1) / 2;
    let single_pole = &logarithm - &double_pole * &log_mass;
    let finite_seed = Atom::num(1) / 24 * Symbol::PI.to_atom() * Symbol::PI.to_atom()
        + &logarithm * &logarithm / 2
        - ratio.dilog();
    let factor = -Atom::num(1) / difference;
    LaurentSeries::new(
        &factor * (finite_seed + (&double_pole * &log_mass / 2 - logarithm) * log_mass),
        &factor * single_pole,
        factor * double_pole,
    )
}

pub(crate) fn triangle_ir_three(
    momentum_2: &Atom,
    momentum_3: &Atom,
    mass: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let difference_13 = mass - momentum_3;
    let difference_23 = mass - momentum_2;
    let q13 = SheetAtom::lower(difference_13.clone());
    let q23 = SheetAtom::lower(difference_23);
    let q_mass = SheetAtom::lower(mass.clone());
    let x1 = &q23 / &q_mass;
    let x2 = &q13 / &q_mass;
    let single_pole = -(&q23 / &q13).log_over_one_minus() / &difference_13;
    LaurentSeries::new(
        -sheet_dilog_divided_difference(&x1, &x2) / mass
            - &single_pole
                * (x1.log() + x2.log() + q_mass.scale_positive(&(Atom::num(1) / mu_squared)).log()),
        single_pole,
        Atom::num(0),
    )
}

pub(crate) fn triangle_ir_four(
    momentum: &Atom,
    mass_2: &Atom,
    mass_3: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let sqrt_2 = physical_sqrt(mass_2);
    let sqrt_3 = physical_sqrt(mass_3);
    let (r23, _) = r_function(&((mass_2 + mass_3 - momentum) / (&sqrt_2 * &sqrt_3)));
    let q23 = SheetAtom::lower(r23.clone());
    let coefficient = q23.log_over_one_minus() * &r23 / ((Atom::num(1) + &r23) * &sqrt_2 * &sqrt_3);
    let ratio = &SheetAtom::upper(sqrt_3.clone()) / &SheetAtom::upper(sqrt_2.clone());
    let finite = &coefficient * (physical_log(&(mass_3 / mu_squared)) - q23.log())
        - sheet_dilog_divided_difference(&(&ratio * &q23), &(&ratio / &q23)) / mass_2
        + sheet_dilog_divided_difference(&(&q23 * &q23), &SheetAtom::lower(Atom::num(1))) * &r23
            / (&sqrt_2 * &sqrt_3);
    LaurentSeries::new(finite, -coefficient, Atom::num(0))
}

pub(crate) fn triangle_finite_massless(p1: &Atom, p2: &Atom, p3: &Atom) -> Atom {
    let r23 = -p1;
    let r24 = -p3;
    let r34 = -p2;
    let a = &r34 * &r24;
    let b = &r24 + &r34 - &r23;
    let (root_1, root_2) = negative_quadratic_roots(&a, &b, &Atom::num(1));
    // The real-axis lips follow r23=-p1 (as in the original trif0
    // continuation), rather than being tied to the quadratic root label.
    // This matters when the real roots swap order across mixed-sign sectors.
    let r23_sign = sheet_exact::sign_nonnegative(&r23);
    let qx1 = SheetAtom::with_sign(root_1.clone(), -&r23_sign);
    let qx2 = SheetAtom::with_sign(root_2.clone(), r23_sign);
    let q23 = SheetAtom::lower(r23);
    let q24 = SheetAtom::lower(r24.clone());
    let q34 = SheetAtom::lower(r34.clone());
    let x1_34 = &qx1 * &q34;
    let x2_34 = &qx2 * &q34;
    let x1_24 = &qx1 * &q24;
    let x2_24 = &qx2 * &q24;
    let ratio = &qx1 / &qx2;
    let root_product = &qx1 * &qx2;
    let log_ratio = ratio.log_over_one_minus();
    (sheet_dilog_divided_difference(&x1_34, &x2_34) * &r34
        + sheet_dilog_divided_difference(&x1_24, &x2_24) * &r24
        - &log_ratio * root_product.log() / (Atom::num(2) * &root_2)
        - log_ratio * q23.log() / root_2)
        / a
}

pub(crate) fn square_root_magnitude(mass: &Atom) -> Atom {
    // sqrt(|m²|) is manifestly real and avoids a complex fourth root.
    Symbol::ABS.call((mass,)).sqrt()
}

fn triangle_one_mass_quadratic_limit(
    p2: &Atom,
    p3: &Atom,
    mass: &Atom,
    quadratic_coefficient: &Atom,
) -> Atom {
    let qmass = SheetAtom::lower(mass.clone());
    let a = &SheetAtom::lower(mass - p3) / &qmass;
    let b = &SheetAtom::lower(mass - p2) / &qmass;
    let inverse_a = &qmass / &SheetAtom::lower(mass - p3);
    let inverse_b = &qmass / &SheetAtom::lower(mass - p2);
    let av = a.value();
    let bv = b.value();
    // At p1=-m*A*B, with A=(m-p3)/m and B=(m-p2)/m, the
    // projective Feynman denominator factors as m*(1+A*x)*(1+B*y).
    // Integrating against 1/(1+x+y) gives this finite dilogarithm limit.
    let numerator =
        a.log() * b.log() + inverse_a.dilog() + inverse_b.dilog() - Symbol::PI.to_atom().pow(2) / 6;
    let denominator = &av + &bv - &av * &bv;
    let regular = numerator / (mass * &denominator);
    // When A+B=A*B, the apparent partial-fraction pole is removable.
    let coincident = -(a.log() + (&av - 1) * b.log()) / (mass * &av);
    let finite = if_nonzero_else(&denominator, regular, coincident);
    // Two negative real factors pinch inside the integration domain. The
    // regulated original can print a large finite number there, but the
    // exact integral is singular and has no regulator-independent limit.
    let both_negative =
        sheet_exact::negative(&sheet_exact::re(&av)) * sheet_exact::negative(&sheet_exact::re(&bv));
    let on_axis = if_nonzero_else(&sheet_exact::im(mass), Atom::num(0), both_negative);
    if_nonzero_else(&on_axis, Atom::num(1) / quadratic_coefficient, finite)
}

pub(crate) fn triangle_finite_one_mass(p1: &Atom, p2: &Atom, p3: &Atom, mass: &Atom) -> Atom {
    let p2_internal = p1;
    let p3_internal = p2;
    let p23 = p3;
    let sqrt_mass = physical_sqrt(mass);
    let sqrt_real = square_root_magnitude(mass);
    let r23 = -p2_internal / (&sqrt_real * &sqrt_real);
    let r24 = (mass - p23) / (&sqrt_real * &sqrt_mass);
    let r34 = (mass - p3_internal) / (&sqrt_real * &sqrt_mass);
    let a = &r34 * &r24 - &r23;
    let b = &r24 / &sqrt_real + &r34 / &sqrt_real - &r23 / &sqrt_mass;
    let c = Atom::num(1) / (&sqrt_real * &sqrt_real);
    let (root_1, root_2) = negative_quadratic_roots(&a, &b, &c);
    // Real roots inherit their lips from OneLOop's p -> p(1 + i eps).
    let dr23 = -p2_internal / (&sqrt_real * &sqrt_real);
    let dr24 = -p23 / (&sqrt_real * &sqrt_mass);
    let dr34 = -p3_internal / (&sqrt_real * &sqrt_mass);
    let da = &dr34 * &r24 + &r34 * &dr24 - &dr23;
    let db = &dr24 / &sqrt_real + &dr34 / &sqrt_real - &dr23 / &sqrt_mass;
    let lip = |root: &Atom| root_lip(root, &a, &b, [&da, &db, &Atom::num(0)]);
    let qx1 = SheetAtom::with_sign(root_1.clone(), lip(&root_1));
    let qx2 = SheetAtom::with_sign(root_2.clone(), lip(&root_2));
    let q23 = SheetAtom::lower(r23.clone());
    let q24 = SheetAtom::lower(r24.clone());
    let q34 = SheetAtom::lower(r34.clone());
    let qmass = SheetAtom::lower(sqrt_mass.clone());
    let root_product = &qx1 * &qx2;
    let mass_product = &qmass * &qmass;
    let root_over_mass = &root_product / &mass_product;
    let ratio = &qx1 / &qx2;
    let mut finite =
        root_over_mass.log() / 2 + q23.scale_positive(&(&sqrt_real * &sqrt_real)).log();
    let x1_mass = &qx1 * &qmass;
    let x2_mass = &qx2 * &qmass;
    finite = -finite * ratio.log_over_one_minus() / &root_2
        - sheet_dilog_divided_difference(&x1_mass, &x2_mass) * &sqrt_mass;
    let scaled_34 = q34.scale_positive(&sqrt_real);
    let x1_34 = &qx1 * &scaled_34;
    let x2_34 = &qx2 * &scaled_34;
    finite += if_nonzero(
        &r34,
        sheet_dilog_divided_difference(&x1_34, &x2_34) * &r34 * &sqrt_real,
    );
    let scaled_24 = q24.scale_positive(&sqrt_real);
    let x1_24 = &qx1 * &scaled_24;
    let x2_24 = &qx2 * &scaled_24;
    finite += if_nonzero(
        &r24,
        sheet_dilog_divided_difference(&x1_24, &x2_24) * &r24 * &sqrt_real,
    );
    // Test the polynomial before its square-root normalization, preserving
    // exact cancellation at rational inputs during runtime substitution.
    let quadratic_coefficient = (mass - p2) * (mass - p3) + p1 * mass;
    if_nonzero_else(
        &quadratic_coefficient,
        finite / (&a * &sqrt_real * &sqrt_real * sqrt_mass),
        triangle_one_mass_quadratic_limit(p2, p3, mass, &quadratic_coefficient),
    )
}

pub(crate) fn triangle_finite_two_masses(
    p1: &Atom,
    p2: &Atom,
    p3: &Atom,
    mass_2_input: &Atom,
    mass_3_input: &Atom,
) -> Atom {
    let p2_internal = p3;
    let p3_internal = p1;
    let p23 = p2;
    let mass_2 = mass_3_input;
    let mass_4 = mass_2_input;
    let sqrt_2 = physical_sqrt(mass_2);
    let sqrt_3 = square_root_magnitude(mass_2);
    let sqrt_4 = physical_sqrt(mass_4);
    let r23 = (mass_2 - p2_internal) / (&sqrt_2 * &sqrt_3);
    let k24 = (mass_2 + mass_4 - p23) / (&sqrt_2 * &sqrt_4);
    let r34 = (mass_4 - p3_internal) / (&sqrt_3 * &sqrt_4);
    let (r24, d24) = r_function_of(p23, mass_2, mass_4, &sqrt_2, &sqrt_4);
    let a = &r34 / &r24 - &r23;
    let b = -&d24 / &sqrt_3 + &r34 / &sqrt_2 - &r23 / &sqrt_4;
    let c = (&sqrt_4 / &sqrt_2 - &r24) / (&sqrt_3 * &sqrt_4);
    // Algebraically b² - 4ac = lambda/(m2 m4 |m2|). Test that polynomial for
    // an exact zero, keeping Gram-degenerate roots exactly coincident. Away
    // from it the rounded b² - 4ac splits nearly coincident roots by about
    // sqrt(eps) instead of eps, which conditions the divided differences far
    // better than the exact but tiny split.
    let lambda = (p3 - p1 - p2).pow(2) - Atom::num(4) * p1 * p2;
    let discriminant = if_nonzero(&lambda, physical_sqrt(&(&b * &b - Atom::num(4) * &a * &c)));
    let (root_1, root_2) = ordered_negative_quadratic_roots(&a, &b, &c, &discriminant);
    // Real roots inherit their lips from OneLOop's p -> p(1 + i eps).
    let dk24 = -p23 / (&sqrt_2 * &sqrt_4);
    let [dr24, dd24, weight] = r_function_velocity(&k24, &r24, &d24, &dk24);
    let dr23 = -p2_internal * &weight / (&sqrt_2 * &sqrt_3);
    let dr34 = -p3_internal * weight / (&sqrt_3 * &sqrt_4);
    let da = &dr34 / &r24 - &r34 * &dr24 / (&r24 * &r24) - &dr23;
    let db = -dd24 / &sqrt_3 + &dr34 / &sqrt_2 - &dr23 / &sqrt_4;
    let dc = -dr24 / (&sqrt_3 * &sqrt_4);
    let lip = |root: &Atom| root_lip(root, &a, &b, [&da, &db, &dc]);
    let qx1 = SheetAtom::with_sign(root_1.clone(), lip(&root_1));
    let qx2 = SheetAtom::with_sign(root_2.clone(), lip(&root_2));
    let q23 = SheetAtom::lower(r23.clone());
    let q24 = SheetAtom::lower(r24.clone());
    let q34 = SheetAtom::lower(r34.clone());
    let qm2 = SheetAtom::lower(sqrt_2.clone());
    let qm3 = SheetAtom::lower(sqrt_3.clone());
    let qm4 = SheetAtom::lower(sqrt_4.clone());
    let y1 = &qx1 / &q24;
    let y2 = &qx2 / &q24;
    let y1m2 = &y1 * &qm2;
    let y2m2 = &y2 * &qm2;
    let mut finite = sheet_dilog_divided_difference(&y1m2, &y2m2) * &sqrt_2 / &r24;
    let y_ratio = &y1 / &y2;
    let y_product = &y1 * &y2;
    let m2_product = &qm2 * &qm2;
    let y_scaled = &y_product / &m2_product;
    let x_ratio = &qx1 / &qx2;
    let x_product = &qx1 * &qx2;
    let m4_product = &qm4 * &qm4;
    let x_scaled = &x_product / &m4_product;
    let root_term = (y_ratio.log_over_one_minus() * y_scaled.log()
        - x_ratio.log_over_one_minus() * x_scaled.log())
        / (Atom::num(2) * &root_2);
    finite += if_nonzero(&root_2, root_term);
    let x1m4 = &qx1 * &qm4;
    let x2m4 = &qx2 * &qm4;
    finite -= sheet_dilog_divided_difference(&x1m4, &x2m4) * &sqrt_4;
    let q23m3 = &q23 * &qm3;
    let scale_23 = &q23m3 / &q24;
    let x1_23 = &qx1 * &scale_23;
    let x2_23 = &qx2 * &scale_23;
    finite -= if_nonzero(
        &r23,
        sheet_dilog_divided_difference(&x1_23, &x2_23) * &r23 * &sqrt_3 / &r24,
    );
    let scale_34 = &q34 * &qm3;
    let x1_34 = &qx1 * &scale_34;
    let x2_34 = &qx2 * &scale_34;
    finite += if_nonzero(
        &r34,
        sheet_dilog_divided_difference(&x1_34, &x2_34) * &r34 * &sqrt_3,
    );
    if_nonzero(&a, finite / (&a * sqrt_2 * sqrt_3 * sqrt_4))
}

pub(crate) fn triangle_finite_three_masses(momenta: [&Atom; 3], masses: [&Atom; 3]) -> Atom {
    three_mass_triangle_call("select", momenta, masses)
}

fn three_mass_triangle_call(stage: &str, p: [&Atom; 3], m: [&Atom; 3]) -> Atom {
    let arguments = p.into_iter().chain(m).cloned().collect::<Vec<_>>();
    symbolica::symbol!(format!("__olo_triangle_three_mass_{stage}")).call(&arguments)
}

fn register_three_mass_triangle(
    map: &mut FunctionMap,
    symbols: &[Symbol],
    p: [&Atom; 3],
    m: [&Atom; 3],
) {
    let rotated_p = [p[2], p[0], p[1]];
    let rotated_m = [m[2], m[0], m[1]];
    let u = m[2] - p[1];
    let v = m[0] - p[0];
    // Eliminate r13 using m1*R²-(m1+m3-p3)*R+m3=0 for
    // R=r13*sqrt(m3/m1). This polynomial detects a zero leading
    // coefficient before rounded square roots can lift its exact zero.
    let resultant = m[0] * &u * &u - (m[0] + m[2] - p[2]) * &u * &v + m[2] * &v * &v;
    // At a vanishing Gram determinant with p3 = 0 (then p1 = p2), c and b
    // can both vanish: the quadratic degenerates to a double root at zero,
    // where the dilogarithm divided differences are singular. Such a chart
    // is skipped like one with a vanishing leading coefficient.
    let lambda = (p[2] - p[0] - p[1]).pow(2) - Atom::num(4) * p[0] * p[1];
    let usable = if_nonzero_else(
        &resultant,
        any_nonzero([lambda, p[2].clone()]),
        Atom::num(0),
    );
    let ordered = three_mass_triangle_call("ordered", p, m);
    let ordered_rotated = three_mass_triangle_call("ordered", rotated_p, rotated_m);
    let second = if_nonzero_else(&usable, ordered.clone(), ordered_rotated);
    let first = if_nonzero_else(
        &usable,
        ordered,
        three_mass_triangle_call("second", rotated_p, rotated_m),
    );
    // Match trif3's cyclic relabeling when the middle propagator has the
    // largest width, including an equality of widths on the real axis.
    let middle_largest = (1 - sheet_exact::negative(&sheet_exact::im(&(m[0] - m[1]))))
        * (1 - sheet_exact::negative(&sheet_exact::im(&(m[2] - m[1]))));
    let select = if_nonzero_else(
        &middle_largest,
        three_mass_triangle_call("first", rotated_p, rotated_m),
        three_mass_triangle_call("first", p, m),
    );
    // A zero leading coefficient can be a degeneration of this particular
    // projective parametrization. Try the remaining cyclic labels at exact
    // zero. Each stage has formal arguments and calls the next native
    // definition, avoiding recursive expansion of selected mass expressions.
    for (stage, body) in [
        ("ordered", triangle_finite_three_masses_ordered(p, m)),
        ("second", second),
        ("first", first),
        ("select", select),
    ] {
        map.add_function_with_options(
            symbolica::symbol!(format!("__olo_triangle_three_mass_{stage}")),
            symbols.to_vec(),
            body,
            native_function_options(),
        )
        .expect("unique finite three-mass triangle stage");
    }
}

fn triangle_finite_three_masses_ordered(momenta: [&Atom; 3], masses: [&Atom; 3]) -> Atom {
    let sqrt_1 = physical_sqrt(masses[0]);
    let sqrt_2 = physical_sqrt(masses[1]);
    let sqrt_3 = physical_sqrt(masses[2]);
    let k12 = (masses[0] + masses[1] - momenta[0]) / (&sqrt_1 * &sqrt_2);
    let k13 = (masses[0] + masses[2] - momenta[2]) / (&sqrt_1 * &sqrt_3);
    let k23 = (masses[1] + masses[2] - momenta[1]) / (&sqrt_2 * &sqrt_3);
    let (r12, _) = r_function_of(momenta[0], masses[0], masses[1], &sqrt_1, &sqrt_2);
    let (r13, d13) = r_function_of(momenta[2], masses[0], masses[2], &sqrt_1, &sqrt_3);
    let (r23, _) = r_function_of(momenta[1], masses[1], masses[2], &sqrt_2, &sqrt_3);
    let a = &sqrt_2 / &sqrt_3 - &k23 + &r13 * (&k12 - &sqrt_2 / &sqrt_1);
    let b = &d13 / &sqrt_2 + &k12 / &sqrt_3 - &k23 / &sqrt_1;
    // Algebraically b²-4ac = lambda/(m1*m2*m3). Evaluating this
    // polynomial identity keeps exact Gram-zero roots coincident instead
    // of introducing an artificial split through roundoff in r13.
    let lambda =
        (momenta[2] - momenta[0] - momenta[1]).pow(2) - Atom::num(4) * momenta[0] * momenta[1];
    let discriminant = physical_sqrt(&(lambda / (masses[0] * masses[1] * masses[2])));
    // The logarithmic root term assumes that an exactly zero root is second,
    // as in OneLOop's solabc. c vanishes identically at p3 = 0 when
    // r13 = sqrt(m3/m1); r_function_inverse_gap keeps that zero exact (the
    // direct difference is rounding noise, and log(0) of a noisy second root
    // is NaN) while keeping c's dependence on every input.
    let c = r_function_inverse_gap(momenta[2], masses[0], masses[2], &sqrt_1, &sqrt_3, &r13)
        / (&sqrt_1 * &sqrt_2);
    let (root_1, root_2) = ordered_negative_quadratic_roots(&a, &b, &c, &discriminant);

    use sheet_exact::{SheetAtom as Q, divided_difference as dd};
    // Real roots inherit their lips from OneLOop's p -> p(1 + i eps).
    let dk12 = -momenta[0] / (&sqrt_1 * &sqrt_2);
    let dk13 = -momenta[2] / (&sqrt_1 * &sqrt_3);
    let dk23 = -momenta[1] / (&sqrt_2 * &sqrt_3);
    let [dr13, dd13, weight] = r_function_velocity(&k13, &r13, &d13, &dk13);
    let dk12 = dk12 * &weight;
    let dk23 = dk23 * weight;
    let da = -&dk23 + &dr13 * (&k12 - &sqrt_2 / &sqrt_1) + &r13 * &dk12;
    let db = dd13 / &sqrt_2 + &dk12 / &sqrt_3 - &dk23 / &sqrt_1;
    let dc = &dr13 / (&r13 * &r13 * &sqrt_1 * &sqrt_2);
    let lip = |root: &Atom| root_lip(root, &a, &b, [&da, &db, &dc]);
    let qx1 = Q::with_sign(root_1.clone(), lip(&root_1));
    let qx2 = Q::with_sign(root_2.clone(), lip(&root_2));
    let [q12, q13, q23, qm1, qm2, qm3] =
        [&r12, &r13, &r23, &sqrt_1, &sqrt_2, &sqrt_3].map(|r| Q::lower(r.clone()));
    let z1 = &qx1 * &qm2;
    let z2 = &qx2 * &qm2;
    let mut finite = (dd(&(&z1 * &q12), &(&z2 * &q12)) * &r12
        + dd(&(&z1 / &q12), &(&z2 / &q12)) / &r12)
        * &sqrt_2;
    let temp = &q13 * &qm2;
    let t1 = &qx1 * &temp;
    let t2 = &qx2 * &temp;
    finite -= (dd(&(&t1 * &q23), &(&t2 * &q23)) * &r23 + dd(&(&t1 / &q23), &(&t2 / &q23)) / &r23)
        * &r13
        * &sqrt_2;
    let y1 = &qx1 * &q13;
    let y2 = &qx2 * &q13;
    finite += dd(&(&y1 * &qm3), &(&y2 * &qm3)) * &r13 * &sqrt_3
        - dd(&(&qx1 * &qm1), &(&qx2 * &qm1)) * &sqrt_1;
    let root_term = ((&y1 / &y2).log_over_one_minus()
        * (y1.log() + y2.log() - Atom::num(2) * qm3.log())
        - (&qx1 / &qx2).log_over_one_minus() * (qx1.log() + qx2.log() - Atom::num(2) * qm1.log()))
        / (Atom::num(2) * &root_2);
    finite += if_nonzero(&root_2, root_term);
    if_nonzero(&a, finite / (&a * sqrt_1 * sqrt_2 * sqrt_3))
}

pub(crate) fn triangle_massless(momenta: [&Atom; 3], mu_squared: &Atom) -> LaurentSeries {
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    let one_scale = |momentum: &Atom| {
        let logarithm = physical_log(&(-momentum / mu_squared));
        LaurentSeries::new(
            (&logarithm * &logarithm / 2 - &pi_squared / 12) / momentum,
            -&logarithm / momentum,
            Atom::num(1) / momentum,
        )
    };
    let two_scale = |first: &Atom, second: &Atom| {
        let first_log = physical_log(&(-first / mu_squared));
        let second_log = physical_log(&(-second / mu_squared));
        let difference = first - second;
        // A logarithm of the ratio loses the individual Feynman lips when
        // the real invariants have opposite signs.
        let single_pole = if_nonzero_else(
            &difference,
            (&second_log - &first_log) / &difference,
            -Atom::num(1) / first,
        );
        LaurentSeries::new(
            -&single_pole * (second_log + first_log) / 2,
            single_pole,
            Atom::num(0),
        )
    };
    let finite = LaurentSeries::new(
        triangle_finite_massless(momenta[0], momenta[1], momenta[2]),
        Atom::num(0),
        Atom::num(0),
    );
    let zero = LaurentSeries::new(Atom::num(0), Atom::num(0), Atom::num(0));

    let branches = [
        zero,
        one_scale(momenta[2]),
        one_scale(momenta[1]),
        two_scale(momenta[1], momenta[2]),
        one_scale(momenta[0]),
        two_scale(momenta[0], momenta[2]),
        two_scale(momenta[0], momenta[1]),
        finite,
    ];
    select_three_series(momenta, branches)
}

pub(crate) fn select_three_series(
    conditions: [&Atom; 3],
    branches: [LaurentSeries; 8],
) -> LaurentSeries {
    let coefficients = branches.map(LaurentSeries::into_coefficients);
    let select = |index: usize| {
        let with_first = if_nonzero_else(
            conditions[1],
            if_nonzero_else(
                conditions[2],
                coefficients[7][index].clone(),
                coefficients[6][index].clone(),
            ),
            if_nonzero_else(
                conditions[2],
                coefficients[5][index].clone(),
                coefficients[4][index].clone(),
            ),
        );
        let without_first = if_nonzero_else(
            conditions[1],
            if_nonzero_else(
                conditions[2],
                coefficients[3][index].clone(),
                coefficients[2][index].clone(),
            ),
            if_nonzero_else(
                conditions[2],
                coefficients[1][index].clone(),
                coefficients[0][index].clone(),
            ),
        );
        if_nonzero_else(conditions[0], with_first, without_first)
    };
    LaurentSeries::new(select(0), select(1), select(2))
}

pub(crate) fn if_series(
    condition: &Atom,
    nonzero: LaurentSeries,
    zero: LaurentSeries,
) -> LaurentSeries {
    let nonzero = nonzero.into_coefficients();
    let zero = zero.into_coefficients();
    LaurentSeries::new(
        if_nonzero_else(condition, nonzero[0].clone(), zero[0].clone()),
        if_nonzero_else(condition, nonzero[1].clone(), zero[1].clone()),
        if_nonzero_else(condition, nonzero[2].clone(), zero[2].clone()),
    )
}

pub(crate) fn finite_series(value: Atom) -> LaurentSeries {
    LaurentSeries::new(value, Atom::num(0), Atom::num(0))
}

pub(crate) fn triangle_one_mass(
    momenta: [&Atom; 3],
    mass: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let finite = finite_series(triangle_finite_one_mass(
        momenta[0], momenta[1], momenta[2], mass,
    ));
    let difference_2 = mass - momenta[1];
    let difference_3 = mass - momenta[2];
    let at_zero_first = if_series(
        &difference_2,
        if_series(
            &difference_3,
            triangle_ir_three(momenta[1], momenta[2], mass, mu_squared),
            triangle_ir_two(momenta[1], mass, mu_squared),
        ),
        if_series(
            &difference_3,
            triangle_ir_two(momenta[2], mass, mu_squared),
            triangle_ir_one(mass, mu_squared),
        ),
    );
    if_series(momenta[0], finite, at_zero_first)
}

pub(crate) fn triangle_two_masses(
    momenta: [&Atom; 3],
    mass_2: &Atom,
    mass_3: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let finite = finite_series(triangle_finite_two_masses(
        momenta[0], momenta[1], momenta[2], mass_2, mass_3,
    ));
    let infrared = triangle_ir_four(momenta[1], mass_2, mass_3, mu_squared);
    let difference_1 = momenta[0] - mass_2;
    let difference_3 = momenta[2] - mass_3;
    if_series(
        &difference_1,
        finite.clone(),
        if_series(&difference_3, finite, infrared),
    )
}

pub(crate) fn with_triangle_normalization(series: LaurentSeries) -> LaurentSeries {
    let result = series.into_coefficients();
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    LaurentSeries::new(
        &result[0] + &pi_squared * &result[2] / 12,
        result[1].clone(),
        result[2].clone(),
    )
}

/// Constructs a scalar triangle with its transparent native function definitions.
/// Momenta are [p1²,p2²,(p1+p2)²]; masses follow propagator order.
///
/// Exact one-mass quadratic limits are supported. Other quadratic/Gram
/// degeneracies are not generally supported.
/// See STATUS.md for the development implementation's coverage limitations.
pub fn c0(p: [&Atom; 3], m: [&Atom; 3], mu_squared: &Atom) -> MappedLaurentSeries {
    let expressions = OneLoopExpressions::new();
    MappedLaurentSeries {
        series: expressions.c0(p, m, mu_squared),
        function_map: expressions.into_definitions(),
    }
}

/// Explicit alias for [c0].
pub fn c0_mapped(p: [&Atom; 3], m: [&Atom; 3], mu_squared: &Atom) -> MappedLaurentSeries {
    c0(p, m, mu_squared)
}

pub(crate) fn register_triangles(map: &mut FunctionMap) {
    let symbols =
        core::array::from_fn::<_, 7, _>(|i| symbolica::symbol!(format!("__olo_c0_arg_{i}")));
    let x = symbols.map(Symbol::to_atom);
    let p = [&x[0], &x[1], &x[2]];
    let m = [&x[3], &x[4], &x[5]];
    register_three_mass_triangle(map, &symbols[..6], p, m);
    let lambda = (p[2] - p[0] - p[1]).pow(2) - Atom::num(4) * p[0] * p[1];
    let generic = if_nonzero_else(
        &sheet_exact::negative(&lambda),
        triangle_hv::expression(p, m),
        triangle_finite_three_masses(p, m),
    );
    let at_zero = triangle_vacuum(m);
    let nonzero_p = any_nonzero(p.into_iter().cloned());
    let definitions = [
        (0, triangle_massless(p, &x[6])),
        (1, triangle_one_mass(p, m[2], &x[6])),
        (3, triangle_two_masses(p, m[1], m[2], &x[6])),
        (
            7,
            finite_series(if_nonzero_else(&nonzero_p, generic, at_zero)),
        ),
    ];
    for (sector, series) in definitions {
        for (coefficient, body) in series.into_coefficients().into_iter().enumerate() {
            map.add_function_with_options(
                symbolica::symbol!(format!("__olo_c0_sector_{sector}_{coefficient}")),
                symbols.to_vec(),
                body,
                native_function_options(),
            )
            .expect("unique triangle sector");
        }
    }
    for (sector, canonical, shift) in [(2, 1, 2), (4, 1, 1), (5, 3, 1), (6, 3, 2)] {
        let args = (0..3)
            .map(|i| x[(i + shift) % 3].clone())
            .chain((0..3).map(|i| x[3 + (i + shift) % 3].clone()))
            .chain([x[6].clone()])
            .collect::<Vec<_>>();
        for coefficient in 0..3 {
            let body = symbolica::symbol!(format!("__olo_c0_sector_{canonical}_{coefficient}"))
                .call(&args);
            map.add_function(
                symbolica::symbol!(format!("__olo_c0_sector_{sector}_{coefficient}")),
                symbols.to_vec(),
                body,
            )
            .unwrap();
        }
    }
}

fn triangle_vacuum(m: [&Atom; 3]) -> Atom {
    let pair = |a: &Atom, b: &Atom, c: &Atom| {
        // Continue each mass before forming a ratio. A principal ratio log
        // loses the lower lip when its denominator is negative real.
        let [qa, qb, qc] = [a, b, c].map(|mass| SheetAtom::lower(mass.clone()));
        ((&qc / &qa).log_over_one_minus() - (&qc / &qb).log_over_one_minus()) / (a - b)
    };
    let unequal = if_nonzero_else(
        &(m[0] - m[1]),
        pair(m[0], m[1], m[2]),
        pair(m[1], m[2], m[0]),
    );
    let difference = any_nonzero([m[0] - m[1], m[1] - m[2]]);
    if_nonzero_else(&difference, unequal, -Atom::num(1) / (Atom::num(2) * m[0]))
}
