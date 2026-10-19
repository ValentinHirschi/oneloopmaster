//! Box mass sectors and infrared formulas.
use super::*;
use sheet_exact::sqrt_lower as physical_sqrt;
use sheet_exact::{
    SheetAtom as Q, complement_difference, divided_difference as dd, im, negative, re,
    sign_nonnegative,
};

pub(crate) fn box_finite_massless(p: [&Atom; 6]) -> Atom {
    let r12 = -p[0];
    let r13 = -p[4];
    let r14 = -p[3];
    let r23 = -p[1];
    let r24 = -p[5];
    let r34 = -p[2];
    let a = &r34 * &r24;
    let b = &r13 * &r24 + &r12 * &r34 - &r14 * &r23;
    let c = &r12 * &r13;
    // OneLOop's boxf0 gives the real root (b - d)/(2a) the lip of Re(r23)
    // and (b + d)/(2a) the opposite one, with d >= 0 for a real positive
    // discriminant; the formula divides by the latter root.
    // Both roots are cancellation-free (Vieta), keeping their labels.
    let gap = channel_gap([&r12, &r13, &r24, &r34]);
    let sum = &r13 * &r24 + &r12 * &r34;
    let square = collinear_square(&gap, sum, &r14 * &r23, Atom::num(0));
    let discriminant = physical_sqrt(&square);
    let (root_2, root_1) = labelled_quadratic_roots(&a, &(-&b), &c, &discriminant);
    let qx1 = Q::with_sign(root_1.clone(), sign_nonnegative(&r23));
    let qx2 = Q::with_sign(root_2.clone(), -sign_nonnegative(&r23));
    let [q12, q13, q14, q23, q24, q34] =
        [&r12, &r13, &r14, &r23, &r24, &r34].map(|r| Q::lower(r.clone()));
    // 1 - X/u for u = r13/r34 and r12/r24, where the quadratic is r14 r23 u:
    // next to a collinear limit (r14 or r23 small) one argument of each pair
    // is next to one ([root_complements]). The offsets X - u are the roots
    // of the shifted quadratic, with the roots' labels.
    let complements = |numerator: &Atom, denominator: &Atom, shifted_b: Atom| {
        let shifted_c = &r14 * &r23 * numerator / denominator;
        let (offset_2, offset_1) =
            labelled_quadratic_roots(&a, &(-shifted_b), &shifted_c, &discriminant);
        let scale = denominator / numerator;
        [-offset_1 * &scale, -offset_2 * scale]
    };
    let [w1, w2] = complements(&r13, &r34, -&gap - &r14 * &r23);
    let [y1, y2] = [&qx1, &qx2].map(|x| x * &(&q34 / &q13));
    let mut finite = complement_difference(&y1, &w1, &y2, &w2) * &r34 / &r13;
    let [w1, w2] = complements(&r12, &r24, &gap - &r14 * &r23);
    let [y1, y2] = [&qx1, &qx2].map(|x| x * &(&q24 / &q12));
    finite += complement_difference(&y1, &w1, &y2, &w2) * &r24 / &r12;
    let logarithmic = -(&qx1 / &qx2).log_over_one_minus() / &root_2;
    finite +=
        logarithmic * ((qx1.log() + qx2.log()) / 2 - q12.log() - q13.log() + q14.log() + q23.log());
    if_nonzero(&(&r13 * &a), -finite / a)
}

// OneLOop selects its contour-based representation for complex masses and
// these timelike configurations. Exact inequalities replace its numerical
// proximity tests. Keep the boxf representation for kinematics boxc rejects.
//
// Four real masses keep this selection (the contour box costs 1.5 to 4 times
// as much on the recorded NinjaRS streams); their boxf4 chart is corrected in
// `box_four_masses_chart`.
fn finite_box_representation(p: [&Atom; 6], m: [&Atom; 4], fallback: Atom) -> Atom {
    let complex_mass = any_nonzero(m.map(sheet_exact::im));
    let all_external_nonnegative = p[..4]
        .iter()
        .fold(Atom::num(1), |c, v| c * (1 - negative(v)));
    let channels_nonnegative = (1 - negative(p[4])) * (1 - negative(p[5]));
    let use_contour = any_nonzero([complex_mass, all_external_nonnegative, channels_nonnegative])
        * box_complex::supported_kinematics(p);
    if_nonzero_else(&use_contour, box_complex::expression(p, m), fallback)
}

// Three massive lines take the contour representation (at the vanishing-width
// limit for real masses) wherever it supports the kinematics. Outside
// OneLOop's timelike configurations its real-mass boxf33 chart depends on the
// labelling: in a sweep of 154 random three-mass boxes, 39 had a wrong value
// in some of their 24 relabellings, as had the NinjaRS box of
// tests/d0_real_mass_orderings.rs in 8 (OneLOop 3.6 included).
fn finite_three_mass_representation(p: [&Atom; 6], m: [&Atom; 4], fallback: Atom) -> Atom {
    if_nonzero_else(
        &box_complex::supported_kinematics(p),
        box_complex::expression(p, m),
        fallback,
    )
}

/// OneLOop's boxf4 as it stands, for complex masses (where the contour box
/// rejects the kinematics) and for real masses without any invariant at or
/// below its pseudo-threshold, where [box_four_masses_chart] does not apply
/// (all invariants positive and all four Källén functions negative: no real
/// momenta, and OneLOop returns its error value).
fn box_four_masses_dns(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    let roots = m.map(physical_sqrt);
    let k12 = (m[0] + m[1] - p[0]) / (&roots[0] * &roots[1]);
    let k13 = (m[0] + m[2] - p[4]) / (&roots[0] * &roots[2]);
    let k14 = (m[0] + m[3] - p[3]) / (&roots[0] * &roots[3]);
    let k23 = (m[1] + m[2] - p[1]) / (&roots[1] * &roots[2]);
    let k24 = (m[1] + m[3] - p[5]) / (&roots[1] * &roots[3]);
    let k34 = (m[2] + m[3] - p[2]) / (&roots[2] * &roots[3]);
    let (r12, _) = r_function(&k12);
    let (r13, d13) = r_function(&k13);
    let (r14, _) = r_function(&k14);
    let (r23, _) = r_function(&k23);
    let (r24, d24) = r_function(&k24);
    let (r34, _) = r_function(&k34);
    let a = &k34 / &r24 + &r13 * &k12 - &k14 * &r13 / &r24 - &k23;
    let b = &d13 * &d24 + &k12 * &k34 - &k14 * &k23;
    let c = &k12 / &r13 + &r24 * &k34 - &k14 * &r24 / &r13 - &k23;
    let (x1, x2) = negative_quadratic_roots(&a, &b, &c);
    // Real roots inherit their lips from OneLOop's p -> p(1 + i eps); its
    // explicit sign applies only when that first-order shift vanishes.
    let dk = |p: &Atom, i: usize, j: usize| -p / (&roots[i] * &roots[j]);
    let (dk12, dk13, dk14) = (dk(p[0], 0, 1), dk(p[4], 0, 2), dk(p[3], 0, 3));
    let (dk23, dk24, dk34) = (dk(p[1], 1, 2), dk(p[5], 1, 3), dk(p[2], 2, 3));
    let ([[dr13, dd13], [dr24, dd24]], weight) =
        r_function_velocities([[&k13, &r13, &d13, &dk13], [&k24, &r24, &d24, &dk24]]);
    let [dk12, dk14, dk23, dk34] = [dk12, dk14, dk23, dk34].map(|v| v * &weight);
    let da = &dk34 / &r24 - &k34 * &dr24 / (&r24 * &r24) + &dr13 * &k12 + &r13 * &dk12
        - &dk14 * &r13 / &r24
        - &k14 * &dr13 / &r24
        + &k14 * &r13 * &dr24 / (&r24 * &r24)
        - &dk23;
    let db =
        &dd13 * &d24 + &d13 * &dd24 + &dk12 * &k34 + &k12 * &dk34 - &dk14 * &k23 - &k14 * &dk23;
    let dc = &dk12 / &r13 - &k12 * &dr13 / (&r13 * &r13) + &dr24 * &k34 + &r24 * &dk34
        - &dk14 * &r24 / &r13
        - &k14 * &dr24 / &r13
        + &k14 * &r24 * &dr13 / (&r13 * &r13)
        - &dk23;
    let base = re(&(&k23 - &r13 * &k12 - &r24 * &k34 + &r13 * &r24 * &k14)) * re(&a);
    let lip = |root: &Atom, other: &Atom| {
        let fallback = sign_nonnegative(&(&base * re(other)));
        root_lip_or(root, &a, &b, [&da, &db, &dc], fallback)
    };
    let qx1 = Q::with_sign(x1.clone(), lip(&x1, &x2));
    let qx2 = Q::with_sign(x2.clone(), lip(&x2, &x1));
    let [q12, q13, q14, q23, q24, q34] =
        [&r12, &r13, &r14, &r23, &r24, &r34].map(|r| Q::lower(r.clone()));
    let y1 = &qx1 / &q24;
    let y2 = &qx2 / &q24;
    let mut finite =
        (dd(&(&y1 * &q12), &(&y2 * &q12)) * &r12 + dd(&(&y1 / &q12), &(&y2 / &q12)) / &r12) / &r24;
    let ratio = &r13 / &r24;
    let qr = Q::with_sign(ratio.clone(), -sign_nonnegative(&re(&r24)));
    let z1 = &qx1 * &qr;
    let z2 = &qx2 * &qr;
    finite -= (dd(&(&z1 * &q23), &(&z2 * &q23)) * &r23 + dd(&(&z1 / &q23), &(&z2 / &q23)) / &r23)
        * &ratio;
    let w1 = &qx1 * &q13;
    let w2 = &qx2 * &q13;
    finite +=
        (dd(&(&w1 * &q34), &(&w2 * &q34)) * &r34 + dd(&(&w1 / &q34), &(&w2 / &q34)) / &r34) * &r13;
    finite -= dd(&(&qx1 * &q14), &(&qx2 * &q14)) * &r14 + dd(&(&qx1 / &q14), &(&qx2 / &q14)) / &r14;
    if_nonzero(
        &a,
        -finite / (&a * &roots[0] * &roots[1] * &roots[2] * &roots[3]),
    )
}

/// The boxf4 chart for four real masses with `r24 >= 1`, i.e. an invariant
/// `p23 <= (m2 - m4)^2` in the `p23` slot: the representation of Denner and
/// Dittmaier (arXiv:1005.2076, eq. (3.61), their `r20`), which is OneLOop's
/// boxf4 term by term.
///
/// OneLOop's boxf4 depended on the labelling (wrong in about 2% of the
/// relabellings of random real four-mass boxes, OneLOop 3.6 included). The
/// representation requires that r-function to be real and positive, which
/// OneLOop does not arrange, and exactly real: its infinitesimal imaginary
/// part moves `x r12 / r24` across the cut where `r12` is real and negative.
/// It also forms `x r13` as a principal value and restores the winding with
/// the explicit term `eta(x, r13) Log R`; OneLOop's phase-tracked three-factor
/// products differ from that by multiples of `(2 pi i)^2`. The other
/// infinitesimals keep OneLOop's `p -> p (1 + i eps)` on every invariant except
/// `p23`.
///
/// Every infinitesimal decision is exact: `D = b^2 - 4 a c` is real (a
/// multiple of the modified Cayley determinant), so with a real `r13` the
/// roots are real exactly when `D >= 0`, and with `r13` on the unit circle
/// (complex, or `-1` or `1` at a threshold) `R` is real when `D < 0` as well.
fn box_four_masses_chart(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    let roots = m.map(physical_sqrt);
    let k12 = (m[0] + m[1] - p[0]) / (&roots[0] * &roots[1]);
    let k13 = (m[0] + m[2] - p[4]) / (&roots[0] * &roots[2]);
    let k14 = (m[0] + m[3] - p[3]) / (&roots[0] * &roots[3]);
    let k23 = (m[1] + m[2] - p[1]) / (&roots[1] * &roots[2]);
    let k24 = (m[1] + m[3] - p[5]) / (&roots[1] * &roots[3]);
    let k34 = (m[2] + m[3] - p[2]) / (&roots[2] * &roots[3]);
    let (r12, _) = r_function(&k12);
    let (r13, d13) = r_function(&k13);
    let (r14, _) = r_function(&k14);
    let (r23, _) = r_function(&k23);
    let (r24, d24) = r_function(&k24);
    let (r34, _) = r_function(&k34);
    let a = &k34 / &r24 + &r13 * &k12 - &k14 * &r13 / &r24 - &k23;
    let b = &d13 * &d24 + &k12 * &k34 - &k14 * &k23;
    let c = &k12 / &r13 + &r24 * &k34 - &k14 * &r24 / &r13 - &k23;
    let (x1, x2) = negative_quadratic_roots(&a, &b, &c);
    // First-order velocities of p -> p (1 + i eps) with p23 static.
    let dk = |p: &Atom, i: usize, j: usize| -p / (&roots[i] * &roots[j]);
    let [dr13, dd13, weight] = r_function_velocity(&k13, &r13, &d13, &dk(p[4], 0, 2));
    // At a threshold r13 = r0 = -1 or 1 the root moves along (r0 - i) sqrt(eps),
    // off the unit circle as well, where R is real: D(r13) is proportional to
    // -(1 + i r0). Real decisions take only its real part, the velocity above.
    let tilt = if_nonzero_else(&weight, Atom::num(1), 1 + sheet_exact::i() * &r13);
    let [dr13, dd13] = [dr13 * &tilt, dd13 * &tilt];
    let [dk12, dk14, dk23, dk34] = [(p[0], 0, 1), (p[3], 0, 3), (p[1], 1, 2), (p[2], 2, 3)]
        .map(|(v, i, j)| dk(v, i, j) * &weight);
    let r13_2 = &r13 * &r13;
    let da = &dk34 / &r24 + &dr13 * &k12 + &r13 * &dk12
        - &dk14 * &r13 / &r24
        - &k14 * &dr13 / &r24
        - &dk23;
    let db = &dd13 * &d24 + &dk12 * &k34 + &k12 * &dk34 - &dk14 * &k23 - &k14 * &dk23;
    let dc = &dk12 / &r13 - &k12 * &dr13 / &r13_2 + &r24 * &dk34 - &dk14 * &r24 / &r13
        + &k14 * &r24 * &dr13 / &r13_2
        - &dk23;
    let base = re(&(&k23 - &r13 * &k12 - &r24 * &k34 + &r13 * &r24 * &k14)) * re(&a);
    let lip = |root: &Atom, other: &Atom| {
        let fallback = sign_nonnegative(&(&base * re(other)));
        root_lip_or(root, &a, &b, [&da, &db, &dc], fallback)
    };
    let lips = [lip(&x1, &x2), lip(&x2, &x1)];
    // r13 is real for |k13| >= 2 and on the unit circle for |k13| <= 2.
    let r13_discriminant = &k13 * &k13 - 4;
    let r13_real = 1 - negative(&r13_discriminant);
    let r13_circle = 1 - negative(&-&r13_discriminant);
    let below = negative(&re(&(&b * &b - Atom::num(4) * &a * &c)));
    let x_real = &r13_real * (1 - &below);
    let ratio_real = &x_real + r13_circle * &below;
    let sign_or =
        |value: &Atom, fallback: Atom| if_nonzero_else(value, sign_nonnegative(value), fallback);
    // The direction of r13 (zero for a static p12, where eta vanishes).
    let s13 = if_nonzero_else(
        &r13_real,
        sign_or(&re(&dr13), Atom::num(0)),
        sign_nonnegative(&im(&r13)),
    );
    let two_pi_i = Atom::num(2) * Symbol::PI.to_atom() * sheet_exact::i();
    let mut principal = Vec::new();
    let mut winding = Vec::new();
    let mut log_slope = Atom::num(0);
    for (index, (x, lip)) in [(&x1, &lips[0]), (&x2, &lips[1])].into_iter().enumerate() {
        let slope = Atom::num(2) * &a * x - &b;
        let dx = if_nonzero(&slope, (&db * x - &da * x * x - &dc) / &slope);
        let v = x * &r13;
        let dv = &dx * &r13 + x * &dr13;
        // The lip of the principal value x r13, real exactly with x.
        let lv = if_nonzero_else(
            &x_real,
            sign_or(&re(&dv), lip * sign_nonnegative(&re(&r13))),
            sign_nonnegative(&im(&v)),
        );
        // eta(x, r13) / (2 pi i) from the three imaginary-part signs.
        let sx = if_nonzero_else(&x_real, lip.clone(), sign_nonnegative(&im(x)));
        let n = -&sx * &s13 * &s13 * (1 + &sx * &s13) * (1 - &sx * &lv) / 4;
        // R = P(x) of the r13 terms; its side, when real, from its velocity.
        let u = &v / &r24;
        let numerator = 1 - &u * &k23 + &u * &u;
        let denominator = 1 - &v * &k34 + &v * &v;
        let du = &dv / &r24;
        let dnumerator = -&du * &k23 - &u * &dk23 + Atom::num(2) * &u * &du;
        let ddenominator = -&dv * &k34 - &v * &dk34 + Atom::num(2) * &v * &dv;
        let ratio = &numerator / &denominator;
        let dratio = &ratio * (&dnumerator / &numerator - &ddenominator / &denominator);
        let side = sign_or(&re(&dratio), Atom::num(1));
        let value = if_nonzero_else(&ratio_real, re(&ratio), ratio);
        if index == 0 {
            // d Log R / dx at a double root (x1 = x2).
            let slope = (Atom::num(2) * &u - &k23) / (&r24 * &numerator)
                - (Atom::num(2) * &v - &k34) / &denominator;
            log_slope = if_nonzero(&n, &n * &r13 * slope);
        }
        winding.push(if_nonzero(&n, &n * Q::with_sign(value, side).log()));
        principal.push(Q::with_sign(v, lv));
    }
    let qx1 = Q::with_sign(x1.clone(), lips[0].clone());
    let qx2 = Q::with_sign(x2.clone(), lips[1].clone());
    let [q12, q14, q23, q24, q34] = [&r12, &r14, &r23, &r24, &r34].map(|r| Q::lower(r.clone()));
    let y1 = &qx1 / &q24;
    let y2 = &qx2 / &q24;
    let mut finite =
        (dd(&(&y1 * &q12), &(&y2 * &q12)) * &r12 + dd(&(&y1 / &q12), &(&y2 / &q12)) / &r12) / &r24;
    let ratio = &r13 / &r24;
    let inverse24 = Atom::num(1) / &r24;
    let z1 = principal[0].scale_positive(&inverse24);
    let z2 = principal[1].scale_positive(&inverse24);
    finite -= (dd(&(&z1 * &q23), &(&z2 * &q23)) * &r23 + dd(&(&z1 / &q23), &(&z2 / &q23)) / &r23)
        * &ratio;
    let [w1, w2] = [&principal[0], &principal[1]];
    finite +=
        (dd(&(w1 * &q34), &(w2 * &q34)) * &r34 + dd(&(w1 / &q34), &(w2 / &q34)) / &r34) * &r13;
    finite -= dd(&(&qx1 * &q14), &(&qx2 * &q14)) * &r14 + dd(&(&qx1 / &q14), &(&qx2 / &q14)) / &r14;
    // The explicit term, -sum_k (-1)^(k+1) eta_k Log R_k / (x1 - x2).
    let gap = &x1 - &x2;
    finite -= two_pi_i * if_nonzero_else(&gap, (&winding[0] - &winding[1]) / &gap, log_slope);
    if_nonzero(
        &a,
        -finite / (&a * &roots[0] * &roots[1] * &roots[2] * &roots[3]),
    )
}

fn box_four_masses_chart_symbol() -> Symbol {
    symbolica::symbol!("__olo_box_four_masses_chart")
}

/// Relabellings of the propagators that move each invariant into the `p23`
/// slot (momenta and masses in the new order): the channels `p23`, `p12`, then
/// the external `p1` to `p4`.
pub(crate) const FOUR_MASS_CHART_ORDERS: [([usize; 6], [usize; 4]); 6] = [
    ([0, 1, 2, 3, 4, 5], [0, 1, 2, 3]),
    ([1, 2, 3, 0, 5, 4], [1, 2, 3, 0]),
    ([4, 3, 5, 1, 2, 0], [2, 0, 3, 1]),
    ([0, 5, 2, 4, 3, 1], [0, 1, 3, 2]),
    ([4, 1, 5, 3, 0, 2], [0, 2, 1, 3]),
    ([0, 4, 2, 5, 1, 3], [1, 0, 2, 3]),
];

pub(crate) fn box_finite_four_masses(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    let roots = m.map(physical_sqrt);
    let real_masses = 1 - any_nonzero(m.map(sheet_exact::im));
    // The first invariant whose r-function is real and at least one (real
    // k >= 2, at or below its pseudo-threshold) takes the p23 slot of
    // box_four_masses_chart; k is formed as the chart forms k24.
    let generic = FOUR_MASS_CHART_ORDERS.iter().rev().fold(
        box_four_masses_dns(p, m),
        |otherwise, (order, masses)| {
            let k = (m[masses[1]] + m[masses[3]] - p[order[5]])
                / (&roots[masses[1]] * &roots[masses[3]]);
            let positive = &real_masses * (1 - negative(&(&k * &k - 4))) * (1 - negative(&k));
            let args: Vec<Atom> = order
                .map(|i| p[i].clone())
                .into_iter()
                .chain(masses.map(|i| m[i].clone()))
                .collect();
            if_nonzero_else(
                &positive,
                box_four_masses_chart_symbol().call(&args),
                otherwise,
            )
        },
    );
    // For equal masses at zero external momentum, the Feynman-parameter
    // denominator is constant and the simplex volume is exactly 1/6.
    let nonvacuum = any_nonzero(
        p.into_iter()
            .cloned()
            .chain(m[1..].iter().map(|mass| *mass - m[0])),
    );
    if_nonzero_else(
        &nonvacuum,
        finite_box_representation(p, m, generic),
        Atom::num(1) / (Atom::num(6) * m[0] * m[0]),
    )
}

/// `N = P - R` for the channel product `P = r13 r24` and the external
/// product `R = r12 r34` of the finite boxes' quadratics, formed once.
///
/// Next to a collinear limit (an invariant between two massless lines that
/// is small but nonzero) and to `P = R`, the discriminant is next to `N²`,
/// and the two dilogarithm terms grow as `1/N` and cancel each other (by 560
/// in a two-mass box with `N = 2e-3 P`, which binary64 kept to 1e-11 for
/// collinear invariants of 1e-5 to 1e-12). Their shifted coefficients
/// `b - 2 a u` are `-N` and `N` up to the collinear terms; summed with those
/// in the evaluator's own order, each term rounded `N` differently, and the
/// residue of their cancellation was the result's error. Taken as one
/// value, `N` cancels between them as it does analytically.
fn channel_gap([r12, r13, r24, r34]: [&Atom; 4]) -> Atom {
    grouped(&(r13 * r24 - r12 * r34))
}

/// `b² - 4ac` of the finite boxes' quadratics with `b = P + R - T`, the
/// collinear term `T` (`r14 r23`, `k14 r23` with two masses) and
/// `4ac = 4PR - rest`, as `N² - T (2 (P + R) - T) + rest` with `N = P - R`
/// ([channel_gap]) and `sum = P + R`.
///
/// Next to a collinear limit `T` and `rest` are small, and `b² - 4ac`
/// cancels by `b²/d²` where `P` is next to `R`. There the dilogarithm terms
/// grow as `1/d` and cancel each other, and the rounding residue of
/// `b² - 4ac` became the result's: of 300 one- and two-adjacent-mass boxes
/// with p1² from 1e-14 to 1e-4 and `P/R - 1` from 1e-6 to 0.1, 16 lost 1e-10
/// to 3e-5 in binary64 (a massless box with `P/R - 1 = 2e-3` 1e-7).
fn collinear_square(gap: &Atom, sum: Atom, collinear: Atom, rest: Atom) -> Atom {
    gap.pow(2) - &collinear * (Atom::num(2) * sum - &collinear) + rest
}

/// Roots `X1, X2` of `a X² - b X + c` as OneLOop's `solabc` returns them
/// (negated): the larger directly, the other from Vieta's product
/// ([ordered_negative_quadratic_roots]); and the square root of the
/// discriminant `square`, `b² - 4ac` in a form exact for the caller.
///
/// Next to a collinear limit (an invariant between two massless lines that
/// is small but nonzero) `a` or `c` is small and one root is small; formed
/// directly, `(b - d)/(2a)` lost it to cancellation (a one-mass box with
/// p2² = -1.3e-10 kept 4e-9 in binary64, OneLOop 8e-14).
fn solabc_roots(a: &Atom, b: &Atom, c: &Atom, square: &Atom) -> (Atom, Atom, Atom) {
    let discriminant = physical_sqrt(square);
    let (x1, x2) = ordered_negative_quadratic_roots(a, b, c, &discriminant);
    (x1, x2, discriminant)
}

/// The complements `1 - s X1, 1 - s X2` of the dilogarithm arguments `s X`,
/// from the offsets of the roots of `a X² - b X + c` from `u = 1/s`
/// ([ordered_root_offsets]), given `b - 2 a u` and the quadratic at `u`.
///
/// When the invariant `r23` between the two massless lines of a one-mass or
/// two-adjacent-mass box is small, the quadratic tends to `(r34 X - r13)
/// (r24 X - r12)`, and `X r34/r13` and `X r24/r12` approach one: their
/// differences from one are proportional to `r23`, which the closed forms of
/// the quadratic at `u` keep ([complement_difference]).
fn root_complements(
    [a, b, discriminant]: [&Atom; 3],
    scale: &Atom,
    shifted_b: &Atom,
    shifted_c: &Atom,
) -> [Atom; 2] {
    let (first, second) = ordered_root_offsets(a, b, discriminant, shifted_b, shifted_c);
    [-first * scale, -second * scale]
}

pub(crate) fn box_finite_one_mass(p: [&Atom; 6], mass: &Atom) -> Atom {
    let sm = physical_sqrt(mass);
    let sr = square_root_magnitude(mass);
    let r12 = (mass - p[3]) / (&sr * &sm);
    let r13 = (mass - p[5]) / (&sr * &sm);
    let r14 = (mass - p[2]) / (&sr * &sm);
    let r23 = -p[0] / (&sr * &sr);
    let r24 = -p[4] / (&sr * &sr);
    let r34 = -p[1] / (&sr * &sr);
    let a = &r34 * &r24;
    let b = &r13 * &r24 + &r12 * &r34 - &r14 * &r23;
    let c = &r12 * &r13 - &r23;
    let gap = channel_gap([&r12, &r13, &r24, &r34]);
    let square = collinear_square(
        &gap,
        &r13 * &r24 + &r12 * &r34,
        &r14 * &r23,
        Atom::num(4) * &r23 * &r24 * &r34,
    );
    let (x1, x2, discriminant) = solabc_roots(&a, &b, &c, &square);
    // Real roots inherit their lips from OneLOop's p -> p(1 + i eps).
    let massive = &sr * &sm;
    let massless = &sr * &sr;
    let [d12, d13, d14] = [p[3], p[5], p[2]].map(|v| -v / &massive);
    let [d23, d24, d34] = [p[0], p[4], p[1]].map(|v| -v / &massless);
    let da = &d34 * &r24 + &r34 * &d24;
    let db = &d13 * &r24 + &r13 * &d24 + &d12 * &r34 + &r12 * &d34 - &d14 * &r23 - &r14 * &d23;
    let dc = &d12 * &r13 + &r12 * &d13 - &d23;
    let lip = |root: &Atom| root_lip(root, &a, &b, [&da, &db, &dc]);
    let qx1 = Q::with_sign(x1.clone(), lip(&x1));
    let qx2 = Q::with_sign(x2.clone(), lip(&x2));
    let [q12, q13, q14, q23, q24, q34] =
        [&r12, &r13, &r14, &r23, &r24, &r34].map(|r| Q::lower(r.clone()));
    let lr = (&qx1 / &qx2).log_over_one_minus();
    let root_log = qx1.log() + qx2.log();
    let both_zero = &lr * (&root_log + q34.log() + q24.log() - q23.log()) / &x2;
    let r13_zero = &lr * (&root_log / 2 + q34.log() + q12.log() - q23.log()) / &x2;
    // At u = r13/r34 and r12/r24 the quadratic is r23 (r14 u - 1).
    let quadratic = [&a, &b, &discriminant];
    let w13 = root_complements(
        quadratic,
        &(&r34 / &r13),
        &(-&gap - &r14 * &r23),
        &(&r23 * (&r13 * &r14 / &r34 - 1)),
    );
    let [y1, y2] = [&qx1, &qx2].map(|x| x * &(&q34 / &q13));
    let r13_nonzero = complement_difference(&y1, &w13[0], &y2, &w13[1]) * &r34 / &r13;
    let r12_zero = &lr * (&root_log / 2 + q24.log() + q13.log() - q23.log()) / &x2;
    let w12 = root_complements(
        quadratic,
        &(&r24 / &r12),
        &(&gap - &r14 * &r23),
        &(&r23 * (&r12 * &r14 / &r24 - 1)),
    );
    let [y1, y2] = [&qx1, &qx2].map(|x| x * &(&q24 / &q12));
    let r12_nonzero = complement_difference(&y1, &w12[0], &y2, &w12[1]) * &r24 / &r12;
    // As in OneLOop's boxf1, a single vanishing r12 or r13 replaces only its
    // own dilogarithm term by the logarithmic limit; the other term remains.
    let term_13 = if_nonzero_else(&r13, r13_nonzero, r13_zero);
    let term_12 = if_nonzero_else(&r12, r12_nonzero, r12_zero);
    let both_logarithm = &lr * (q12.log() + q13.log() - q23.log()) / &x2;
    let ordinary = term_13 + term_12 + if_nonzero(&(&r12 * &r13), both_logarithm);
    let neither_both = if_nonzero_else(
        &r12,
        ordinary.clone(),
        if_nonzero_else(&r13, ordinary, both_zero),
    );
    let r14_term = if_nonzero(&r14, dd(&(&qx1 * &q14), &(&qx2 * &q14)) * &r14);
    let complex_mass = one_mass_complex_sheets(
        p,
        mass,
        [&x1, &x2],
        &(&sm / &sr),
        [&w13, &w12],
        [&r34 / &r13, &r24 / &r12, r14],
    );
    let scale = &a * &sr * &sr * &sr * sm;
    // Each branch is scaled on its own: real masses keep their arithmetic.
    let finite = if_nonzero_else(
        &im(mass),
        -complex_mass / &scale,
        -(neither_both - r14_term) / scale,
    );
    if_nonzero(&a, finite)
}

/// The terms of [box_finite_one_mass] for a mass with a nonzero width, on
/// the sheets of Denner and Dittmaier's one-mass box (arXiv:1005.2076,
/// eqs. (3.82) and (3.83); their line 0 is the massive line 4 here), which
/// hold for any real momenta and complex mass.
///
/// Their variables are the roots `u = X sqrt(m²)/|sqrt(m²)|` of the
/// unnormalized quadratic and the ratios `Y23/Y02 = -p2/(m² - p23)`,
/// `Y13/Y01 = -p12/(m² - p4)` and `Y03/m² = (m² - p3)/m²`, each with its
/// principal logarithm: `log(X r34/r13)` is `log u + log(Y23/Y02)`.
/// OneLOop's boxf1 takes `log X + log r34 - log r13` instead: equal for real
/// masses, but the normalization by the complex `sqrt(m²)` turns `X` against
/// `u`, across the negative real axis for a root next to it. In `4 1.0 P
/// 14.78 12.26 -18.32 15.39 -17.48 0 0 0 0 0 0 3.018 -0.302`
/// (tests/d0_one_mass_width.rs) that happened at P = -0.055: OneLOop (both
/// precisions) and these expressions jumped there from -0.0250 to 0.0254 +
/// 1.87i, and at P = -1e-6 gave 6.11 + 383i and 3.03 + 192i for the -0.0599 -
/// 8.2e-4i of COLLIER and of a direct Feynman-parameter integral. A root
/// exactly on the axis takes its lower lip, the limit of the side that keeps
/// the formula.
///
/// The terms are `(F(u1) - F(u2))/(X1 - X2)` with `F(u) = Li2(u, Y23/Y02) +
/// Li2(u, Y13/Y01) - Li2(u, Y03/m²) - log(u) log(Y01 Y02/(Y12 m²))`. When
/// the roots lie on opposite lips of the negative real axis (next to each
/// other near a collinear limit, where the dilogarithm pairs cancel), their
/// principal sheets differ by a winding and the terms grow as `1/(X1 - X2)`
/// and cancel each other (in binary64 a box next to a collinear limit went
/// from 1e-16 to 6e-8 on these sheets).
/// `u2` then takes `u1`'s lip (`F'`), which keeps each pair on one sheet,
/// and the change `F'(u2) - F(u2)` is added back. It is a multiple of `4 pi²`
/// (the dilogarithms' windings cancel against the logarithm's, since
/// `(1 - u v13) (1 - u v12)/(1 - u v14)` is `Y12 m²/(Y01 Y02)` at a root),
/// exact once rounded; its rounding residue would otherwise grow as
/// `1/(X1 - X2)` too.
///
/// With a width, `r12`, `r13` and `r14` cannot vanish, so only the ordinary
/// terms remain. The ratios are `r34/r13`, `r24/r12` and `r14` (the
/// coefficients of the three dilogarithm pairs), the complements the pairs'
/// `1 - X r34/r13` and `1 - X r24/r12`.
fn one_mass_complex_sheets(
    p: [&Atom; 6],
    mass: &Atom,
    [x1, x2]: [&Atom; 2],
    rotation: &Atom,
    [w13, w12]: [&[Atom; 2]; 2],
    [r34_over_r13, r24_over_r12, r14]: [Atom; 3],
) -> Atom {
    let [u1, u2] = [x1, x2].map(|x| x * rotation);
    let upper = |u: &Atom| negative(&-im(u));
    let turns = (upper(&u1) - upper(&u2)) * negative(&re(&u1)) * negative(&re(&u2));
    let first = Q::lower(u1);
    let principal = Q::lower(u2);
    let second = principal.wound(&turns);
    let [v13, v12, v14] = [
        -p[1] / (mass - p[5]),
        -p[4] / (mass - p[3]),
        (mass - p[2]) / mass,
    ]
    .map(Q::lower);
    let pair = |v: &Q| [&first * v, &second * v];
    let [y1, y2] = pair(&v13);
    let channel_13 = complement_difference(&y1, &w13[0], &y2, &w13[1]) * r34_over_r13;
    let [y1, y2] = pair(&v12);
    let channel_12 = complement_difference(&y1, &w12[0], &y2, &w12[1]) * r24_over_r12;
    let invariants = Q::lower((mass - p[3]) / -p[0]).log() + Q::lower((mass - p[5]) / mass).log();
    let logarithm = (&first / &second).log_over_one_minus() * &invariants / x2;
    let [y1, y2] = pair(&v14);
    let one_sheet = channel_13 + channel_12 + logarithm - dd(&y1, &y2) * r14;
    let change =
        |v: &Q, w: &Atom| (&second * v).dilog_complement(w) - (&principal * v).dilog_complement(w);
    let one_minus_14 = 1 - (&principal * &v14).value();
    let windings = change(&v13, &w13[1]) + change(&v12, &w12[1])
        - change(&v14, &one_minus_14)
        - Atom::num(2) * sheet_exact::i() * Symbol::PI.to_atom() * &turns * invariants;
    let four_pi_squared = Atom::num(4) * Symbol::PI.to_atom().pow(2);
    let correction = nearest_integer(&(re(&windings) / &four_pi_squared)) * four_pi_squared;
    one_sheet + if_nonzero(&turns, correction / (x1 - x2))
}

/// The integer nearest to a real `x` with `|x| < 8.5`, from sign tests.
fn nearest_integer(x: &Atom) -> Atom {
    (1..=8).fold(Atom::num(0), |count, j| {
        let half = Atom::num(2 * j - 1) / 2;
        count + (1 - negative(&(x - &half))) - negative(&(x + half))
    })
}

/// The leading coefficient `r34/r24 - r23` of [box_finite_three_masses_zero_third];
/// that representation is degenerate (OneLOop returns zero) when it vanishes.
fn box_three_masses_zero_third_leading(
    p: [&Atom; 6],
    mass_1: &Atom,
    mass_2: &Atom,
    mass_4: &Atom,
) -> Atom {
    let _ = mass_1;
    let sm2 = physical_sqrt(mass_2);
    let sm4 = physical_sqrt(mass_4);
    let sr = square_root_magnitude(mass_2);
    let r23 = (mass_2 - p[1]) / (&sm2 * &sr);
    let k24 = (mass_2 + mass_4 - p[5]) / (&sm2 * &sm4);
    let r34 = (mass_4 - p[2]) / (&sr * &sm4);
    let (r24, _) = r_function(&k24);
    &r34 / &r24 - &r23
}

pub(crate) fn box_finite_three_masses_zero_third(
    p: [&Atom; 6],
    mass_1: &Atom,
    mass_2: &Atom,
    mass_4: &Atom,
) -> Atom {
    let sm1 = physical_sqrt(mass_1);
    let sm2 = physical_sqrt(mass_2);
    let sm4 = physical_sqrt(mass_4);
    let sr = square_root_magnitude(mass_2);
    let k12 = (mass_1 + mass_2 - p[0]) / (&sm1 * &sm2);
    let r13 = (mass_1 - p[4]) / (&sm1 * &sr);
    let k14 = (mass_1 + mass_4 - p[3]) / (&sm1 * &sm4);
    let r23 = (mass_2 - p[1]) / (&sm2 * &sr);
    let k24 = (mass_2 + mass_4 - p[5]) / (&sm2 * &sm4);
    let r34 = (mass_4 - p[2]) / (&sr * &sm4);
    let (r12, _) = r_function(&k12);
    let (r14, _) = r_function(&k14);
    let (r24, d24) = r_function(&k24);
    let a = &r34 / &r24 - &r23;
    let b = -&r13 * &d24 + &k12 * &r34 - &k14 * &r23;
    let c = &k12 * &r13 + &r24 * &r34 - &k14 * &r24 * &r13 - &r23;
    let (x1, x2, _) = solabc_roots(&a, &b, &c, &(&b * &b - Atom::num(4) * &a * &c));
    // Real roots inherit their lips from OneLOop's p -> p(1 + i eps).
    let dk24 = -p[5] / (&sm2 * &sm4);
    let [dr24, dd24, weight] = r_function_velocity(&k24, &r24, &d24, &dk24);
    let dk12 = -p[0] * &weight / (&sm1 * &sm2);
    let dr13 = -p[4] * &weight / (&sm1 * &sr);
    let dk14 = -p[3] * &weight / (&sm1 * &sm4);
    let dr23 = -p[1] * &weight / (&sm2 * &sr);
    let dr34 = -p[2] * weight / (&sr * &sm4);
    let da = &dr34 / &r24 - &r34 * &dr24 / (&r24 * &r24) - &dr23;
    let db =
        -&dr13 * &d24 - &r13 * dd24 + &dk12 * &r34 + &k12 * &dr34 - &dk14 * &r23 - &k14 * &dr23;
    let dc = &dk12 * &r13 + &k12 * &dr13 + &dr24 * &r34 + &r24 * &dr34
        - &dk14 * &r24 * &r13
        - &k14 * &dr24 * &r13
        - &k14 * &r24 * &dr13
        - &dr23;
    let lip = |root: &Atom| root_lip(root, &a, &b, [&da, &db, &dc]);
    let qx1 = Q::with_sign(x1.clone(), lip(&x1));
    let qx2 = Q::with_sign(x2.clone(), lip(&x2));
    let [q12, q13, q14, q23, q24, q34] =
        [&r12, &r13, &r14, &r23, &r24, &r34].map(|r| Q::lower(r.clone()));
    let y1 = &qx1 / &q24;
    let y2 = &qx2 / &q24;
    let mut finite = dd(&(&y1 * &q12), &(&y2 * &q12)) * &r12 / &r24
        + dd(&(&y1 / &q12), &(&y2 / &q12)) / (&r24 * &r12)
        - dd(&(&qx1 * &q14), &(&qx2 * &q14)) * &r14
        - dd(&(&qx1 / &q14), &(&qx2 / &q14)) / &r14;
    let scale23 = &(&q23 / &q13) / &q24;
    let scale34 = &q34 / &q13;
    let nonzero_r13 = -if_nonzero(
        &r23,
        dd(&(&qx1 * &scale23), &(&qx2 * &scale23)) * &r23 / (&r13 * &r24),
    ) + if_nonzero(
        &r34,
        dd(&(&qx1 * &scale34), &(&qx2 * &scale34)) * &r34 / &r13,
    );
    let zero_r13 = -(&qx1 / &qx2).log_over_one_minus() * (q23.log() - q24.log() - q34.log()) / &x2;
    finite += if_nonzero_else(&r13, nonzero_r13, zero_r13);
    if_nonzero(&a, -finite / (&a * sm1 * sm2 * sr * sm4))
}

pub(crate) fn permute_six(values: [&Atom; 6], order: [usize; 6]) -> [&Atom; 6] {
    order.map(|index| values[index])
}

pub(crate) fn box_finite_two_adjacent(p: [&Atom; 6], mass_3: &Atom, mass_4: &Atom) -> Atom {
    let sm3 = physical_sqrt(mass_3);
    let sm4 = physical_sqrt(mass_4);
    let sr = square_root_magnitude(mass_3);
    let r12 = (mass_4 - p[3]) / (&sr * &sm4);
    let r13 = (mass_4 - p[5]) / (&sr * &sm4);
    let k14 = (mass_3 + mass_4 - p[2]) / (&sm3 * &sm4);
    let r23 = -p[0] / (&sr * &sr);
    let r24 = (mass_3 - p[4]) / (&sr * &sm3);
    let r34 = (mass_3 - p[1]) / (&sr * &sm3);
    let (r14, _) = r_function(&k14);
    let a = &r34 * &r24 - &r23;
    let b = &r13 * &r24 + &r12 * &r34 - &k14 * &r23;
    let c = &r12 * &r13 - &r23;
    let gap = channel_gap([&r12, &r13, &r24, &r34]);
    let square = collinear_square(
        &gap,
        &r13 * &r24 + &r12 * &r34,
        &k14 * &r23,
        Atom::num(4) * &r23 * (&r24 * &r34 + &r12 * &r13 - &r23),
    );
    let (x1, x2, discriminant) = solabc_roots(&a, &b, &c, &square);
    // Real roots inherit their lips from OneLOop's p -> p(1 + i eps).
    let [d12, d13] = [p[3], p[5]].map(|v| -v / (&sr * &sm4));
    let dk14 = -p[2] / (&sm3 * &sm4);
    let d23 = -p[0] / (&sr * &sr);
    let [d24, d34] = [p[4], p[1]].map(|v| -v / (&sr * &sm3));
    let da = &d34 * &r24 + &r34 * &d24 - &d23;
    let db = &d13 * &r24 + &r13 * &d24 + &d12 * &r34 + &r12 * &d34 - &dk14 * &r23 - &k14 * &d23;
    let dc = &d12 * &r13 + &r12 * &d13 - &d23;
    let lip = |root: &Atom| root_lip(root, &a, &b, [&da, &db, &dc]);
    let qx1 = ExactSheetAtom::with_sign(x1.clone(), lip(&x1));
    let qx2 = ExactSheetAtom::with_sign(x2.clone(), lip(&x2));
    let q12 = ExactSheetAtom::lower(r12.clone());
    let q13 = ExactSheetAtom::lower(r13.clone());
    let q14 = ExactSheetAtom::lower(r14.clone());
    let q23 = ExactSheetAtom::lower(r23.clone());
    let q24 = ExactSheetAtom::lower(r24.clone());
    let q34 = ExactSheetAtom::lower(r34.clone());
    let x1_14 = &qx1 * &q14;
    let x2_14 = &qx2 * &q14;
    let x1_i14 = &qx1 / &q14;
    let x2_i14 = &qx2 / &q14;
    let scale_13 = &q34 / &q13;
    let scale_12 = &q24 / &q12;
    let x1_13 = &qx1 * &scale_13;
    let x2_13 = &qx2 * &scale_13;
    let x1_12 = &qx1 * &scale_12;
    let x2_12 = &qx2 * &scale_12;
    let ratio = &qx1 / &qx2;
    let q12q13 = &q12 * &q13;
    let log_scale = &q12q13 / &q23;
    let lr = ratio.log_over_one_minus();
    // Continued logarithms add exactly, without expanding phase bookkeeping
    // for products that appear only inside a logarithm.
    let root_log = qx1.log() + qx2.log();
    let limit_13 = &lr * (&root_log / 2 + q34.log() + q12.log() - q23.log()) / &x2;
    let limit_12 = &lr * (&root_log / 2 + q24.log() + q13.log() - q23.log()) / &x2;
    let both_zero = &lr * (root_log + q34.log() + q24.log() - q23.log()) / &x2;
    // At u = r13/r34 and r12/r24 the quadratic is r23 (k14 u - u² - 1).
    let quadratic = [&a, &b, &discriminant];
    let complements = |u: Atom, shifted_b: Atom| {
        let shifted_b = shifted_b - &k14 * &r23 + Atom::num(2) * &r23 * &u;
        let shifted_c = &r23 * (&k14 * &u - &u * &u - 1);
        root_complements(quadratic, &(Atom::num(1) / &u), &shifted_b, &shifted_c)
    };
    let [w1, w2] = complements(&r13 / &r34, -&gap);
    let term_13 = if_nonzero_else(
        &r13,
        if_nonzero(&r34, complement_difference(&x1_13, &w1, &x2_13, &w2) * &r34 / &r13),
        limit_13,
    );
    let [w1, w2] = complements(&r12 / &r24, gap.clone());
    let term_12 = if_nonzero_else(
        &r12,
        if_nonzero(&r24, complement_difference(&x1_12, &w1, &x2_12, &w2) * &r24 / &r12),
        limit_12,
    );
    let ordinary = term_13 + term_12 + if_nonzero(&(&r12 * &r13), &lr * log_scale.log() / &x2);
    let finite = -exact_sheet_difference(&x1_14, &x2_14) * &r14
        - exact_sheet_difference(&x1_i14, &x2_i14) / &r14
        + if_nonzero_else(
            &r12,
            ordinary.clone(),
            if_nonzero_else(&r13, ordinary, both_zero),
        );
    if_nonzero(&a, -finite / (&a * &sr * &sr * sm3 * sm4))
}

pub(crate) fn box_finite_two_opposite(p: [&Atom; 6], mass_2: &Atom, mass_4: &Atom) -> Atom {
    box_finite_two_adjacent(permute_six(p, [4, 1, 5, 3, 0, 2]), mass_2, mass_4)
}

pub(crate) fn with_box_normalization(series: LaurentSeries) -> LaurentSeries {
    let [finite, pole, double] = series.into_coefficients();
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    LaurentSeries::new(finite + &double * pi_squared / 12, pole, double)
}

pub(crate) fn box_one_mass_06(
    p12: &Atom,
    p23: &Atom,
    mass: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let r13 = -p12;
    let r24 = mass - p23;
    // The quotient inherits both Feynman lips. Applying a lower-lip log
    // directly to the quotient loses the continuation for negative masses.
    let q13 = Q::lower(r13.clone());
    let q24 = Q::lower(r24.clone());
    let qm = Q::lower(mass.clone());
    let log_mass = qm.scale_positive(&(Atom::num(1) / mu_squared)).log();
    let log_1 = (&q13 / &qm).log();
    let log_2 = (&q24 / &qm).log();
    let z2 = Atom::num(2);
    let z1 = -Atom::num(2) * &log_2 - &log_1;
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    let z0 = Atom::num(2) * (&log_2 * &log_1 - pi_squared / 3);
    let factor = Atom::num(1) / (&r13 * &r24);
    with_box_normalization(LaurentSeries::new(
        &factor * (z0 + (&z2 * &log_mass / 2 - &z1) * &log_mass),
        &factor * (z1 - &z2 * &log_mass),
        factor * z2,
    ))
}

// The one-mass infrared sectors follow OneLOop's box07-box10 in its continued
// (qmplx) arithmetic: every invariant carries its own Feynman lip, and ratios
// and products keep the accumulated phase. A principal logarithm or
// dilogarithm of a plain ratio loses that phase when invariants of opposite
// signs meet, e.g. for timelike channels with spacelike external masses.
pub(crate) fn box_one_mass_07(
    p4: &Atom,
    p12: &Atom,
    p23: &Atom,
    mass: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let r13 = -p12;
    let r14 = mass - p4;
    let r24 = mass - p23;
    let [q13, q14, q24, qm] = [&r13, &r14, &r24, mass].map(|r| Q::lower(r.clone()));
    let inverse_mu = Atom::num(1) / mu_squared;
    let log_mass = qm.scale_positive(&inverse_mu).log();
    let log_12 = (&q13 / &qm).log();
    let log_23 = (&q24 / &qm).log();
    let log_4 = (&q14 / &qm).log();
    let z2 = Atom::num(3) / 2;
    let z1 = -Atom::num(2) * &log_23 - &log_12 + &log_4;
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    let z0 = Atom::num(2) * (&log_12 * &log_23 - (&q14 / &q24).dilog())
        - &log_4 * &log_4
        - Atom::num(13) * pi_squared / 24;
    let factor = Atom::num(1) / (&r13 * &r24);
    with_box_normalization(LaurentSeries::new(
        &factor * (z0 + (&z2 * &log_mass / 2 - &z1) * &log_mass),
        &factor * (z1 - &z2 * &log_mass),
        factor * z2,
    ))
}

pub(crate) fn box_one_mass_08(
    p3: &Atom,
    p4: &Atom,
    p12: &Atom,
    p23: &Atom,
    mass: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let r13 = -p12;
    let r34 = mass - p3;
    let r14 = mass - p4;
    let r24 = mass - p23;
    let [q13, q34, q14, q24, qm] = [&r13, &r34, &r14, &r24, mass].map(|r| Q::lower(r.clone()));
    let inverse_mu = Atom::num(1) / mu_squared;
    let x1 = &q34 / &q24;
    let x2 = &q14 / &q24;
    let x3 = q13.scale_positive(&inverse_mu);
    let z1 = (&(&x1 * &x2) / &x3).log();
    let mut z0 =
        Atom::num(2) * (q24.scale_positive(&inverse_mu).log() * x3.log() - x1.dilog() - x2.dilog());
    let scaled_1 = q34.scale_positive(&inverse_mu);
    let scaled_2 = q14.scale_positive(&inverse_mu);
    let xx = &(&scaled_1 * &scaled_2) / &x3;
    let log_1 = scaled_1.log();
    let log_2 = scaled_2.log();
    let log_xx = xx.log();
    z0 += -&log_1 * &log_1 - &log_2 * &log_2
        + &log_xx * &log_xx / 2
        + (&qm / &xx).scale_positive(&inverse_mu).dilog();
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    let factor = Atom::num(1) / (&r13 * &r24);
    with_box_normalization(LaurentSeries::new(
        &factor * (z0 - pi_squared / 4),
        &factor * z1,
        factor,
    ))
}

pub(crate) fn box_one_mass_09(
    p2: &Atom,
    p3: &Atom,
    p12: &Atom,
    p23: &Atom,
    mass: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let r23 = -p2;
    let r13 = -p12;
    let r34 = mass - p3;
    let r24 = mass - p23;
    let [q23, q13, q34, q24, qm] = [&r23, &r13, &r34, &r24, mass].map(|r| Q::lower(r.clone()));
    let log_mass = qm.scale_positive(&(Atom::num(1) / mu_squared)).log();
    let q12 = &q13 / &q23;
    let q23_mass = &q24 / &qm;
    let z2 = Atom::num(1) / 2;
    let z1 = -q12.log() - q23_mass.log();
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    let z0 = (&(&q23_mass * &q34) / &q23).dilog()
        + Atom::num(2) * q12.dilog()
        + &z1 * &z1
        + pi_squared / 24;
    let factor = Atom::num(1) / (&r13 * &r24);
    with_box_normalization(LaurentSeries::new(
        &factor * (z0 + (&z2 * &log_mass / 2 - &z1) * &log_mass),
        &factor * (z1 - &z2 * &log_mass),
        factor * z2,
    ))
}

pub(crate) fn box_one_mass_10(
    p2: &Atom,
    p3: &Atom,
    p4: &Atom,
    p12: &Atom,
    p23: &Atom,
    mass: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    // OneLOop's box10 exchanges (p2, p4) with (p12, p23), a relabeling of the
    // two massless propagators, when that enlarges |p12 (m - p23)|.
    let first = p12 * (mass - p23);
    let second = p2 * (mass - p4);
    let swap = negative(&(Symbol::ABS.call((&first,)) - Symbol::ABS.call((&second,))));
    let original = [p2.clone(), p4.clone(), p12.clone(), p23.clone()];
    let [p2, p4, p12, p23] = [2, 3, 0, 1]
        .map(|k| if_nonzero_else(&swap, original[k].clone(), original[(k + 2) % 4].clone()));
    let r23 = -&p2;
    let r13 = -&p12;
    let r34 = mass - p3;
    let r14 = mass - &p4;
    let r24 = mass - &p23;
    let [q23, q13, q34, q14, q24, qm] =
        [&r23, &r13, &r34, &r14, &r24, mass].map(|r| Q::lower(r.clone()));
    let base = &q34 / &qm;
    let extra = -dd(&(&(&base * &q24) / &q23), &(&(&base * &q14) / &q13)) * &r34
        / (Atom::num(2) * mass * &r23);
    let x1 = &q23 / &q13;
    let x2 = &q24 / &q14;
    let xx = &x1 / &x2;
    let z1 = -xx.log_over_one_minus() / &r24;
    let z0 = if_nonzero(&r34, extra) - dd(&x1, &x2) / &r14
        + dd(&xx, &Q::lower(Atom::num(1))) / &r24
        + &z1 * ((&qm / &q24).log() - qm.scale_positive(&(Atom::num(1) / mu_squared)).log() / 2);
    with_box_normalization(LaurentSeries::new(
        -Atom::num(2) * z0 / &r13,
        -z1 / r13,
        Atom::num(0),
    ))
}

pub(crate) fn box_one_mass_ir(p: [&Atom; 6], mass: &Atom, mu_squared: &Atom) -> LaurentSeries {
    let d3 = p[2] - mass;
    let d4 = p[3] - mass;
    let case_06 = box_one_mass_06(p[4], p[5], mass, mu_squared);
    let case_07_d4 = box_one_mass_07(p[3], p[4], p[5], mass, mu_squared);
    let case_07_d3 = box_one_mass_07(p[2], p[4], p[5], mass, mu_squared);
    let case_08 = box_one_mass_08(p[2], p[3], p[4], p[5], mass, mu_squared);
    let case_09_p2 = box_one_mass_09(p[1], p[2], p[4], p[5], mass, mu_squared);
    let case_10_p2 = box_one_mass_10(p[1], p[2], p[3], p[4], p[5], mass, mu_squared);
    let case_09_p1 = box_one_mass_09(p[0], p[3], p[4], p[5], mass, mu_squared);
    let case_10_p1 = box_one_mass_10(p[0], p[3], p[2], p[4], p[5], mass, mu_squared);
    if_series(
        p[0],
        if_series(&d3, case_10_p1, case_09_p1),
        if_series(
            p[1],
            if_series(&d4, case_10_p2, case_09_p2),
            if_series(
                &d3,
                if_series(&d4, case_08, case_07_d3),
                if_series(&d4, case_07_d4, case_06),
            ),
        ),
    )
}

pub(crate) fn box_one_mass(p: [&Atom; 6], mass: &Atom, mu_squared: &Atom) -> LaurentSeries {
    let infrared = box_one_mass_ir(p, mass, mu_squared);
    if_series(
        p[0],
        if_series(
            p[1],
            finite_series(box_finite_one_mass(p, mass)),
            infrared.clone(),
        ),
        infrared,
    )
}

pub(crate) fn box_two_adjacent_11(
    p3: &Atom,
    p12: &Atom,
    p23: &Atom,
    mass_3: &Atom,
    mass_4: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let sqrt_3 = physical_sqrt(mass_3);
    let sqrt_4 = physical_sqrt(mass_4);
    let mu = mu_squared.sqrt();
    let r13 = (mass_3 - p12) / (&mu * &sqrt_3);
    let r24 = (mass_4 - p23) / (&mu * &sqrt_4);
    let (r34, _) = r_function(&((mass_3 + mass_4 - p3) / (&sqrt_3 * &sqrt_4)));
    let l13 = physical_log(&r13);
    let l24 = physical_log(&r24);
    let l34 = physical_log(&r34);
    let denominator = (mass_3 - p12) * (mass_4 - p23);
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    with_box_normalization(LaurentSeries::new(
        (Atom::num(2) * &l13 * &l24 - &l34 * &l34 - Atom::num(7) * pi_squared / 12) / &denominator,
        (-l13 - l24) / &denominator,
        Atom::num(1) / denominator,
    ))
}

pub(crate) fn box_two_adjacent_12(
    p3: &Atom,
    p4: &Atom,
    p12: &Atom,
    p23: &Atom,
    mass_3: &Atom,
    mass_4: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let sqrt_3 = physical_sqrt(mass_3);
    let sqrt_4 = physical_sqrt(mass_4);
    let mu = mu_squared.sqrt();
    let r13 = (mass_3 - p12) / (&mu * &sqrt_3);
    let r14 = (mass_4 - p4) / (&mu * &sqrt_4);
    let r24 = (mass_4 - p23) / (&mu * &sqrt_4);
    let (r34, _) = r_function(&((mass_3 + mass_4 - p3) / (&sqrt_3 * &sqrt_4)));
    let [q13, q14, q24, q34] = [&r13, &r14, &r24, &r34].map(|r| Q::lower(r.clone()));
    let l13 = q13.log();
    let l14 = q14.log();
    let l24 = q24.log();
    let l34 = q34.log();
    let y = &q14 / &q13;
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    let finite = Atom::num(2) * &l13 * &l24
        - &l14 * &l14
        - &l34 * &l34
        - Atom::num(2) * (&q14 / &q24).dilog()
        - (&y * &q34).dilog()
        - (&y / &q34).dilog()
        - pi_squared / 8;
    let denominator = (mass_3 - p12) * (mass_4 - p23);
    with_box_normalization(LaurentSeries::new(
        finite / &denominator,
        (l14 - l24 - l13) / &denominator,
        Atom::num(1) / (Atom::num(2) * denominator),
    ))
}

#[allow(clippy::too_many_arguments)] // OneLOop sector 13's named invariants.
pub(crate) fn box_two_adjacent_13(
    p2: &Atom,
    p3: &Atom,
    p4: &Atom,
    p12: &Atom,
    p23: &Atom,
    mass_3: &Atom,
    mass_4: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    // Match box13's invariant ordering before constructing the continued
    // ratios; products and quotients below retain their individual sheets.
    let first = (mass_3 - p12) * (mass_4 - p23);
    let second = (mass_3 - p2) * (mass_4 - p4);
    let swap = negative(&(Symbol::ABS.call((&first,)) - Symbol::ABS.call((&second,))));
    let original = [p2.clone(), p4.clone(), p12.clone(), p23.clone()];
    let [p2, p4, p12, p23] = [2, 3, 0, 1]
        .map(|k| if_nonzero_else(&swap, original[k].clone(), original[(k + 2) % 4].clone()));
    let sqrt_3 = physical_sqrt(mass_3);
    let sqrt_4 = physical_sqrt(mass_4);
    let mu = mu_squared.sqrt();
    let r13 = (mass_3 - &p12) / (&mu * &sqrt_3);
    let r14 = (mass_4 - &p4) / (&mu * &sqrt_4);
    let r23 = (mass_3 - &p2) / (&mu * &sqrt_3);
    let r24 = (mass_4 - &p23) / (&mu * &sqrt_4);
    let (r34, _) = r_function(&((mass_3 + mass_4 - p3) / (&sqrt_3 * &sqrt_4)));
    let [q13, q14, q23, q24, q34] = [&r13, &r14, &r23, &r24, &r34].map(|r| Q::lower(r.clone()));
    let y = &(&(&q14 * &q23) / &q13) / &q24;
    let scale = &r13 * &r24;
    let logd = y.log_over_one_minus() / &scale;
    let li2d = dd(&y, &Q::lower(Atom::num(1))) / &scale;
    let y1 = &q23 / &q24;
    let y2 = &q13 / &q14;
    let li2f = dd(&(&y1 * &q34), &(&y2 * &q34)) * &r34 / (&r14 * &r24);
    let li2b = dd(&(&y1 / &q34), &(&y2 / &q34)) / (&r34 * &r14 * &r24);
    let li2e = dd(&(&q14 / &q24), &(&q13 / &q23)) / (&r23 * &r24);
    let denominator = &mu * &mu * sqrt_3 * sqrt_4;
    with_box_normalization(LaurentSeries::new(
        (li2f + li2b + Atom::num(2) * li2e
            - Atom::num(2) * li2d
            - Atom::num(2) * &logd * q13.log())
            / &denominator,
        logd / denominator,
        Atom::num(0),
    ))
}

pub(crate) fn box_two_adjacent_ir(
    p: [&Atom; 6],
    mass_3: &Atom,
    mass_4: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let d2 = p[1] - mass_3;
    let d4 = p[3] - mass_4;
    let case_11 = box_two_adjacent_11(p[2], p[4], p[5], mass_3, mass_4, mu_squared);
    let case_12_d4 = box_two_adjacent_12(p[2], p[3], p[4], p[5], mass_3, mass_4, mu_squared);
    let case_12_d2 = box_two_adjacent_12(p[2], p[1], p[5], p[4], mass_4, mass_3, mu_squared);
    let case_13 = box_two_adjacent_13(p[1], p[2], p[3], p[4], p[5], mass_3, mass_4, mu_squared);
    let infrared = if_series(
        &d2,
        if_series(&d4, case_13, case_12_d2),
        if_series(&d4, case_12_d4, case_11),
    );
    if_series(
        p[0],
        finite_series(box_finite_two_adjacent(p, mass_3, mass_4)),
        infrared,
    )
}

pub(crate) fn box_two_opposite_14(
    p12: &Atom,
    p23: &Atom,
    mass_2: &Atom,
    mass_4: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let sqrt_2 = physical_sqrt(mass_2);
    let sqrt_4 = physical_sqrt(mass_4);
    let (r24, _) = r_function(&((mass_2 + mass_4 - p23) / (&sqrt_2 * &sqrt_4)));
    let coefficient = -Atom::num(2) * log_over_one_minus(&r24) * &r24
        / ((Atom::num(1) + &r24) * sqrt_2 * sqrt_4 * p12);
    with_box_normalization(LaurentSeries::new(
        -&coefficient * physical_log(&(-p12 / mu_squared)),
        coefficient,
        Atom::num(0),
    ))
}

pub(crate) fn box_two_opposite_15(
    p2: &Atom,
    p3: &Atom,
    p12: &Atom,
    p23: &Atom,
    mass_2: &Atom,
    mass_4: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    // Match OneLOop's ``abs(m2-p2) > abs(m4-p3)`` ordering before forming
    // the qmplx variables.  Keeping the choice symbolic preserves the
    // ordering when this expression is evaluated with substituted inputs.
    let left = mass_2 - p2;
    let right = mass_4 - p3;
    let swap = negative(&(&left * &left.conj() - &right * &right.conj()));
    let orig_p2 = p2.clone();
    let orig_p3 = p3.clone();
    let orig_m2 = mass_2.clone();
    let orig_m4 = mass_4.clone();
    let p2 = if_nonzero_else(&swap, orig_p3.clone(), orig_p2.clone());
    let p3 = if_nonzero_else(&swap, orig_p2.clone(), orig_p3.clone());
    let mass_2 = if_nonzero_else(&swap, orig_m4.clone(), orig_m2.clone());
    let mass_4 = if_nonzero_else(&swap, orig_m2, orig_m4);
    let mu = mu_squared.sqrt();
    let sqrt_2 = physical_sqrt(&mass_2);
    let sqrt_4 = physical_sqrt(&mass_4);
    let sqrt_real = square_root_magnitude(&mass_2);
    let r13 = -p12 / (&mu * &sqrt_real);
    let r23 = (mass_2.clone() - p2) / (&sqrt_2 * &sqrt_real);
    let r34 = (mass_4.clone() - p3) / (&sqrt_real * &sqrt_4);
    let (r24, _) = r_function(&((mass_2 + mass_4 - p23) / (&sqrt_2 * &sqrt_4)));
    let [q13, q23, q24, q34] = [&r13, &r23, &r24, &r34].map(|r| Q::lower(r.clone()));
    let qratio = &q13 / &q23;
    let qsquare = &qratio * &qratio;
    let qss = &qsquare / &q24;
    let factor = &r24 / (&sqrt_2 * &sqrt_4 * p12);
    let log24 = q24.log_over_one_minus() / (Atom::num(1) + &r24);
    let qz = &q34 / &q23;
    let qz1 = &qz * &q24;
    let qz2 = &qz / &q24;
    let correction = dd(&qz1, &qz2) * &r34 / (&r23 * &r24);
    let finite = &log24 * qss.log() + dd(&(&q24 * &q24), &Q::lower(Atom::num(1)))
        - if_nonzero(&r34, correction);
    with_box_normalization(LaurentSeries::new(
        &factor * finite,
        -factor * log24,
        Atom::num(0),
    ))
}

pub(crate) fn box_two_opposite_ir(
    p: [&Atom; 6],
    mass_2: &Atom,
    mass_4: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let d1 = p[0] - mass_2;
    let d4 = p[3] - mass_4;
    let d2 = p[1] - mass_2;
    let d3 = p[2] - mass_4;
    let finite = finite_series(box_finite_two_opposite(p, mass_2, mass_4));
    let case_15_first = box_two_opposite_15(p[0], p[3], p[4], p[5], mass_2, mass_4, mu_squared);
    let case_15_second = box_two_opposite_15(p[1], p[2], p[4], p[5], mass_2, mass_4, mu_squared);
    let case_14 = box_two_opposite_14(p[4], p[5], mass_2, mass_4, mu_squared);
    let first_true = if_series(
        &d2,
        finite.clone(),
        if_series(&d3, finite.clone(), case_15_first),
    );
    let first_false = if_series(
        &d2,
        case_15_second.clone(),
        if_series(&d3, case_15_second, case_14),
    );
    if_series(
        &d1,
        first_true.clone(),
        if_series(&d4, first_true, first_false),
    )
}

#[allow(clippy::too_many_arguments)] // OneLOop sector 16's named invariants.
pub(crate) fn box_three_masses_16(
    p2: &Atom,
    p3: &Atom,
    p12: &Atom,
    p23: &Atom,
    mass_2: &Atom,
    mass_3: &Atom,
    mass_4: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    // The reference implementation orders the two massive endpoints by
    // magnitude, exchanging the adjacent invariants at the same time.
    let swap = negative(&(mass_2 * &mass_2.conj() - mass_4 * &mass_4.conj()));
    let orig_p2 = p2.clone();
    let orig_p3 = p3.clone();
    let orig_m2 = mass_2.clone();
    let orig_m4 = mass_4.clone();
    let p2 = if_nonzero_else(&swap, orig_p3.clone(), orig_p2.clone());
    let p3 = if_nonzero_else(&swap, orig_p2, orig_p3);
    let mass_2 = if_nonzero_else(&swap, orig_m4.clone(), orig_m2.clone());
    let mass_4 = if_nonzero_else(&swap, orig_m2, orig_m4);
    let mu = mu_squared.sqrt();
    let sqrt_2 = physical_sqrt(&mass_2);
    let sqrt_3 = physical_sqrt(mass_3);
    let sqrt_4 = physical_sqrt(&mass_4);
    let r13 = (mass_3 - p12) / (&mu * &sqrt_3);
    let (r23, _) = r_function(&((mass_2.clone() + mass_3 - p2) / (&sqrt_2 * &sqrt_3)));
    let (r24, _) = r_function(&((mass_2 + mass_4.clone() - p23) / (&sqrt_2 * &sqrt_4)));
    let r34_arg = (mass_3 + mass_4) - p3;
    let (r34, _) = r_function(&(&r34_arg / (&sqrt_3 * &sqrt_4)));
    let [q13, q23, q24, q34] = [&r13, &r23, &r24, &r34].map(|r| Q::lower(r.clone()));
    let qproduct = &q13 * &q23;
    let qsquare = &qproduct * &qproduct;
    let qss = &qsquare / &q24;
    let y1 = &(&q23 * &q34) * &q24;
    let y2 = &(&q23 * &q34) / &q24;
    let z1 = &(&q23 / &q34) * &q24;
    let z2 = &(&q23 / &q34) / &q24;
    let factor = Atom::num(1) / (&sqrt_2 * &sqrt_4 * (p12 - mass_3));
    let log24 = q24.log_over_one_minus() * &r24 / (Atom::num(1) + &r24);
    let finite = &log24 * qss.log() + dd(&(&q24 * &q24), &Q::lower(Atom::num(1))) * &r24
        - dd(&y1, &y2) * &r23 * &r34
        - dd(&z1, &z2) * &r23 / &r34;
    with_box_normalization(LaurentSeries::new(
        &factor * finite,
        -factor * log24,
        Atom::num(0),
    ))
}

pub(crate) fn box_three_masses_ir(
    p: [&Atom; 6],
    masses: [&Atom; 4],
    mu_squared: &Atom,
) -> LaurentSeries {
    let d1 = p[0] - masses[1];
    let d4 = p[3] - masses[3];
    // OneLOop's boxf3 relabels the propagators (m1 = 0, m2, m3, m4) as
    // (m2, m3, m1, m4) before applying boxf33. Other relabelings are equally
    // valid for the integral, but not for boxf33's continuation of its real
    // roots and r-functions: a cyclic relabeling fails in some regions.
    // When that chart degenerates (its leading coefficient vanishes, e.g. for
    // a vanishing channel where the contour form is unavailable), use the
    // cyclic relabeling (m3, m4, m1, m2) instead.
    let finite_p = permute_six(p, [1, 4, 3, 5, 0, 2]);
    let leading = box_three_masses_zero_third_leading(finite_p, masses[1], masses[2], masses[3]);
    let cyclic_p = permute_six(p, [2, 3, 0, 1, 4, 5]);
    let ordinary = if_nonzero_else(
        &leading,
        box_finite_three_masses_zero_third(finite_p, masses[1], masses[2], masses[3]),
        box_finite_three_masses_zero_third(cyclic_p, masses[2], masses[3], masses[1]),
    );
    let finite = finite_series(finite_three_mass_representation(p, masses, ordinary));
    let infrared = box_three_masses_16(
        p[1], p[2], p[4], p[5], masses[1], masses[2], masses[3], mu_squared,
    );
    if_series(&d1, finite.clone(), if_series(&d4, finite, infrared))
}

// OneLOop's box03 and box05 in continued (qmplx) arithmetic. Each -p carries
// its own lower lip; quotients and products keep the accumulated phase.
fn massless_invariants<const N: usize>(p: [&Atom; N]) -> [ExactSheetAtom; N] {
    p.map(|v| Q::lower(-v))
}

pub(crate) fn box_massless_two_opposite(p: [&Atom; 6], mu_squared: &Atom) -> LaurentSeries {
    let [q2, q4, q5, q6] = massless_invariants([p[1], p[3], p[4], p[5]]);
    let q26 = &q2 / &q6;
    let q54 = &q5 / &q4;
    let y = &q26 / &q54;
    let log_y = y.log_over_one_minus() / (p[4] * p[5]);
    let scale = &(&q54 * &q2) * &q6;
    let finite = Atom::num(2)
        * (dd(&(&q6 / &q4), &(&q2 / &q5)) / (p[3] * p[4]) + dd(&q54, &q26) / (p[3] * p[5])
            - dd(&y, &Q::lower(Atom::num(1))) / (p[4] * p[5])
            - &log_y
                * scale
                    .scale_positive(&(Atom::num(1) / (mu_squared * mu_squared)))
                    .log()
                / 2);
    with_box_normalization(LaurentSeries::new(
        finite,
        Atom::num(2) * log_y,
        Atom::num(0),
    ))
}

pub(crate) fn box_massless_three_nonzero(p: [&Atom; 6], mu_squared: &Atom) -> LaurentSeries {
    let [q2, q3, q4, q5, q6] = massless_invariants([p[1], p[2], p[3], p[4], p[5]]);
    let q25 = &q2 / &q5;
    let q64 = &q6 / &q4;
    let y = &q25 / &q64;
    let log_y = y.log_over_one_minus() / (p[4] * p[5]);
    let z = &(&(&(&(&q64 * &q2) * &q5) * &q6) * &q6) / &(&q3 * &q3);
    let finite = Atom::num(2)
        * (dd(&q64, &q25) / (p[3] * p[4])
            - dd(&y, &Q::lower(Atom::num(1))) / (p[4] * p[5])
            - &log_y
                * z.scale_positive(&(Atom::num(1) / (mu_squared * mu_squared)))
                    .log()
                / 4);
    with_box_normalization(LaurentSeries::new(finite, log_y, Atom::num(0)))
}

pub(crate) fn box_massless_low_case(
    p: [&Atom; 6],
    case: usize,
    mu_squared: &Atom,
) -> LaurentSeries {
    if case == 5 {
        return box_massless_two_opposite(p, mu_squared);
    }
    if case == 3 {
        return box_massless_three_nonzero(p, mu_squared);
    }
    if case == 4 {
        return finite_series(box_finite_massless(p));
    }

    let factor = Atom::num(1) / (p[4] * p[5]);
    // OneLOop's box00 continues Li2(1 - x f), f = (p5 + p6 - p4)/(p5 p6),
    // for negative real x f on the lip given by the sign of f.
    let scaled_dilogs = |factor: &Atom| {
        let scaled = factor * (p[4] + p[5] - p[3]);
        let lip = sheet_exact::sign_nonnegative(&scaled);
        [p[3], p[4], p[5]].map(|v| Q::with_sign(v * &scaled, lip.clone()).dilog())
    };
    let l3 = physical_log(&(-p[2] / mu_squared));
    let l4 = physical_log(&(-p[3] / mu_squared));
    let l5 = physical_log(&(-p[4] / mu_squared));
    let l6 = physical_log(&(-p[5] / mu_squared));
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    let (finite, pole, double) = match case {
        0 => {
            // Preserve the two Feynman lips for mixed-sign real invariants.
            let ratio = &l5 - &l6;
            (
                &factor
                    * (&l5 * &l5 + &l6 * &l6 - &ratio * &ratio - Atom::num(4) * &pi_squared / 3),
                -Atom::num(2) * &factor * (&l5 + &l6),
                Atom::num(4) * &factor,
            )
        }
        1 => {
            let [li4, li5, li6] = scaled_dilogs(&factor);
            (
                &factor
                    * (&l5 * &l5 + &l6 * &l6 - &l4 * &l4 - &pi_squared / 2
                        + Atom::num(2) * (li5 + li6 - li4)),
                Atom::num(2) * &factor * (&l4 - &l5 - &l6),
                Atom::num(2) * &factor,
            )
        }
        2 => {
            let [li4, li5, li6] = scaled_dilogs(&factor);
            let [q3, q4, q5, q6] = massless_invariants([p[2], p[3], p[4], p[5]]);
            let li54 = (&q4 / &q5).dilog();
            let li63 = (&q3 / &q6).dilog();
            (
                &factor
                    * (&l5 * &l5 + &l6 * &l6 - &l3 * &l3 - &l4 * &l4
                        + (&l3 + &l4 - &l5) * (&l3 + &l4 - &l5) / 2
                        - &pi_squared / 12
                        + Atom::num(2) * (li54 - li63 + li5 + li6 - li4)),
                &factor * (&l4 + &l3 - &l5 - Atom::num(2) * &l6),
                factor,
            )
        }
        _ => (Atom::num(0), Atom::num(0), Atom::num(0)),
    };
    with_box_normalization(LaurentSeries::new(finite, pole, double))
}

pub(crate) fn select_four_series(
    conditions: [&Atom; 4],
    branches: [LaurentSeries; 16],
) -> LaurentSeries {
    let coefficients = branches.map(LaurentSeries::into_coefficients);
    let select = |coefficient: usize| {
        let pair = |offset: usize| {
            if_nonzero_else(
                conditions[3],
                coefficients[offset + 1][coefficient].clone(),
                coefficients[offset][coefficient].clone(),
            )
        };
        if_nonzero_else(
            conditions[0],
            if_nonzero_else(
                conditions[1],
                if_nonzero_else(conditions[2], pair(14), pair(12)),
                if_nonzero_else(conditions[2], pair(10), pair(8)),
            ),
            if_nonzero_else(
                conditions[1],
                if_nonzero_else(conditions[2], pair(6), pair(4)),
                if_nonzero_else(conditions[2], pair(2), pair(0)),
            ),
        )
    };
    LaurentSeries::new(select(0), select(1), select(2))
}

pub(crate) fn box_massless(p: [&Atom; 6], mu_squared: &Atom) -> LaurentSeries {
    const PERMUTATIONS: [[usize; 6]; 16] = [
        [0, 1, 2, 3, 4, 5],
        [0, 1, 2, 3, 4, 5],
        [3, 0, 1, 2, 5, 4],
        [0, 1, 2, 3, 4, 5],
        [2, 3, 0, 1, 4, 5],
        [0, 1, 2, 3, 4, 5],
        [3, 0, 1, 2, 5, 4],
        [0, 1, 2, 3, 4, 5],
        [1, 2, 3, 0, 5, 4],
        [1, 2, 3, 0, 5, 4],
        [3, 0, 1, 2, 5, 4],
        [1, 2, 3, 0, 5, 4],
        [2, 3, 0, 1, 4, 5],
        [2, 3, 0, 1, 4, 5],
        [3, 0, 1, 2, 5, 4],
        [0, 1, 2, 3, 4, 5],
    ];
    const CASES: [usize; 16] = [0, 1, 1, 2, 1, 5, 2, 3, 1, 2, 5, 3, 2, 3, 3, 4];
    let branches = core::array::from_fn(|index| {
        box_massless_low_case(
            permute_six(p, PERMUTATIONS[index]),
            CASES[index],
            mu_squared,
        )
    });
    select_four_series([p[0], p[1], p[2], p[3]], branches)
}

/// Constructs a scalar box with its transparent native function definitions.
/// Momenta are [p1²,p2²,p3²,p4²,(p1+p2)²,(p2+p3)²].
///
/// The equal-mass vacuum limit and native boxc reductions are implemented, but
/// arbitrary Gram/Cayley degeneracies and threshold limits are not fully validated.
/// See STATUS.md before using degenerate kinematics.
pub fn d0(p: [&Atom; 6], m: [&Atom; 4], mu_squared: &Atom) -> MappedLaurentSeries {
    let expressions = OneLoopExpressions::new();
    MappedLaurentSeries {
        series: expressions.d0(p, m, mu_squared),
        function_map: expressions.into_definitions(),
    }
}

/// Explicit alias for [d0].
pub fn d0_mapped(p: [&Atom; 6], m: [&Atom; 4], mu_squared: &Atom) -> MappedLaurentSeries {
    d0(p, m, mu_squared)
}

pub(crate) fn register_boxes(map: &mut FunctionMap) {
    let symbols =
        core::array::from_fn::<_, 11, _>(|i| symbolica::symbol!(format!("__olo_d0_arg_{i}")));
    let x = symbols.map(Symbol::to_atom);
    // One definition for every relabelling of box_finite_four_masses. It binds
    // the sector parameters: new symbols would renumber later ones and reorder
    // the other families' evaluators.
    map.add_function_with_options(
        box_four_masses_chart_symbol(),
        symbols[..10].to_vec(),
        box_four_masses_chart(
            [&x[0], &x[1], &x[2], &x[3], &x[4], &x[5]],
            [&x[6], &x[7], &x[8], &x[9]],
        ),
        native_function_options(),
    )
    .expect("unique four-mass chart");
    let p = [&x[0], &x[1], &x[2], &x[3], &x[4], &x[5]];
    let m = [&x[6], &x[7], &x[8], &x[9]];
    let definitions = [
        (0, box_massless(p, &x[10])),
        (1, box_one_mass(p, m[3], &x[10])),
        (3, box_two_adjacent_ir(p, m[2], m[3], &x[10])),
        (5, box_two_opposite_ir(p, m[1], m[3], &x[10])),
        (7, box_three_masses_ir(p, m, &x[10])),
        (15, finite_series(box_finite_four_masses(p, m))),
    ];
    for (sector, series) in definitions {
        for (coefficient, body) in series.into_coefficients().into_iter().enumerate() {
            map.add_function_with_options(
                symbolica::symbol!(format!("__olo_d0_sector_{sector}_{coefficient}")),
                symbols.to_vec(),
                body,
                native_function_options(),
            )
            .expect("unique box sector");
        }
    }
    for (sector, canonical, shift) in [
        (2, 1, 3),
        (4, 1, 2),
        (6, 3, 3),
        (8, 1, 1),
        (9, 3, 1),
        (10, 5, 3),
        (11, 7, 1),
        (12, 3, 2),
        (13, 7, 2),
        (14, 7, 3),
    ] {
        let p_order = match shift {
            1 => [1, 2, 3, 0, 5, 4],
            2 => [2, 3, 0, 1, 4, 5],
            _ => [3, 0, 1, 2, 5, 4],
        };
        let args = p_order
            .map(|i| x[i].clone())
            .into_iter()
            .chain((0..4).map(|i| x[6 + (i + shift) % 4].clone()))
            .chain([x[10].clone()])
            .collect::<Vec<_>>();
        for coefficient in 0..3 {
            let body = symbolica::symbol!(format!("__olo_d0_sector_{canonical}_{coefficient}"))
                .call(&args);
            map.add_function(
                symbolica::symbol!(format!("__olo_d0_sector_{sector}_{coefficient}")),
                symbols.to_vec(),
                body,
            )
            .unwrap();
        }
    }
}
