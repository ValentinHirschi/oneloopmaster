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
    // The mass pair's arguments `u v14 = X r14` reach one where a root meets
    // `1/r14`, at which the quadratic is `(r24 - r12 r14)(r34 - r13 r14)/r14²`:
    // next to a collinear limit (`a` small) with a mass next to its shell,
    // `X2 r14 - 1 = 7e-8` in a 0C00 box, and `1 - X2 r14` formed from the
    // rounded root lost 1.2e-10 of the value on the wound sheet's
    // `log(1 - x)`, where OneLOop's sheets kept 8e-16. The offsets of the
    // roots from `1/r14` keep it, as the channels' complements do.
    let w14 = root_complements(
        quadratic,
        &r14,
        &(&b - Atom::num(2) * &a / &r14),
        &((&r24 - &r12 * &r14) * (&r34 - &r13 * &r14) / (&r14 * &r14)),
    );
    // Denner and Dittmaier's one-mass box on its own sheets ([one_mass_sheets]).
    let [w13a, w13b] = w13;
    let [w12a, w12b] = w12;
    let [w14a, w14b] = w14;
    let complex_mass = one_mass_sheets_symbol().call(
        [
            x1.clone(),
            x2.clone(),
            &sm / &sr,
            -p[1] / (mass - p[5]),
            w13a,
            w13b,
            -p[4] / (mass - p[3]),
            w12a,
            w12b,
            (mass - p[2]) / mass,
            w14a,
            w14b,
            Q::lower((mass - p[3]) / -p[0]).log() + Q::lower((mass - p[5]) / mass).log(),
        ]
        .as_slice(),
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

/// The terms of [box_finite_one_mass] and [box_finite_two_adjacent] for
/// masses with a nonzero width, on the sheets of Denner and Dittmaier's
/// one- and two-mass boxes (arXiv:1005.2076, eqs. (3.82)-(3.83) and
/// (3.78)-(3.80); their line 0 is the massive line 4 here, their line 3 the
/// line 3), which hold for any real momenta and complex masses.
///
/// Their variables are the roots `u = X rotation` of the unnormalized
/// quadratic (`rotation` is `sqrt(m4²)/|sqrt(m4²)|` with one mass and
/// `sqrt(m4²)/sqrt(m3²)` with two) and the ratios `v` of the dilogarithm
/// pairs, each with its principal logarithm: `log(X r34/r13)` is
/// `log u + log(Y23/Y02)`. OneLOop's boxf1 and boxf2 take
/// `log X + log r34 - log r13` instead: equal for real masses, but the
/// normalization by the complex square roots turns `X` against `u`, across
/// the negative real axis for a root next to it. In `4 1.0 P 14.78 12.26
/// -18.32 15.39 -17.48 0 0 0 0 0 0 3.018 -0.302`
/// (tests/d0_one_mass_width.rs) that happened at P = -0.055: OneLOop (both
/// precisions) and these expressions jumped there from -0.0250 to 0.0254 +
/// 1.87i, and at P = -1e-6 gave 6.11 + 383i and 3.03 + 192i for the -0.0599 -
/// 8.2e-4i of COLLIER and of a direct Feynman-parameter integral. A root
/// exactly on the axis takes its lower lip, the limit of the side that keeps
/// the formula.
///
/// The terms are `(F(u1) - F(u2))/(X1 - X2)` with `F(u) = Li2(u, Y23/Y02) +
/// Li2(u, Y13/Y01) - sum_k Li2(u, r03k) - log(u) log(Y01 Y02/(Y12 m4²))`,
/// where the `r03k` are `Y03/m4²` with one mass and the roots `r031, r032` of
/// `m4² + Y03 x + m3² x²` (`-1/r03k`) with two. When the roots lie on
/// opposite lips of the negative real axis (next to each other near a
/// collinear limit, where the dilogarithm pairs cancel), their principal
/// sheets differ by a winding and the terms grow as `1/(X1 - X2)` and cancel
/// each other (in binary64 a one-mass box next to a collinear limit went from
/// 1e-16 to 6e-8 on these sheets).
/// `u2` then takes `u1`'s lip (`F'`), which keeps each pair on one sheet,
/// and the change `F'(u2) - F(u2)` is added back. It is a multiple of `4 pi²`
/// (the dilogarithms' windings cancel against the logarithm's, since
/// `(1 - u Y23/Y02) (1 - u Y13/Y01)/prod_k (1 - u r03k)` is
/// `Y12 m4²/(Y01 Y02)` at a root), exact once rounded; its rounding residue
/// would otherwise grow as `1/(X1 - X2)` too. The one- and two-mass terms
/// label the roots so that `u2` is the one nearer the axis
/// ([one_mass_sheets], [two_mass_sheets]).
///
/// `channels` are the pairs with the ratios `Y23/Y02` and `Y13/Y01`, their
/// complements `1 - X r34/r13` and `1 - X r24/r12`, and their coefficients
/// `(y1 - y2)/(X1 - X2)`; `masses` the pairs of the `r03k` with their
/// coefficients and their complements `1 - u r03k` at the two roots;
/// `invariants` is `log(Y01/Y12) + log(Y02/m4²)`. With
/// `vanishing_channels` a channel pair whose coefficient vanishes (its ratio
/// is zero: a real mass on shell with a massless line) is dropped.
fn complex_mass_sheets<const N: usize>(
    [x1, x2]: [&Atom; 2],
    rotation: &Atom,
    channels: [(Q, &[Atom; 2], Atom); 2],
    masses: [(Q, Atom, [Atom; 2]); N],
    invariants: Atom,
    vanishing_channels: bool,
) -> Atom {
    let [u1, u2] = [x1, x2].map(|x| x * rotation);
    let upper = |u: &Atom| negative(&-im(u));
    let turns = (upper(&u1) - upper(&u2)) * negative(&re(&u1)) * negative(&re(&u2));
    let first = Q::lower(u1);
    let principal = Q::lower(u2);
    let second = principal.wound(&turns);
    let pair = |v: &Q| [&first * v, &second * v];
    let change =
        |v: &Q, w: &Atom| (&second * v).dilog_complement(w) - (&principal * v).dilog_complement(w);
    let mut one_sheet = (&first / &second).log_over_one_minus() * &invariants / x2;
    let mut windings =
        -Atom::num(2) * sheet_exact::i() * Symbol::PI.to_atom() * &turns * invariants;
    for (v, w, coefficient) in &channels {
        let [y1, y2] = pair(v);
        let term = complement_difference(&y1, &w[0], &y2, &w[1]) * coefficient;
        let winding = change(v, &w[1]);
        if vanishing_channels {
            one_sheet += if_nonzero(coefficient, term);
            windings += if_nonzero(coefficient, winding);
        } else {
            one_sheet += term;
            windings += winding;
        }
    }
    for (v, coefficient, w) in &masses {
        let [y1, y2] = pair(v);
        one_sheet -= complement_difference(&y1, &w[0], &y2, &w[1]) * coefficient;
        windings -= change(v, &w[1]);
    }
    let four_pi_squared = Atom::num(4) * Symbol::PI.to_atom().pow(2);
    let correction = nearest_integer(&(re(&windings) / &four_pi_squared)) * four_pi_squared;
    one_sheet + if_nonzero(&turns, correction / (x1 - x2))
}

/// 1 where [complex_mass_sheets] should take the roots `X1, X2` in the other
/// order, so that `u2 = X2 rotation` is the one nearer the negative real axis
/// (in angle): where `u1, u2` lie on its opposite lips within an eighth of
/// `|X1| + |X2|` of each other and `u1` is the nearer ([two_mass_sheets],
/// [one_mass_sheets]); 0 otherwise.
fn nearer_root_first(x1: &Atom, x2: &Atom, rotation: &Atom) -> Atom {
    let [u1, u2] = [x1, x2].map(|x| x * rotation);
    let upper = |u: &Atom| negative(&-im(u));
    let turns = (upper(&u1) - upper(&u2)) * negative(&re(&u1)) * negative(&re(&u2));
    let abs = |z: &Atom| Symbol::ABS.call((z,));
    let close = negative(&(Atom::num(8) * abs(&(x1 - x2)) - abs(x1) - abs(x2)));
    let nearer_first = negative(&(abs(&im(&u1)) * abs(&u2) - abs(&im(&u2)) * abs(&u1)));
    if_nonzero(&turns, close * nearer_first)
}

/// `[a, b]`, or `[b, a]` where `swap` is nonzero.
fn swapped(swap: &Atom, a: &Atom, b: &Atom) -> [Atom; 2] {
    [
        if_nonzero_else(swap, b.clone(), a.clone()),
        if_nonzero_else(swap, a.clone(), b.clone()),
    ]
}

fn one_mass_sheets_symbol() -> Symbol {
    symbolica::symbol!("__olo_box_one_mass_sheets")
}

fn one_mass_sheets_ordered_symbol() -> Symbol {
    symbolica::symbol!("__olo_box_one_mass_sheets_ordered")
}

/// [one_mass_sheets_ordered] with the roots labelled so that `u2` is the one
/// nearer the negative real axis ([nearer_root_first]), as [two_mass_sheets]
/// does for two masses; the terms are symmetric in the labels, with each
/// root's complements. Next to a collinear limit with `r13 r24` next to
/// `r12 r34` the roots of the one-mass quadratic are next to each other
/// (`X1 - X2` from 1e-4 to 1e-2 of `X`) and, with a width, on opposite lips:
/// taken across from the farther root, the channel pairs carried windings
/// 2e2 to 7e3 times the terms, which cancelled between the two channels. Of
/// the relabellings of 9f663c5's 1000 boxes next to a collinear limit and to
/// `r13 r24 = r12 r34`, 16 lost 1.7e-12 to 6.4e-12 in binary64 that way,
/// where OneLOop's sheets (09b2c74) kept 8e-16 and a one-ulp change of an
/// input moves the value by 1.2e-16 (tests/d0_one_mass_width.rs); from the
/// nearer root they keep 4.5e-16. Traced through the binary64 expressions,
/// 162 of 75,996 labellings of complex one-mass boxes (the verifier's and the
/// streams' sweep rows, 9f663c5's sweeps) take the other order: 21 become
/// more than ten times more accurate (1.3e-13 to 5e-12 before, at most 1e-14
/// after), none less. A registered function of its own, so that the swapped
/// arguments are not repeated in the terms (the one-mass sector shrinks from
/// 4.1 MB to 1.6 MB); it binds the sector parameters and two of
/// box_complex's, so that no parameter symbol is created.
fn one_mass_sheets(
    [
        x1,
        x2,
        rotation,
        v13,
        w13a,
        w13b,
        v12,
        w12a,
        w12b,
        v14,
        w14a,
        w14b,
        invariants,
    ]: [&Atom; 13],
) -> Atom {
    let swap = nearer_root_first(x1, x2, rotation);
    let pick = |a: &Atom, b: &Atom| swapped(&swap, a, b);
    let [x1, x2] = pick(x1, x2);
    let [w13a, w13b] = pick(w13a, w13b);
    let [w12a, w12b] = pick(w12a, w12b);
    let [w14a, w14b] = pick(w14a, w14b);
    one_mass_sheets_ordered_symbol().call(
        [
            x1,
            x2,
            rotation.clone(),
            v13.clone(),
            w13a,
            w13b,
            v12.clone(),
            w12a,
            w12b,
            v14.clone(),
            w14a,
            w14b,
            invariants.clone(),
        ]
        .as_slice(),
    )
}

/// The terms of [box_finite_one_mass] for a mass with a width
/// ([complex_mass_sheets]), as a registered function of the roots `X1, X2`,
/// the `rotation` (`u = X rotation`), the channel ratios `Y23/Y02` and
/// `Y13/Y01` with the complements `1 - X r34/r13` and `1 - X r24/r12`, the
/// mass pair's ratio `Y03/m4²` with the complements `1 - X r14`, and the
/// invariants' logarithm. The coefficients `r34/r13`, `r24/r12` and `r14` are
/// the ratios times `rotation`.
fn one_mass_sheets_ordered(
    [
        x1,
        x2,
        rotation,
        v13,
        w13a,
        w13b,
        v12,
        w12a,
        w12b,
        v14,
        w14a,
        w14b,
        invariants,
    ]: [&Atom; 13],
) -> Atom {
    let w13 = [w13a.clone(), w13b.clone()];
    let w12 = [w12a.clone(), w12b.clone()];
    complex_mass_sheets(
        [x1, x2],
        rotation,
        [
            (Q::lower(v13.clone()), &w13, rotation * v13),
            (Q::lower(v12.clone()), &w12, rotation * v12),
        ],
        [(Q::lower(v14.clone()), rotation * v14, [w14a.clone(), w14b.clone()])],
        invariants.clone(),
        false,
    )
}

fn two_mass_sheets_symbol() -> Symbol {
    symbolica::symbol!("__olo_box_two_mass_sheets")
}

fn two_mass_sheets_ordered_symbol() -> Symbol {
    symbolica::symbol!("__olo_box_two_mass_sheets_ordered")
}

/// [two_mass_sheets_ordered] with the roots labelled so that `u2 = X2
/// rotation` is the one nearer the negative real axis (in angle) when they
/// lie on its opposite lips next to each other; the terms are symmetric in
/// the labels, with each root's complements. [complex_mass_sheets] takes `u2`
/// across the axis to `u1`'s lip and adds back the change of `F` there. Taken
/// from the farther root, the change spread over the pairs as
/// `2 pi i log(1 - y)`-sized windings, up to 1e6 times the result next to a
/// collinear limit (channel arguments `y` next to one), which cancelled
/// against the correction `4 pi² k/(X1 - X2)`: with `u2` taken across
/// whatever its angle (a5da789), two-mass boxes there lost 1e-12 to 1.6e-9,
/// where OneLOop's sheets kept 1e-16 (TM 66556 and 63221 of
/// tests/data/two_mass_conditioning.txt). Taken from the nearer root it keeps
/// 4e-16; of 376 relabellings of the two-mass sweeps where the choice changed
/// the binary64 error tenfold, the nearer root was the better one in 368 (the
/// other 8 by at most 1.6e-12), and in all of them the roots were within 0.08
/// of their size. Farther apart than an eighth of `|X1| + |X2|` the
/// correction is not much larger than the value, and the labels are kept:
/// relabelled there, 128-bit `Float` lost up to 13 more bits (numerica rounds
/// the near-shell differences in the channels' complements to their tracked
/// precision) without a binary64 gain. A registered function of its own, so
/// that the swapped arguments are not repeated in the terms.
fn two_mass_sheets(
    [
        x1,
        x2,
        rotation,
        v13,
        w13a,
        w13b,
        v12,
        w12a,
        w12b,
        r14,
        invariants,
        product,
    ]: [&Atom; 12],
) -> Atom {
    let swap = nearer_root_first(x1, x2, rotation);
    let pick = |a: &Atom, b: &Atom| swapped(&swap, a, b);
    let [x1, x2] = pick(x1, x2);
    let [w13a, w13b] = pick(w13a, w13b);
    let [w12a, w12b] = pick(w12a, w12b);
    two_mass_sheets_ordered_symbol().call(
        [
            x1,
            x2,
            rotation.clone(),
            v13.clone(),
            w13a,
            w13b,
            v12.clone(),
            w12a,
            w12b,
            r14.clone(),
            invariants.clone(),
            product.clone(),
        ]
        .as_slice(),
    )
}

/// The terms of [box_finite_two_adjacent] for masses with a width
/// ([complex_mass_sheets]), as a registered function of the roots `X1, X2`,
/// the `rotation` (`u = X rotation`), the channel ratios `Y23/Y02` and
/// `Y13/Y01` with the complements `1 - X r34/r13` and `1 - X r24/r12`, `r14`
/// (the pairs `r031 = r14/rotation` and `r032 = 1/(r14 rotation)`, with
/// coefficients `r14` and `1/r14`), the invariants' logarithm and their
/// `product = r12 r13/r23` (`Y01 Y02/(Y12 m4²)`, the logarithm's argument).
/// Inlined, every occurrence of a root or ratio carried its whole expression:
/// the terms took 22 MB per copy of the finite box (four in sector 5's body),
/// and the definitions grew from 97 MB to 209 MB
/// (tests/expression_construction.rs). Like the four-mass chart it binds the
/// sector parameters (and one parameter of box_complex's definitions), so
/// that no new parameter symbols renumber later ones.
///
/// The mass pairs' arguments `u r031 = X r14` and `u r032 = X/r14` reach one
/// where a root of the box's quadratic meets a root `-1/r03k` of
/// `m4² + Y03 x + m3² x²`: in Denner and Dittmaier's variables the quadratic
/// is `(Y01 - Y13/r03k)(Y02 - Y23/r03k)` there, small when a real mass is
/// next to its shell (`Y01` or `Y02` small) and its partner's `Y13` or `Y23`
/// is small too (a complex mass whose real part is next to its shell, with a
/// narrow width, or another invariant next to the real mass). Formed from the
/// rounded argument, `1 - X r14` kept only `|1 - X r14|` relative digits,
/// and the dilogarithms on wound sheets (which carry `log(1 - x)`) and their
/// divided differences lost the rest: binary64 lost 3.3e-11 in
/// tests/data/two_mass_shell.txt's first row (a real mass on shell to 2e-16
/// and a complex one whose real part is), 4.8e-9 with a real mass next to its
/// shell in two invariants (1e-6 and 4e-8), and up to 0.16 of the value with
/// both masses on or next to their shells. At a root the product of the two
/// complements, `X² - k14 X + 1`, is `(1 - X r34/r13)(1 - X r24/r12)
/// r12 r13/r23` (the quadratic is `(r34 X - r13)(r24 X - r12) - r23 (X² -
/// k14 X + 1)`): the channels' complements, which the caller forms from the
/// kinematics, times `product`. Both are small next to `r14 = 1`, the
/// threshold `Y03² = 4 m3² m4²` (equal masses and a vanishing invariant
/// between them), where `r14` itself is known only to the square root of its
/// rounding: the direct forms share that error and it cancels between the two
/// pairs, which the product, formed without `r14`, would not (two such sweep
/// boxes lost 1e-9 that way, where the direct forms keep 5e-16 and 1.6e-12).
/// With `d1 = 1 - X r14` small, `X` is next to `1/r14` and the other
/// complement `d2` next to `1 - 1/r14²`: `r14`'s relative error is about
/// `2 eps/|d2|`, the product over `d2` keeps `d1` to about `2 eps/|d2|²`, and
/// the direct difference (whose `r14` error cancels) to `eps/|d1|`. A
/// complement below 1/16 is that product over the other where the other is
/// above 1/4 (the direct difference would keep fewer digits below 1/16 than
/// the product's few roundings), and also where the complement is below half
/// the other's square and the direct difference loses binary64 digits: the
/// complement below 2^-13 (`2^-40 d` negligible against one at the working
/// precision, which `DoubleFloat` and `Float` never reach there). Otherwise
/// both are taken directly. With the other above 1/4 alone, complements
/// `1 - X/r14` of 9e-7 to 3e-5 (the other 0.01 to 0.19) were formed from the
/// rounded root, and boxes with both masses next to their shells lost 5e-12
/// to 2.8e-10 in binary64 where OneLOop's sheets kept 6e-14 to 7e-12
/// (tests/data/two_mass_conditioning.txt). At 128 bits the direct
/// difference keeps them, while numerica rounds the product's near-shell
/// differences to their tracked precision: without the bound the product lost
/// up to 30 bits there (a real mass within 1e-15 of its shell, 123 to 93
/// tracked bits). A channel with a vanishing ratio (a real mass on shell)
/// contributes the complement one; its given complement is not finite then.
fn two_mass_sheets_ordered(
    [
        x1,
        x2,
        rotation,
        v13,
        w13a,
        w13b,
        v12,
        w12a,
        w12b,
        r14,
        invariants,
        product,
    ]: [&Atom; 12],
) -> Atom {
    let w13 = [w13a.clone(), w13b.clone()];
    let w12 = [w12a.clone(), w12b.clone()];
    let abs = |z: &Atom| Symbol::ABS.call((z,));
    // 1 when `x >= 0` is at most half an ulp of 1 at the working precision.
    let negligible =
        |x: Atom| if_nonzero_else(&(grouped(&(1 + x)) - 1), Atom::num(0), Atom::num(1));
    // 1 when `d` is below 1/16 and `other` above 1/4, or `d` below half the
    // square of `other` and below 2^-13 in binary64 (2^-40 d negligible).
    let small = |d: &Atom, other: &Atom| {
        negative(&(Atom::num(16) * abs(d) - 1))
            * if_nonzero_else(
                &negative(&(1 - Atom::num(4) * abs(other))),
                Atom::num(1),
                negative(&(Atom::num(2) * abs(d) - abs(other).pow(2)))
                    * negligible(abs(d) / Atom::num(2).pow(40)),
            )
    };
    let [[c1a, c2a], [c1b, c2b]] =
        [(x1, &w13[0], &w12[0]), (x2, &w13[1], &w12[1])].map(|(x, w13, w12)| {
            let d1 = 1 - x * r14;
            let d2 = 1 - x / r14;
            let both = if_nonzero_else(v13, w13.clone(), Atom::num(1))
                * if_nonzero_else(v12, w12.clone(), Atom::num(1))
                * product;
            [
                if_nonzero_else(&small(&d1, &d2), &both / &d2, d1.clone()),
                if_nonzero_else(&small(&d2, &d1), &both / &d1, d2.clone()),
            ]
        });
    complex_mass_sheets(
        [x1, x2],
        rotation,
        [
            (Q::lower(v13.clone()), &w13, rotation * v13),
            (Q::lower(v12.clone()), &w12, rotation * v12),
        ],
        [
            (Q::lower(r14 / rotation), r14.clone(), [c1a, c1b]),
            (
                Q::lower(Atom::num(1) / (r14 * rotation)),
                Atom::num(1) / r14,
                [c2a, c2b],
            ),
        ],
        invariants.clone(),
        true,
    )
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
    let w13 = complements(&r13 / &r34, -&gap);
    let term_13 = if_nonzero_else(
        &r13,
        if_nonzero(
            &r34,
            complement_difference(&x1_13, &w13[0], &x2_13, &w13[1]) * &r34 / &r13,
        ),
        limit_13,
    );
    let w12 = complements(&r12 / &r24, gap.clone());
    let term_12 = if_nonzero_else(
        &r12,
        if_nonzero(
            &r24,
            complement_difference(&x1_12, &w12[0], &x2_12, &w12[1]) * &r24 / &r12,
        ),
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
    // Denner and Dittmaier's two-mass box (their eq. (3.78), lines 0 and 3
    // the massive lines 4 and 3 here) on its own sheets ([two_mass_sheets]).
    // With a real m4², Y01/Y12 and Y02/m4² are real: their logarithms are the
    // differences of those of the -i0 invariants (exact for complex ones).
    let rotation = &sm4 / &sm3;
    let [y01, y02, y13, y23] = [mass_4 - p[3], mass_4 - p[5], mass_3 - p[4], mass_3 - p[1]];
    let log = |y: &Atom| Q::lower(y.clone()).log();
    let sheets = |arguments: [Atom; 12]| two_mass_sheets_symbol().call(arguments.as_slice());
    let [w13a, w13b] = w13.clone();
    let [w12a, w12b] = w12.clone();
    let direct = sheets([
        x1.clone(),
        x2.clone(),
        rotation.clone(),
        &y23 / &y02,
        w13a,
        w13b,
        &y13 / &y01,
        w12a,
        w12b,
        r14.clone(),
        log(&y01) - log(&-p[0]) + log(&y02) - log(mass_4),
        &r12 * &r13 / &r23,
    ]);
    // A real m4² on shell with a massless line (Y01 or Y02 zero) takes the
    // formula with lines 0 and 3 exchanged, as Denner and Dittmaier state: its
    // quadratic has a and c exchanged, so its roots are 1/u (1/X here), its
    // rotation, ratios Y02/Y23 and Y01/Y13 and r03k the inverses, and the
    // terms above are its terms divided by X1 X2. A real m3² on shell (Y13 or
    // Y23 zero) drops that pair, there or in the exchanged formula. The
    // complements of the inverse arguments 1/y are -w/y, with y = X r34/r13 or
    // X r24/r12 itself: as -w/(1 - w) they lost y's digits for a small y.
    let [w13a, w13b] = [(&x1, &w13[0]), (&x2, &w13[1])].map(|(x, w)| -w / (x * (&r34 / &r13)));
    let [w12a, w12b] = [(&x1, &w12[0]), (&x2, &w12[1])].map(|(x, w)| -w / (x * (&r24 / &r12)));
    let exchanged = sheets([
        Atom::num(1) / &x1,
        Atom::num(1) / &x2,
        &sm3 / &sm4,
        &y02 / &y23,
        w13a,
        w13b,
        &y01 / &y13,
        w12a,
        w12b,
        r14.clone(),
        log(&y13) - log(&-p[0]) + log(&y23) - log(mass_3),
        &r24 * &r34 / &r23,
    ]) / (&x1 * &x2);
    let complex_masses = if_nonzero_else(
        &if_nonzero(&r12, if_nonzero(&r13, Atom::num(1))),
        direct,
        exchanged,
    );
    let width = any_nonzero([im(mass_3), im(mass_4)]);
    let scale = &a * &sr * &sr * sm3 * sm4;
    // Each branch is scaled on its own: real masses keep their arithmetic.
    let finite = if_nonzero_else(&width, -complex_masses / &scale, -finite / scale);
    if_nonzero(&a, finite)
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
    let general = with_box_normalization(LaurentSeries::new(
        -Atom::num(2) * z0 / &r13,
        -z1 / r13,
        Atom::num(0),
    ));
    // The chart exchange can put p4 exactly on the massive shell even when
    // the original labelling selected box10. Reclassify that infrared limit:
    // box10 divides by m-p4 and has no double pole, whereas the existing
    // box09 expression supplies the finite limit and both poles.
    if_series(
        &(&p4 - mass),
        general,
        box_one_mass_09(&p2, p3, &p12, &p23, mass, mu_squared),
    )
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
    let general = with_box_normalization(LaurentSeries::new(
        (li2f + li2b + Atom::num(2) * li2e
            - Atom::num(2) * li2d
            - Atom::num(2) * &logd * q13.log())
            / &denominator,
        logd / denominator,
        Atom::num(0),
    ));
    // The chart exchange can reach one or both exact mass shells even when
    // the original chart selected box13. Reclassify before its zero ratios
    // are evaluated, using the already-normalized infrared limits.
    if_series(
        &(&p2 - mass_3),
        if_series(
            &(&p4 - mass_4),
            general,
            box_two_adjacent_12(p3, &p2, &p23, &p12, mass_4, mass_3, mu_squared),
        ),
        if_series(
            &(&p4 - mass_4),
            box_two_adjacent_12(p3, &p4, &p12, &p23, mass_3, mass_4, mu_squared),
            box_two_adjacent_11(p3, &p12, &p23, mass_3, mass_4, mu_squared),
        ),
    )
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
    // Likewise one definition of the two-mass terms with widths. Its twelfth
    // parameter is box_complex's (registered before), so that no symbol is
    // created here either.
    let sheet_symbols: [Symbol; 12] = core::array::from_fn(|i| {
        if i < 11 {
            symbols[i]
        } else {
            symbolica::symbol!("__olo_boxc_arg_11")
        }
    });
    map.add_function_with_options(
        two_mass_sheets_ordered_symbol(),
        sheet_symbols.to_vec(),
        two_mass_sheets_ordered(sheet_symbols.map(Symbol::to_atom).each_ref()),
        native_function_options(),
    )
    .expect("unique two-mass sheets with ordered roots");
    map.add_function_with_options(
        two_mass_sheets_symbol(),
        sheet_symbols.to_vec(),
        two_mass_sheets(sheet_symbols.map(Symbol::to_atom).each_ref()),
        native_function_options(),
    )
    .expect("unique two-mass sheets");
    // And of the one-mass terms with a width, with two of box_complex's
    // parameters.
    let one_mass_symbols: [Symbol; 13] = core::array::from_fn(|i| {
        if i < 11 {
            symbols[i]
        } else {
            symbolica::symbol!(format!("__olo_boxc_arg_{i}"))
        }
    });
    map.add_function_with_options(
        one_mass_sheets_ordered_symbol(),
        one_mass_symbols.to_vec(),
        one_mass_sheets_ordered(one_mass_symbols.map(Symbol::to_atom).each_ref()),
        native_function_options(),
    )
    .expect("unique one-mass sheets with ordered roots");
    map.add_function_with_options(
        one_mass_sheets_symbol(),
        one_mass_symbols.to_vec(),
        one_mass_sheets(one_mass_symbols.map(Symbol::to_atom).each_ref()),
        native_function_options(),
    )
    .expect("unique one-mass sheets");
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
