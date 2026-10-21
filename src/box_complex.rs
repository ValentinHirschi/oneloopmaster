//! Native complex-mass box representation from OneLOop's `avh_olo_boxc`.
//!
//! The logarithms, dilogarithms, contour residues and permutation choices are
//! transparent Symbolica definitions. Infinitesimal root lips are carried by
//! `SheetAtom`; no finite regulator or numerical integral callback is used.
#![allow(clippy::too_many_arguments)] // Match the source's named polynomial coefficients.
use super::*;
use sheet_exact::{
    SheetAtom as Q, complement_difference, divided_difference as dd, i, im, negative, re,
    sign_nonnegative,
};

fn function(name: &str) -> Symbol {
    symbolica::symbol!(format!("__olo_boxc_{name}"))
}

fn call<const N: usize>(name: &str, args: [&Atom; N]) -> Atom {
    function(name).call(args.map(Clone::clone).as_slice())
}

/// The squared norm `re(z)² + im(z)²`, as a call of the registered `norm`.
///
/// Crossed products in `z conj(z)` need not cancel when the two components
/// carry different tracked precisions, and a product `s z conj(z)` that
/// Symbolica flattens (a real factor `s` applied to `z` first) keeps a
/// rounding residue in its imaginary part even in binary64. The body
/// ([register_norm]) groups the operand and its real components so that
/// conjugation and expansion cannot restore those products, and is exactly
/// real on every route. It is a call so that each use holds one copy of its
/// operand: written out, the grouping holds 16 copies, and nested norms
/// multiply them.
fn norm(z: &Atom) -> Atom {
    call("norm", [z])
}

fn norm_body(z: &Atom) -> Atom {
    let z = grouped(z);
    let real = grouped(&re(&z));
    let imaginary = grouped(&im(&z));
    &real * &real + &imaginary * &imaginary
}

fn less(a: &Atom, b: &Atom) -> Atom {
    negative(&(a - b))
}

fn le(a: &Atom, b: &Atom) -> Atom {
    1 - less(b, a)
}

fn min(a: &Atom, b: &Atom) -> Atom {
    if_nonzero_else(&less(a, b), a.clone(), b.clone())
}

/// Nonzero where `2^bits x < y` for the real `x, y >= 0` (squared norms):
/// the chart tests compare at binary64's half precision with fixed powers of
/// two, so that the hand-written, generated, interpreted, SymJIT and
/// `DoubleFloat` routes take the same chart on the same binary64 inputs
/// ([tfun_chart], [zero_external]).
fn below(bits: u32, x: &Atom, y: &Atom) -> Atom {
    less(&(Atom::num(1u64 << bits) * x), y)
}

/// `|v - w| < 2^-26 |w|`: the two values agree to binary64's half precision.
fn agree(v: &Atom, w: &Atom) -> Atom {
    below(52, &norm(&(v - w)), &norm(w))
}

/// `|v - w| < 2^-bits |w|`: a chart moved within the rounding of its line or
/// of its K keeps `bits` bits ([tfun_chart]).
fn close(bits: u32, v: &Atom, w: &Atom) -> Atom {
    below(2 * bits, &norm(&(v - w)), &norm(w))
}

/// Zero where `z` is infinite or NaN (`z - z` is not zero there).
fn finite(z: &Atom) -> Atom {
    if_nonzero_else(&(grouped(z) - z), Atom::num(0), Atom::num(1))
}

/// NaN, formed at run time from `z` (a quotient of two zeros that Symbolica
/// cannot cancel): the value of a chart that the working precision cannot
/// resolve, which the caller detects (and NinjaRS rescues in higher
/// precision) instead of a silent wrong value.
fn undetermined(z: &Atom) -> Atom {
    let twice = z * Atom::num(2);
    (grouped(z) - z) / (grouped(&twice) - &twice)
}

/// Nonzero where the arithmetic rounds like binary64: `2^60 t + t` rounds to
/// `2^60 t` for a positive real `t`. `DoubleFloat` and `Float` keep the sum
/// (and give zero), so that their charts are not tested at binary64's
/// rounding (they resolve what binary64 does not; NinjaRS's rescue relies on
/// them).
fn rounds_like_binary64(t: &Atom) -> Atom {
    let scaled = grouped(&(t * Atom::num(1i64 << 60)));
    let sum = grouped(&(&scaled + t));
    if_nonzero_else(&(sum - scaled), Atom::num(0), Atom::num(1))
}

/// Nonzero where the arithmetic is `DoubleFloat`: it keeps `t + 2^-1000 t`
/// exactly for the exact `t = 1 + first` (1 or 2, with `first` 0 or 1), as
/// its sums of two doubles never round, while binary64 and `Float` of fewer
/// than 1000 bits round it to `t`. (A test that follows the working precision,
/// `1 + x` rounding to 1, is therefore never true in `DoubleFloat`.)
fn double_double(first: &Atom) -> Atom {
    let t = grouped(&(1 + first));
    let power = Atom::num(symbolica::domains::integer::Integer::from(2).pow(1000));
    // The addend is formed at run time: a constant `1 + 2^-1000` would reach
    // `DoubleFloat` rounded to 1 by its conversion.
    let small = grouped(&(&t / &power));
    let moved = grouped(&(&t + &small));
    if_nonzero_else(&(moved - &t), Atom::num(1), Atom::num(0))
}

/// The sum of `terms`, added left to right as in OneLOop's `boxc`.
///
/// Symbolica orders the terms of a sum canonically, so `m2 - m3 - p2` need
/// not subtract the masses first. For equal masses and a tiny `p2`, forming
/// `(m2 - p2) - m3` instead keeps the value `-p2` in `Float`, but tracks it with
/// only the bits left after the cancellation, and every later operation on it
/// runs at that precision (a D0 with `p2 = 2^-114` lost 129 bits from 132
/// bits up). Each partial sum is kept by [grouped].
fn in_order(terms: &[Atom]) -> Atom {
    terms[1..]
        .iter()
        .fold(terms[0].clone(), |sum, term| grouped(&sum) + term)
}

/// The coefficient sum `terms[0] + terms[1] + ...` of a quadratic, or `exact`,
/// the invariant it equals, when the sum cancels by more than 2^16.
///
/// A tiny invariant (`p2 = 2^-114` in a D0 with masses 29929, computed as
/// `(2^-114 + 40000) - 40000`) is lost to cancellation, in `Float` with the
/// tracked precision of the difference. Otherwise the sum is formed from the
/// coefficients themselves, as in OneLOop: the edge quadratics then stay
/// consistent with the conic they restrict, which kinematics within an ulp of
/// a threshold (an exact double root of one and not the other) rely on.
fn coefficient_sum(terms: &[&Atom], exact: &Atom) -> Atom {
    let numeric = terms.iter().fold(Atom::num(0), |sum, term| sum + *term);
    let scale = terms.iter().fold(Atom::num(0), |sum, term| sum + norm(term));
    let cancels = less(&(Atom::num(1u64 << 32) * norm(exact)), &scale);
    if_nonzero_else(&cancels, exact.clone(), numeric)
}

/// The Källén function `λ(p, mi, mj) = (p - (√mi + √mj)²)(p - (√mi - √mj)²)`,
/// the discriminant of the edge quadratic `p x² + (mi - mj - p) x + mj` (and
/// of `p x² + (mj - mi - p) x + mi`, the same edge from its other end).
///
/// Each factor is formed from the invariants: at a threshold or
/// pseudo-threshold it vanishes exactly, and next to one it keeps the
/// distance to it. The discriminant of the rounded edge coefficients instead
/// leaves a rounding residue there, whose square root moves the roots by
/// sqrt(eps). With a massless line (`mj = 0`) both thresholds are `mi` itself.
fn kallen([p, mi, mj]: [&Atom; 3]) -> Atom {
    let plus = grouped(&in_order(&[
        mi.clone(),
        mj.clone(),
        Atom::num(2) * principal_sqrt(&(mi * mj)),
    ]));
    let difference = grouped(&(mi - mj));
    // (√mi - √mj)² = (mi - mj)²/(√mi + √mj)², without cancellation.
    let minus = if_nonzero_else(
        mj,
        if_nonzero_else(mi, &difference * &difference / &plus, mj.clone()),
        mi.clone(),
    );
    grouped(&(p - &plus)) * grouped(&(p - minus))
}

fn kallen_call(p: &Atom, mi: &Atom, mj: &Atom) -> Atom {
    call("kallen", [p, mi, mj])
}

/// `[λ_a, λ_abc, λ_bc]`: the discriminants of the edge quadratics of the
/// conic `a x² + b y² + c x y + d x + e y + f` on the triangle `tfun`
/// integrates over: `(a, d, f)` on `y = 0`, `(a + b + c, d + e, f)` on `y = x`
/// and `(b, c + e, a + d + f)` on `x = 1`. The conic takes the values `f`,
/// `q10` and `q11` (masses) at the corners (0, 0), (1, 0) and (1, 1), and
/// `abc` is `a + b + c`; the caller forms all of them from the invariants
/// they stand for.
fn edge_kallens(a: &Atom, b: &Atom, abc: &Atom, f: &Atom, q10: &Atom, q11: &Atom) -> [Atom; 3] {
    [
        kallen_call(a, q10, f),
        kallen_call(abc, q11, f),
        kallen_call(b, q11, q10),
    ]
}

fn imaginary_sign(z: &Atom, fallback: &Atom) -> Atom {
    let imaginary = im(z);
    if_nonzero_else(&imaginary, sign_nonnegative(&imaginary), fallback.clone())
}

/// A value and the first two coefficients of its infinitesimal continuation.
/// The coefficients are symbolic directions, never finite numerical offsets.
#[derive(Clone)]
struct Continued {
    value: Atom,
    first: Atom,
    second: Atom,
}

impl Continued {
    fn from_args(x: [&Atom; 3]) -> Self {
        Self {
            value: x[0].clone(),
            first: x[1].clone(),
            second: x[2].clone(),
        }
    }

    fn sign(&self) -> Atom {
        imaginary_sign(
            &self.value,
            &imaginary_sign(&self.first, &imaginary_sign(&self.second, &Atom::num(1))),
        )
    }

    fn sheet(&self) -> Q {
        Q::with_sign(self.value.clone(), self.sign())
    }

    fn complement(&self, constant: i64) -> Self {
        Self {
            value: Atom::num(constant) - &self.value,
            first: -&self.first,
            second: -&self.second,
        }
    }

    fn subtract(&self, other: &Self) -> Self {
        Self {
            value: &self.value - &other.value,
            first: &self.first - &other.first,
            second: &self.second - &other.second,
        }
    }
}

fn mass_directions(m: [&Atom; 4]) -> [Atom; 4] {
    // Original boxc supplies -i*EPSN*abs(Re m) to exactly real masses.
    m.map(|mass| {
        if_nonzero_else(
            &im(mass),
            Atom::num(0),
            -i() * Symbol::ABS.call((re(mass),)),
        )
    })
}

fn continued_root(y: Atom, k: &Atom, l: &Atom, dl: &Atom, dn: &Atom, dn2: &Atom) -> Continued {
    let derivative = Atom::num(2) * k * &y + l;
    continued_root_with_slope(y, derivative, k, dl, dn, dn2)
}

/// [continued_root] given the slope `2 K y + L` of the polynomial at its root.
fn continued_root_with_slope(
    y: Atom,
    derivative: Atom,
    k: &Atom,
    dl: &Atom,
    dn: &Atom,
    dn2: &Atom,
) -> Continued {
    // Insert y + eps*y' + eps^2*y'' into K*y^2+L*y+N=0. The source
    // has mass widths at order EPSN and its explicit -IEPS*abs(Re f)
    // at order EPSN^2. Retaining both coefficients fixes cancellations
    // between equal mass widths without ever choosing a finite epsilon.
    let first = -(dl * &y + dn) / &derivative;
    let second = -(k * &first * &first + dl * &first + dn2) / derivative;
    Continued {
        value: y,
        first,
        second,
    }
}

fn r0_call(y1: &Continued, y2: &Continued) -> Atom {
    call(
        "r0",
        [
            &y1.value, &y1.first, &y1.second, &y2.value, &y2.first, &y2.second,
        ],
    )
}

/// [r0_call] given `delta = y1 - y2` ([r0_line]).
fn r0_line_call(y1: &Continued, y2: &Continued, delta: &Atom) -> Atom {
    call(
        "r0_line",
        [
            &y1.value, &y1.first, &y1.second, &y2.value, &y2.first, &y2.second, delta,
        ],
    )
}

/// [r1_call] of a line edge, given the product `(y1 - z)(y2 - z)` ([r1_line]).
fn r1_line_call(z: &Continued, y1: &Continued, y2: &Continued, fy: &Atom, product: &Atom) -> Atom {
    call(
        "r1_line",
        [
            &z.value, &z.first, &z.second, &y1.value, &y1.first, &y1.second, &y2.value, &y2.first,
            &y2.second, fy, product,
        ],
    )
}

fn r1_call(z: &Continued, y1: &Continued, y2: &Continued, fy: &Atom) -> Atom {
    call(
        "r1",
        [
            &z.value, &z.first, &z.second, &y1.value, &y1.first, &y1.second, &y2.value, &y2.first,
            &y2.second, fy,
        ],
    )
}

/// The two-point subgraph whose Feynman polynomial `p x² + (far - near - p) x
/// + near` an edge quadratic `a x² + b x + c` of `s3` is, from the invariants:
/// `lambda = λ(p, near, far)` (see [kallen]), `p` (`a`), and the masses
/// `near` (`c`, at x = 0) and `far` (`a + b + c`, at x = 1).
struct Edge<'a> {
    lambda: &'a Atom,
    p: &'a Atom,
    near: &'a Atom,
    far: &'a Atom,
}

fn s3_call(
    y1: &Continued,
    y2: &Continued,
    d: &Atom,
    e: &Atom,
    a: &Atom,
    b: &Atom,
    c: &Atom,
    db: &Atom,
    dc: &Atom,
    edge: Edge,
) -> Atom {
    call(
        "s3",
        [
            &y1.value, &y1.first, &y1.second, &y2.value, &y2.first, &y2.second, d, e, a, b, c, db,
            dc, edge.lambda, edge.p, edge.near, edge.far,
        ],
    )
}

/// [s3_call] for a linear `b x + c`, with `far = b + c`.
fn s3_linear_call(
    y1: &Continued,
    y2: &Continued,
    d: &Atom,
    e: &Atom,
    b: &Atom,
    c: &Atom,
    dc: &Atom,
    far: &Atom,
) -> Atom {
    let zero = Atom::num(0);
    let edge = Edge {
        lambda: &zero,
        p: &zero,
        near: c,
        far,
    };
    s3_call(y1, y2, d, e, &zero, b, c, &zero, dc, edge)
}

/// [s3_linear_call] for a line edge of [tfun_selected], given the product
/// `(y1 - z)(y2 - z)` of the mapped poles' offsets from its root `z`
/// ([line_product]).
fn s3_line_call(
    y1: &Continued,
    y2: &Continued,
    d: &Atom,
    e: &Atom,
    b: &Atom,
    c: &Atom,
    dc: &Atom,
    far: &Atom,
    product: &Atom,
) -> Atom {
    call(
        "s3_line",
        [
            &y1.value, &y1.first, &y1.second, &y2.value, &y2.first, &y2.second, d, e, b, c, dc,
            far, product,
        ],
    )
}

fn eta(a: &Atom, sa: &Atom, b: &Atom, sb: &Atom, c: &Atom, sc: &Atom) -> Atom {
    let sa = imaginary_sign(a, sa);
    let sb = imaginary_sign(b, sb);
    let sc = imaginary_sign(c, sc);
    i() * Symbol::PI.to_atom() * &sc * (1 + &sa * &sb) * (1 - sa * sc) / 2
}

fn roots(a: &Atom, b: &Atom, c: &Atom) -> [Atom; 2] {
    [call("root_1", [a, b, c]), call("root_2", [a, b, c])]
}

fn principal_sqrt(value: &Atom) -> Atom {
    // The generic complex sqrt uses polar trigonometry. On the negative axis
    // that can leave a tiny real part, which changes subsequent contour tests.
    let off_axis = off_real_axis(value);
    let value = grouped(value);
    let absolute = Symbol::ABS.call((&value,));
    let axis = if_nonzero_else(&(&value - &absolute), i() * (-&value).sqrt(), value.sqrt());
    // Recover the small rectangular component by the product 2*u*v=Im(z),
    // rather than a cosine near pi/2 or a difference of nearly equal norms.
    let real = re(&value);
    let imaginary = im(&value);
    let u = ((&absolute + &real) / 2).sqrt();
    let v = sign_nonnegative(&imaginary) * ((&absolute - &real) / 2).sqrt();
    let positive = &u + i() * &imaginary / (Atom::num(2) * &u);
    let negative_real = &imaginary / (Atom::num(2) * &v) + i() * &v;
    if_nonzero_else(
        &off_axis,
        if_nonzero_else(&negative(&real), negative_real, positive),
        axis,
    )
}

fn root_bodies(a: &Atom, b: &Atom, c: &Atom) -> [Atom; 2] {
    // Same root labels and cancellation-avoiding product relation as solabc.
    let delta = principal_sqrt(&(b * b - Atom::num(4) * a * c));
    let plus = -b + &delta;
    let minus = -b - delta;
    let swap = less(&norm(&plus), &norm(&minus));
    let first = if_nonzero_else(&swap, Atom::num(2) * c / &minus, &plus / (Atom::num(2) * a));
    let second = if_nonzero_else(&swap, &minus / (Atom::num(2) * a), Atom::num(2) * c / &plus);
    let linear = if_nonzero(b, -c / b);
    [first, second]
        .map(|r| if_nonzero_else(a, r, linear.clone()))
        .into_iter()
        .enumerate()
        .map(|(index, r)| {
            if_nonzero_else(
                c,
                r,
                if_nonzero(a, if index == 0 { -b / a } else { Atom::num(0) }),
            )
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}

/// `(log((1-y2)/(-y2)) - log((1-y1)/(-y1)))/(y1 - y2)` on the continued
/// sheets, as `log(q)/(1 - q)` of the two quotients.
///
/// When the two sheets of a quotient differ by a winding (`y1` and `y2`
/// infinitesimally on opposite lips of a cut), its `log(q)/(1 - q)` carries
/// `i pi n/(1 - q)`, which for nearly coincident `y1, y2` is large and cancels
/// against other terms of the contour formula; `1 - q` formed from the rounded
/// quotient then costs digits. The windings are therefore taken out and their
/// contribution `i pi (n0 - n1)/(y1 - y2)` uses the difference of the poles.
fn r0(x: [&Atom; 6]) -> Atom {
    let y1 = Continued::from_args([x[0], x[1], x[2]]);
    let y2 = Continued::from_args([x[3], x[4], x[5]]);
    let (q0, n0) = (&y2.complement(0).sheet() / &y1.complement(0).sheet()).split_winding();
    let (q1, n1) = (&y2.complement(1).sheet() / &y1.complement(1).sheet()).split_winding();
    let winding = n0 - n1;
    q0.log_over_one_minus() / &y1.value
        + q1.log_over_one_minus() / (1 - &y1.value)
        + if_nonzero(
            &winding,
            i() * Symbol::PI.to_atom() * &winding / (&y1.value - &y2.value),
        )
}

/// [r0] given `delta = y1 - y2`, formed by the caller without rounding the
/// poles ([r1_line]). With `1 - q0 = delta/y1` and `1 - q1 = -delta/(1 - y1)`,
/// `r0 = (log(q0) - log(q1) + i pi (n0 - n1))/delta`, the logarithms next to
/// one from the quotients' offsets `-delta/y1` and `delta/(1 - y1)`
/// ([SheetAtom::log_given]). A vanishing `delta` keeps [r0]'s limit.
fn r0_line(x: [&Atom; 7]) -> Atom {
    let delta = x[6];
    let y1 = Continued::from_args([x[0], x[1], x[2]]);
    let y2 = Continued::from_args([x[3], x[4], x[5]]);
    let (q0, n0) = (&y2.complement(0).sheet() / &y1.complement(0).sheet()).split_winding();
    let (q1, n1) = (&y2.complement(1).sheet() / &y1.complement(1).sheet()).split_winding();
    let logarithms = q0.log_given(&(-delta / &y1.value))
        - q1.log_given(&(delta / (1 - &y1.value)))
        + i() * Symbol::PI.to_atom() * (n0 - n1);
    if_nonzero_else(delta, logarithms / delta, r0_call(&y1, &y2))
}

fn r1(x: [&Atom; 10]) -> Atom {
    r1_with(x, None)
}

/// [r1] of a line edge, given the product `(y1 - z)(y2 - z)` of the poles'
/// offsets from the root (the last argument, [line_product]). Next to a
/// vanishing Gram determinant a pole lies within 1e-5 of the root of each line
/// edge, and `y - z` formed from the rounded values lost those digits in every
/// quotient and logarithm of the offset. Where the nearer pole's `y - z`
/// cancels by more than 64 (`4096 |y - z|² < |z|²`), its offset is the product
/// over the farther pole's; elsewhere, where the product vanishes (it is
/// not given, [tfun_lines]) and at coincident poles, both are `y - z`, as in
/// [r1]. At an exactly double real pole (`y1 = y2`, the product the square of
/// one offset) the interpreted route (f2adce4) and the generated code
/// (82c87d1) returned 7.6e18 times the value of four relabellings of a box of
/// the streams' sweep, where 09b2c74 and the other routes were non-finite.
fn r1_line(x: [&Atom; 11]) -> Atom {
    r1_with(core::array::from_fn(|index| x[index]), Some(x[10]))
}

fn r1_with(x: [&Atom; 10], product: Option<&Atom>) -> Atom {
    let zz = Continued::from_args([x[0], x[1], x[2]]);
    let yy1 = Continued::from_args([x[3], x[4], x[5]]);
    let yy2 = Continued::from_args([x[6], x[7], x[8]]);
    let (z, y1, y2, fy) = (&zz.value, &yy1.value, &yy2.value, x[9]);
    let oz: Atom = 1 - z;
    let [direct1, direct2] = [y1 - z, y2 - z];
    let (yz1, yz2, flags) = match product {
        None => (direct1, direct2, None),
        Some(product) => {
            let [n1, n2] = [norm(&direct1), norm(&direct2)];
            let nearer_first = less(&n1, &n2);
            let nearer = if_nonzero_else(&nearer_first, n1, n2);
            let cancels = less(&(Atom::num(4096) * nearer), &norm(z))
                * if_nonzero(product, Atom::num(1))
                * if_nonzero(&(y1 - y2), Atom::num(1));
            let flags = [&nearer_first * &cancels, (1 - nearer_first) * cancels];
            // An exactly real `y - z` (real pole and root) stays real: the
            // lip then follows the infinitesimal directions, not the rounded
            // imaginary part of the quotient.
            let from = |near: &Atom, far: &Atom| {
                let offset = if_nonzero(far, product / far);
                if_nonzero_else(&im(near), offset.clone(), re(&offset))
            };
            (
                if_nonzero_else(&flags[0], from(&direct1, &direct2), direct1.clone()),
                if_nonzero_else(&flags[1], from(&direct2, &direct1), direct2),
                Some(flags),
            )
        }
    };
    let offset = |y: &Continued, value: &Atom| Continued {
        value: value.clone(),
        ..y.subtract(&zz)
    };
    let q1 = offset(&yy1, &yz1).sheet();
    let q2 = offset(&yy2, &yz2).sheet();
    let r0_call = |index: usize, y: &Continued, z: &Continued, delta: &Atom| match &flags {
        None => r0_call(y, z),
        Some(flags) => if_nonzero_else(&flags[index], r0_line_call(y, z, delta), r0_call(y, z)),
    };
    let qz = zz.complement(0).sheet();
    let qo = zz.complement(1).sheet();
    let h12 = norm(&(y1 - y2));
    let hz1 = norm(&yz1);
    let hz2 = norm(&yz2);
    let hzz = norm(z);
    let hoz = norm(&oz);
    let zsmall = less(&hzz, &hz1) * less(&hzz, &hz2) * less(&hzz, &hoz);
    let osmall = less(&hoz, &hz1) * less(&hoz, &hz2);
    let product_log = q1.log() + q2.log();
    let yratio = Continued {
        value: (y2 - 1) / y2,
        first: &yy2.first / (y2 * y2),
        second: &yy2.second / (y2 * y2) - &yy2.first * &yy2.first / (y2 * y2 * y2),
    };
    let ylog = yratio.sheet().log();
    let ratio = (&q1 / &q2).log_over_one_minus();
    let zlimit = fy * q1.log() - (&product_log / 2 + &ylog - qo.log()) * &ratio / &yz2;
    let olimit = fy * q1.log() - (-&product_log / 2 + &ylog + qz.log()) * &ratio / &yz2;
    let near_y2 = fy * q1.log() - r0_call(1, &yy2, &zz, &yz2) * &ratio;
    let near_y1 = fy * q2.log() - r0_call(0, &yy1, &zz, &yz1) * (&q2 / &q1).log_over_one_minus();
    let ordinary = (if_nonzero(&yz1, &yz1 * q1.log() * r0_call(0, &yy1, &zz, &yz1))
        - if_nonzero(&yz2, &yz2 * q2.log() * r0_call(1, &yy2, &zz, &yz2)))
        / (y1 - y2);
    let logarithms = if_nonzero_else(
        &zsmall,
        zlimit,
        if_nonzero_else(
            &osmall,
            olimit,
            if_nonzero_else(
                &(le(&h12, &hz2) * le(&hz2, &hz1)),
                near_y2,
                if_nonzero_else(&(le(&h12, &hz1) * le(&hz1, &hz2)), near_y1, ordinary),
            ),
        ),
    );
    // At z=0 the source subtracts two identical Li2(1) values. The divided
    // difference has a logarithmic derivative there, so take the exact zero
    // prefactor limit before evaluating that derivative (likewise for 1-z).
    //
    // A pole next to an end of the interval puts the dilogarithms' arguments
    // `(y - z)/(-z)`, `-z/(y - z)` (and with 1 - z, 1 - y) next to one, where
    // their sheets' `i pi n log(1 - x)` takes the logarithm of the small
    // `1 - x`: `y/z`, `y/(y - z)`, `(1 - y)/(1 - z)` and `(y - 1)/(y - z)`,
    // formed from the poles ([complement_difference]). From the rounded
    // quotient a pole 7e-16 from x = 0 lost 8% of `log(1 - x)`, and a box
    // with an invariant one ulp from a squared mass was off by 1.2 times its
    // value (tests/d0_contour_regressions.rs).
    let dilog_z = if_nonzero_else(
        &zsmall,
        if_nonzero(
            z,
            complement_difference(&(&qz / &q1), &(y1 / &yz1), &(&qz / &q2), &(y2 / &yz2)) * z
                / (&yz1 * &yz2),
        ),
        complement_difference(&(&q1 / &qz), &(y1 / z), &(&q2 / &qz), &(y2 / z)) / z,
    );
    // The first matching branch in r1fun determines which inversion is used.
    let dilog_o = if_nonzero_else(
        &((1 - &zsmall) * &osmall),
        if_nonzero(
            &oz,
            complement_difference(
                &(&qo / &q1),
                &((y1 - 1) / &yz1),
                &(&qo / &q2),
                &((y2 - 1) / &yz2),
            ) * &oz
                / (&yz1 * &yz2),
        ),
        complement_difference(
            &(&q1 / &qo),
            &((1 - y1) / &oz),
            &(&q2 / &qo),
            &((1 - y2) / &oz),
        ) / &oz,
    );
    if_nonzero_else(
        &far_root(y1, y2, &hzz),
        far_root_series(&zz, &yy1, &yy2, fy, &qz),
        logarithms + dilog_z + dilog_o,
    )
}

/// Nonzero when `z` is so far from the interval and the poles that the
/// eight-term [far_root_series] is exact at the working precision:
/// `16 (max(1, |y1|², |y2|²)/|z|²)⁴` leaves 1 unchanged.
///
/// For a far real root (a nearly vanishing linear coefficient, `|z| = 1e15`
/// in a contour box with `p3 = p5`) the dilogarithm terms of [r1] sit next to
/// the cut on different sheets: each carries `i pi n log(1 - q)` with
/// `1 - q = y/z`, so they are O(1/|y|) and cancel, while `1 - q` formed from
/// the rounded quotient costs eps |z| (1e-12 at |z| = 1e3, 1e-6 at 1e9, O(1)
/// at 1e15 in binary64).
fn far_root(y1: &Atom, y2: &Atom, hzz: &Atom) -> Atom {
    let larger = |a: Atom, b: Atom| if_nonzero_else(&less(&a, &b), b, a);
    let scale = larger(larger(Atom::num(1), re(&norm(y1))), re(&norm(y2)));
    let ratio = grouped(&(scale / re(hzz)));
    let squared = &ratio * &ratio;
    let bound = Atom::num(16) * &squared * &squared;
    if_nonzero_else(&(grouped(&(1 + bound)) - 1), Atom::num(0), Atom::num(1))
}

/// [r1] for a far root: `log(x - z) = log(-z) + log(1 - x/z)` on (0, 1) on
/// the continued sheets, so `r1 = fy log(-z) - sum M_n/(n zⁿ)` with the moments
/// `M_n` of `1/((x - y1)(x - y2))` (`M_0 = fy`, `M_1 = y1 fy + L2` with
/// `L2 = log((1 - y2)/(-y2))`, and
/// `M_n = 1/(n - 1) + (y1 + y2) M_{n-1} - y1 y2 M_{n-2}`).
fn far_root_series(z: &Continued, y1: &Continued, y2: &Continued, fy: &Atom, qz: &Q) -> Atom {
    let l2 = y2.complement(1).sheet().log() - y2.complement(0).sheet().log();
    let sum_y = grouped(&(&y1.value + &y2.value));
    let product_y = grouped(&(&y1.value * &y2.value));
    let inverse = grouped(&(Atom::num(1) / &z.value));
    // Each moment expanded as a polynomial in y1 + y2, y1 y2, M_0 and M_1, so
    // that it does not repeat the earlier moments as subexpressions.
    let mut moments = vec![fy.clone(), grouped(&(&y1.value * fy + l2))];
    for n in 2..=8i64 {
        let next = Atom::num(1) / Atom::num(n - 1) + &sum_y * &moments[n as usize - 1]
            - &product_y * &moments[n as usize - 2];
        moments.push(next.expand());
    }
    // Horner in 1/z: sum_{n=1}^{8} M_n/n z^-n.
    let mut sum = &moments[8] / Atom::num(8);
    for n in (1..8).rev() {
        sum = grouped(&(&moments[n] / Atom::num(n as i64) + &inverse * &sum));
    }
    fy * qz.log() - &inverse * &sum
}

fn s3_double_real([y1, y2, z, a]: [&Atom; 4]) -> Atom {
    // For a>0 and a real double root inside (0,1), the physical limit is
    // log[a*(x-z)^2-i0] = log(a)+2*log|x-z| almost everywhere. The outer
    // poles are off the closed integration interval, so its kernel is
    // bounded and this integrable logarithmic limit can be taken directly.
    // Splitting the integral at z gives J(y)=integral log|x-z|/(x-y) dx.
    // Exact divided differences retain the removable y1=y2 limit without
    // assigning incompatible epsilon and sqrt(epsilon) root directions.
    let a1 = z - y1;
    let a2 = z - y2;
    let oz: Atom = 1 - z;
    let q01 = &Q::upper(-y1) / &Q::upper(a1.clone());
    let q02 = &Q::upper(-y2) / &Q::upper(a2.clone());
    let q11 = &Q::upper(1 - y1) / &Q::upper(a1.clone());
    let q12 = &Q::upper(1 - y2) / &Q::upper(a2.clone());
    let jdd = (z * dd(&q01, &q02) + &oz * dd(&q11, &q12)
        - z * z.log() * (&q01 / &q02).log_over_one_minus() / q02.value()
        - &oz * oz.log() * (&q11 / &q12).log_over_one_minus() / q12.value())
        / (&a1 * &a2);
    let pole = |y: &Atom| Continued {
        value: y.clone(),
        first: Atom::num(0),
        second: Atom::num(0),
    };
    r0_call(&pole(y1), &pole(y2)) * a.log() + Atom::num(2) * jdd
}

fn outside_closed_interval(y: &Atom) -> Atom {
    if_nonzero_else(
        &im(y),
        Atom::num(1),
        negative(&re(y)) + negative(&(1 - re(y))),
    )
}

/// The poles `y1, y2` of `s3`'s arguments mapped by `y -> (d + y)/e`, `r0`
/// of them and the sign that continues the coefficients' logarithms.
fn s3_poles(y: [&Atom; 6], [d, e, b, c, db, dc]: [&Atom; 6]) -> (Continued, Continued, Atom, Atom) {
    let mut y1 = Continued::from_args([y[0], y[1], y[2]]);
    let mut y2 = Continued::from_args([y[3], y[4], y[5]]);
    let minus_one = Atom::num(-1);
    let simc = imaginary_sign(
        c,
        &imaginary_sign(dc, &imaginary_sign(b, &imaginary_sign(db, &minus_one))),
    );
    for y in [&mut y1, &mut y2] {
        y.value = (d + &y.value) / e;
        y.first = &y.first / e;
        y.second = &y.second / e;
    }
    let fy = r0_call(&y1, &y2);
    (y1, y2, fy, simc)
}

/// `s3` of a linear `b x + c` (`far = b + c`) against the mapped poles, with
/// the product of the poles' offsets from its root when the caller has it.
fn s3_linear_term(
    [y1, y2]: [&Continued; 2],
    fy: &Atom,
    simc: &Atom,
    [b, c, db, dc, far]: [&Atom; 5],
    product: Option<&Atom>,
) -> Atom {
    let zero = Atom::num(0);
    let mut z = continued_root(
        if_nonzero_else(far, -c / b, Atom::num(1)),
        &zero,
        b,
        db,
        dc,
        &zero,
    );
    let known = any_nonzero([im(&z.value), im(&z.first), im(&z.second)]);
    z.second = if_nonzero_else(
        &known,
        z.second.clone(),
        -i() * simc * sign_nonnegative(&re(b)) * Symbol::ABS.call((re(&z.value),)),
    );
    let sz = z.sign();
    let sb = imaginary_sign(db, simc);
    let r1 = match product {
        None => r1_call(&z, y1, y2, fy),
        Some(product) => r1_line_call(&z, y1, y2, fy, product),
    };
    fy * (Q::with_sign(b.clone(), sb.clone()).log() + eta(b, &sb, &(-&z.value), &(-&sz), c, simc))
        + r1
}

/// [s3] of a line edge of [tfun_selected] (`a = db = 0`), given the product
/// `(y1 - z)(y2 - z)` of the mapped poles' offsets from the root `z = -c/b`
/// (1 for `far = 0`), formed by [line_product].
fn s3_line(x: [&Atom; 13]) -> Atom {
    let zero = Atom::num(0);
    let [d, e, b, c, dc, far, product] = [x[6], x[7], x[8], x[9], x[10], x[11], x[12]];
    let (y1, y2, fy, simc) = s3_poles(
        [x[0], x[1], x[2], x[3], x[4], x[5]],
        [d, e, b, c, &zero, dc],
    );
    let linear = s3_linear_term(
        [&y1, &y2],
        &fy,
        &simc,
        [b, c, &zero, dc, far],
        Some(product),
    );
    let constant = if_nonzero(c, Q::with_sign(c.clone(), simc).log() * &fy);
    if_nonzero(e, if_nonzero_else(b, linear, constant) / e)
}

/// The integral over `x` in (0, 1) of `log(a x² + b x + c)` against the poles
/// `y1, y2` mapped by `y -> (d + y)/e`.
///
/// Each quadratic (`a` nonzero) is the Feynman polynomial of a two-point
/// subgraph, which the caller passes from the invariants (see [Edge]): its
/// discriminant `lambda = λ(p, near, far)` ([kallen]), `p` and the masses
/// `near` and `far` at x = 0 and x = 1. An exact double root is then decided exactly,
/// next to one the roots keep their separation `sqrt(lambda)/a` and the signs
/// of their slopes `2 a z + b = ±sqrt(lambda)`, and a massless end is an
/// exact root. A linear `b x + c` passes `far = b + c`.
fn s3(x: [&Atom; 17]) -> Atom {
    let [d, e, a, b, c, db, dc] = [x[6], x[7], x[8], x[9], x[10], x[11], x[12]];
    let [lambda, p, near, far] = [x[13], x[14], x[15], x[16]];
    let (y1, y2, fy, simc) = s3_poles([x[0], x[1], x[2], x[3], x[4], x[5]], [d, e, b, c, db, dc]);
    let ra = re(a);
    let rea = sign_nonnegative(&ra);
    // solabc's roots with sqrt(lambda): the first has slope 2 a z + b equal
    // to +sqrt(lambda), the second -sqrt(lambda). The continuation of a root
    // follows the sign of its slope, which next to a double root the rounded
    // coefficients no longer resolve.
    let delta = principal_sqrt(lambda);
    let plus = -b + &delta;
    let minus = -b - &delta;
    let swap = less(&norm(&plus), &norm(&minus));
    let two_a = Atom::num(2) * &ra;
    let upper = if_nonzero_else(&swap, Atom::num(2) * c / &minus, &plus / &two_a);
    let lower = if_nonzero_else(&swap, &minus / &two_a, Atom::num(2) * c / &plus);
    let zero = Atom::num(0);
    let general1 = continued_root_with_slope(upper, delta.clone(), a, db, dc, &zero);
    let general2 = continued_root_with_slope(lower, -&delta, a, db, dc, &zero);
    // A massless end (near = 0 at x = 0, or far = 0 at x = 1) is a root, and
    // the other root lies `slope/p` from it, with the slope 2 a z + b there
    // (p - far, or near - p), exact from the invariants. Next to a
    // pseudo-threshold that root rounds onto the end, and the rounded
    // coefficients no longer give the sign of its slope, nor its first-order
    // numerator db z + dc, which then cancels to rounding (db end + dc is the
    // end value's direction, zero for a massless line): both roots continue
    // from the end value instead.
    let at_end = |general: Atom, zero_near: Atom, zero_far: Atom| {
        if_nonzero_else(far, if_nonzero_else(near, general, zero_near), zero_far)
    };
    let end = at_end(zero.clone(), zero.clone(), Atom::num(1));
    let slope = if_nonzero_else(far, p - far, near - p);
    let direction = if_nonzero_else(far, dc.clone(), db + dc);
    let end_first = &direction / &slope;
    let end_root = Continued {
        value: end.clone(),
        first: end_first.clone(),
        second: (a * &end_first * &end_first + db * &end_first) / &slope,
    };
    let near_first = -db / a - &direction / &slope;
    let near_root = Continued {
        value: &end + &slope / re(p),
        first: near_first.clone(),
        second: &near_first * a * &direction / (&slope * &slope),
    };
    let pick = |general: Continued, special: Continued| Continued {
        value: at_end(general.value, special.value.clone(), special.value),
        first: at_end(general.first, special.first.clone(), special.first),
        second: at_end(general.second, special.second.clone(), special.second),
    };
    let mut z1 = pick(general1, end_root);
    let mut z2 = pick(general2, near_root);
    // An exact double root has no first-order continuation (2 a z + b = 0).
    // Infinitesimal widths split it like sqrt(eps): one root moves up and
    // the other down, which dominates every first-order shift. Only the
    // signs are used, so encode the split as a first-order direction that
    // exceeds those of the poles y1, y2 it is combined with.
    let double = if_nonzero_else(lambda, Atom::num(0), if_nonzero(a, Atom::num(1)));
    let dominant = Atom::num(2)
        * (Atom::num(1)
            + [&y1.first, &y1.second, &y2.first, &y2.second]
                .into_iter()
                .fold(Atom::num(0), |sum, v| sum + Symbol::ABS.call((v,))));
    // Both roots are -b/(2a) there, or the massless end.
    let z0 = at_end(-b / (Atom::num(2) * a), Atom::num(0), Atom::num(1));
    let centre = z0.clone();
    for (z, direction) in [(&mut z1, Atom::num(1)), (&mut z2, Atom::num(-1))] {
        z.value = if_nonzero_else(&double, centre.clone(), z.value.clone());
        z.first = if_nonzero_else(&double, i() * direction * &dominant, z.first.clone());
        z.second = if_nonzero_else(&double, Atom::num(0), z.second.clone());
    }
    let fallback1 = &simc * &rea * sign_nonnegative(&re(&z2.value));
    let fallback2 = &simc * &rea * sign_nonnegative(&re(&z1.value));
    for (z, fallback) in [(&mut z1, fallback1), (&mut z2, fallback2)] {
        let known = any_nonzero([im(&z.value), im(&z.first), im(&z.second)]);
        z.second = if_nonzero_else(
            &known,
            z.second.clone(),
            i() * fallback * Symbol::ABS.call((re(&z.value),)),
        );
    }
    let s1 = z1.sign();
    let s2 = z2.sign();
    let ordinary_quadratic = &fy
        * (Q::with_sign(a.clone(), simc.clone()).log()
            + eta(
                &(-&z1.value),
                &(-&s1),
                &(-&z2.value),
                &(-&s2),
                &Atom::num(0),
                &(&simc * &rea),
            ))
        + r1_call(&z1, &y1, &y2, &fy)
        + r1_call(&z2, &y1, &y2, &fy);
    let real_coefficients = 1 - any_nonzero([im(a), im(b), im(c)]);
    let positive_a = negative(&(-re(a)));
    let interior_z = negative(&(-re(&z0))) * negative(&(re(&z0) - 1));
    let regular_poles = outside_closed_interval(&y1.value) * outside_closed_interval(&y2.value);
    let double_root = if_nonzero_else(
        lambda,
        Atom::num(0),
        real_coefficients * positive_a * interior_z * regular_poles,
    );
    let quadratic = if_nonzero_else(
        &double_root,
        call("s3_double_real", [&y1.value, &y2.value, &z0, a]),
        ordinary_quadratic,
    );
    let linear = s3_linear_term([&y1, &y2], &fy, &simc, [b, c, db, dc, far], None);
    let constant = if_nonzero(c, Q::with_sign(c.clone(), simc).log() * &fy);
    if_nonzero(
        e,
        if_nonzero_else(a, quadratic, if_nonzero_else(b, linear, constant)) / e,
    )
}

/// Residues of the poles enclosed by the contour (OneLOop's `plnr`).
///
/// `real` (the last argument) is nonzero when the logarithm arguments
/// `a/(b y + c)` are exactly real; see [real_intersections_absent].
fn planar(x: [&Atom; 13]) -> Atom {
    let y1 = Continued::from_args([x[0], x[1], x[2]]);
    let y2 = Continued::from_args([x[3], x[4], x[5]]);
    let [p1, p2, a, b, c, dc, real] = [x[6], x[7], x[8], x[9], x[10], x[11], x[12]];
    let x1 = b * &y1.value + c;
    let x2 = b * &y2.value + c;
    let sheet = |y: &Continued, value: &Atom| {
        let first = b * &y.first + dc;
        let second = b * &y.second;
        let ratio = a / value;
        Continued {
            value: if_nonzero_else(real, re(&ratio), ratio),
            first: -a * &first / (value * value),
            second: a * &first * &first / (value * value * value) - a * second / (value * value),
        }
        .sheet()
    };
    let q1 = sheet(&y1, &x1);
    let q2 = sheet(&y2, &x2);
    let both = (&q2 / &q1).log_over_one_minus() * b / x2;
    let only_first = q1.log() / (&y1.value - &y2.value);
    let only_second = q2.log() / (&y2.value - &y1.value);
    Atom::num(2)
        * i()
        * Symbol::PI.to_atom()
        * if_nonzero_else(
            p1,
            if_nonzero_else(p2, both, only_first),
            if_nonzero(p2, only_second),
        )
}

/// `ac` is `a + c`, formed by the caller from the invariants it stands for
/// (see [coefficient_sum]). The edge quadratics `(a + c, dpe, f)` and
/// `(a, d, f)` have the discriminants `lambda_ac`, `lambda_a` and the values
/// `far_ac`, `far_a` at x = 1 (see [Edge]); `gh`, `gj` and `ghj` are the
/// line's `g + h`, `g + j` and `g + h + j` (see [tfun]).
fn t1(
    [
        a,
        c,
        g,
        h,
        d,
        e,
        f,
        j,
        dpe,
        dd,
        de,
        df,
        dj,
        ddpe,
        ac,
        lambda_ac,
        far_ac,
        lambda_a,
        far_a,
        gh,
        gj,
        ghj,
    ]: [&Atom; 22],
) -> Atom {
    let k = h * a - c * g;
    let l = h * d - c * j - e * g;
    let n = h * f - e * j;
    let dl = h * dd - c * dj - de * g;
    let dn = h * df - de * j - e * dj;
    let dn2 = -de * dj - i() * Symbol::ABS.call((re(f),)) * h;
    let [y1, y2] = roots(&k, &l, &n).map(|y| continued_root(y, &k, &l, &dl, &dn, &dn2));
    let zero = Atom::num(0);
    let one = Atom::num(1);
    let sum = coefficient_sum(&[a, c], ac);
    let edge_ac = Edge {
        lambda: lambda_ac,
        p: ac,
        near: f,
        far: far_ac,
    };
    let edge_a = Edge {
        lambda: lambda_a,
        p: a,
        near: f,
        far: far_a,
    };
    (-s3_call(&y1, &y2, &zero, &one, &sum, dpe, f, ddpe, df, edge_ac)
        + s3_linear_call(&y1, &y2, &zero, &one, gh, j, dj, ghj)
        - s3_linear_call(&y1, &y2, &zero, &one, g, j, dj, gj)
        + s3_call(&y1, &y2, &zero, &one, a, d, f, dd, df, edge_a))
        / k
}

/// The edge quadratics have the discriminants `lambda_*` (in the order of the
/// terms) and the values `far_ae` (first, last) and `far_gj` (second,
/// third) at x = 1; `gh` and `ac` are `g + h` and `a + c` (see [Edge]).
fn t13(
    [
        a,
        c,
        g,
        h,
        d,
        e,
        f,
        j,
        dpe,
        dpj,
        dpf,
        dd,
        de,
        df,
        dj,
        ddpe,
        ddpj,
        ddpf,
        lambda_ae,
        lambda_gj,
        lambda_gh,
        lambda_ac,
        far_ae,
        far_gj,
        gh,
        ac,
    ]: [&Atom; 26],
) -> Atom {
    let k = h * a - c * g;
    let l = a * d + h * e - d * g - c * j;
    let n = d * (e - j) + (h - c) * f;
    let dl = (a - g) * dd + h * de - c * dj;
    let dn = dd * (e - j) + d * (de - dj) + (h - c) * df;
    let dn2 = dd * (de - dj) - i() * Symbol::ABS.call((re(f),)) * (h - c);
    let [y1, y2] = roots(&k, &l, &n).map(|y| continued_root(y, &k, &l, &dl, &dn, &dn2));
    let zero = Atom::num(0);
    let one = Atom::num(1);
    let edge = |lambda, p, near, far| Edge {
        lambda,
        p,
        near,
        far,
    };
    (-s3_call(&y1, &y2, &zero, &one, a, &(e + c), dpf, de, ddpf, edge(lambda_ae, a, dpf, far_ae))
        + s3_call(&y1, &y2, &zero, &one, g, &(j + h), dpf, dj, ddpf, edge(lambda_gj, g, dpf, far_gj))
        - s3_call(&y1, &y2, &zero, &one, &(g + h), dpj, f, ddpj, df, edge(lambda_gh, gh, f, far_gj))
        + s3_call(&y1, &y2, &zero, &one, &(a + c), dpe, f, ddpe, df, edge(lambda_ac, ac, f, far_ae)))
        / k
}

fn coordinates(beta: &Atom, y: &Atom) -> [Atom; 2] {
    // Barycentric coordinates in the triangle (0, -beta, 1-beta).
    // This is the algebraically reduced Gram-system solution in tfun.
    let sum = -im(y) / im(beta);
    let x2 = re(y) + re(beta) * &sum;
    let x1 = sum - &x2;
    [x1, x2]
}

fn continued_nonnegative(value: &Atom, first: &Atom, second: &Atom) -> Atom {
    if_nonzero_else(
        value,
        1 - negative(value),
        if_nonzero_else(first, 1 - negative(first), 1 - negative(second)),
    )
}

fn inside([beta, y, dy, ddy]: [&Atom; 4]) -> Atom {
    let [x1, x2] = coordinates(beta, y);
    let [d1, d2] = coordinates(beta, dy);
    let [dd1, dd2] = coordinates(beta, ddy);
    if_nonzero(
        &im(beta),
        continued_nonnegative(&x1, &d1, &dd1)
            * continued_nonnegative(&x2, &d2, &dd2)
            * continued_nonnegative(&(1 - x1 - x2), &(-d1 - d2), &(-dd1 - dd2)),
    )
}

fn distance([beta, y, dy, ddy]: [&Atom; 4]) -> Atom {
    let [x1, x2] = coordinates(beta, y);
    if_nonzero_else(
        &call("inside", [beta, y, dy, ddy]),
        min(&x1, &x2),
        Atom::num(1),
    )
}

struct BetaData {
    a: Atom,
    b: Atom,
    c: Atom,
    dc: Atom,
    k: Atom,
    y: [Continued; 2],
}

/// `w` is `beta - 1`, `dpe` and `ddpe` are `d + e` and `dd + de`.
///
/// A small `a + b + c` (an invariant) puts a root `beta` of
/// `b beta² + c beta + a` next to 1, where `d + beta e` and the like cancel
/// (by 800 in a three-mass box with s23 = 0.0049 next to masses 0.013, which
/// cost binary64 3e-11 in its poles). Nearer 1 than 0 the coefficients are
/// therefore expanded about beta = 1.
///
/// `k`, where given, replaces `K = b a1 - h b1` ([chart_k]).
fn beta_data(
    beta: &Atom,
    w: &Atom,
    b: &Atom,
    c: &Atom,
    g: &Atom,
    h: &Atom,
    d: &Atom,
    e: &Atom,
    f: &Atom,
    j: &Atom,
    dpe: &Atom,
    dd: &Atom,
    de: &Atom,
    df: &Atom,
    dj: &Atom,
    ddpe: &Atom,
    k: Option<&Atom>,
) -> BetaData {
    let near_one = less(&norm(w), &norm(beta));
    let expanded = |at_one: Atom, at_zero: Atom| if_nonzero_else(&near_one, at_one, at_zero);
    let [a1, b1, c1, l, n] = beta_polynomial(beta, w, b, c, g, h, d, e, f, j, dpe);
    let k = k.cloned().unwrap_or_else(|| b * &a1 - h * &b1);
    let dc1 = expanded(ddpe + w * de, dd + beta * de);
    let dl = de * &a1 - h * &dc1 - dj * &b1;
    let dn = df * &a1 - dj * &c1 - j * &dc1;
    let dn2 = -dj * &dc1 - i() * Symbol::ABS.call((re(f),)) * &a1;
    let y = roots(&k, &l, &n).map(|y| continued_root(y, &k, &l, &dl, &dn, &dn2));
    BetaData {
        a: a1,
        b: b1,
        c: c1,
        dc: dc1,
        k,
        y,
    }
}

/// `[a1, b1, c1, L, N]` of [beta_data].
fn beta_polynomial(
    beta: &Atom,
    w: &Atom,
    b: &Atom,
    c: &Atom,
    g: &Atom,
    h: &Atom,
    d: &Atom,
    e: &Atom,
    f: &Atom,
    j: &Atom,
    dpe: &Atom,
) -> [Atom; 5] {
    let near_one = less(&norm(w), &norm(beta));
    let expanded = |at_one: Atom, at_zero: Atom| if_nonzero_else(&near_one, at_one, at_zero);
    let a1 = expanded(grouped(&(g + h)) + w * h, g + beta * h);
    let b1 = expanded(
        grouped(&(c + Atom::num(2) * b)) + Atom::num(2) * w * b,
        c + Atom::num(2) * beta * b,
    );
    let c1 = expanded(dpe + w * e, d + beta * e);
    let l = e * &a1 - h * &c1 - j * &b1;
    let n = f * &a1 - j * &c1;
    [a1, b1, c1, l, n]
}

/// The `K` of [tfun_selected]'s chart: `b a1 - h b1` of [beta_data], or,
/// where that rounds to zero, the one that puts the far pole `-L/K` at
/// `2^52 (1 + |N/L|²)`.
///
/// `K = b (g + beta' h)` with the other root `beta'` vanishes with the Gram
/// determinant: one pole goes to infinity, and the chart's terms vanish with
/// `K` (they are divided by it). A `K` that rounds to zero (0/0 and NaN at
/// 09b2c74) is replaced so that the far pole lies 2^52 beyond the near pole
/// `-N/L` and the triangle; the terms are then within O(2^-52) of their
/// limit. The contour selection ([tfun]) takes the roots of the unreplaced
/// `K`, where a vanishing `K` repeats the near pole.
fn chart_k(
    beta: &Atom,
    w: &Atom,
    b: &Atom,
    c: &Atom,
    g: &Atom,
    h: &Atom,
    d: &Atom,
    e: &Atom,
    f: &Atom,
    j: &Atom,
    dpe: &Atom,
) -> Atom {
    let [a1, b1, _, l, n] = beta_polynomial(beta, w, b, c, g, h, d, e, f, j, dpe);
    let k = b * &a1 - h * &b1;
    if_nonzero_else(
        &k,
        k.clone(),
        if_nonzero(&l, -&l / (Atom::num(1u64 << 52) * (1 + norm(&(&n / &l))))),
    )
}

/// The product `(Y1 - z)(Y2 - z)` of the offsets of the mapped poles
/// `Y = (d + y)/e` of [beta_data] from the root `z` of a line edge `b x + c`
/// (the `[d, e, b, c, far]` of its [s3_line] term), where `h u + j = ratio a1`
/// at the preimage `u = e z - d`.
///
/// The poles are the roots of `K y² + L y + N = a1 Q(y) - (h y + j) B(y)`,
/// with `a1 = g + beta h`, `Q(y) = b y² + e y + f` and `B(y) = b1 y + c1`; at
/// `u` the quadratic is `a1 (Q(u) - ratio B(u))`, exactly proportional to
/// `a1`, so that `(y1 - u)(y2 - u) = a1 (Q(u) - ratio B(u))/K` (Vieta). Next
/// to a vanishing Gram determinant `a1` is small (the line is nearly parallel
/// to the contour direction `(1, beta)`), and one pole lies next to the root
/// of every line edge: at the NinjaRS gg -> HHH D0 within 5e-6 to 7e-5 of
/// roots of order one, with `a1/g = 1.4e-5`. `y - z` formed from the rounded
/// pole and root lost those digits, and each edge lost them differently: the
/// edges' terms then failed to cancel by 1e-12 of themselves (D0 3.6e-11 in
/// binary64). [r1_line] takes the nearer pole's offset from this product and
/// the farther pole's direct offset, with the poles' own `a1` (passed as `a1`,
/// see [tfun_selected]) and `K`, so that every edge sees the offsets from the
/// line the rounded `a1` describes.
///
/// `Q(u) - ratio B(u)` is the conic at the line's intersection with the edge.
/// Where that intersection lies on the conic to within the square root of
/// the rounding (`2^52 |Q - ratio B|²` below `conic_norm`, the sum of the
/// conic coefficients' norms), the conic's own root on the edge lies next to
/// the line's and to the pole, and the line's and the conic's terms on the
/// edge cancel with the errors of their common pole: both keep `y - z`
/// (the product is zero, and `r1_line` keeps `y - z`). With the product there,
/// a box with a massless line next to an on-shell corner lost 1e-5 at 128
/// bits where `y - z` kept 4e-11.
fn line_product(
    data: &BetaData,
    a1: &Atom,
    conic_norm: &Atom,
    [d, e, b, c, far]: [&Atom; 5],
    ratio: &Atom,
    [qb, qe, qf]: [&Atom; 3],
) -> Atom {
    let u = if_nonzero_else(far, -(e * c + d * b) / b, e - d);
    let quadratic = (qb * &u + qe) * &u + qf;
    let value = quadratic - ratio * (&data.b * &u + &data.c);
    let off_conic = less(conic_norm, &(Atom::num(1u64 << 52) * norm(&value)));
    if_nonzero(
        e,
        if_nonzero_else(&off_conic, a1 * value / (&data.k * e * e), Atom::num(0)),
    )
}

/// Nonzero when the real line `g x + h y + j = 0` meets the real conic
/// `a x² + b y² + c x y + d x + e y + f = 0` at complex conjugate points.
///
/// The residue logarithms of `tfun` take `q = a1/(b1 y + c1)` at a pole `y`,
/// where `a1 = L(P + s w)/s` and `b1 y + c1 = Q(P + s w)/s` for the
/// intersection point `P` belonging to `y`, the direction `w = (1, beta)` and
/// every `s` (the conic has no `s²` term along `w`). For a non-real `P` and
/// non-real `beta` one `s` makes `P + s w` real, so `q = L/Q` there is real.
/// For exactly real masses `q` is then exactly real, while its rounded
/// imaginary part is noise that would otherwise choose its lip (and add
/// `±2 pi i` to a residue). The lip must follow the infinitesimal widths.
/// The intersections are complex iff `n^T adj(M) n > 0`, with `M` the conic's
/// symmetric matrix and `n = (g, h, j)`.
fn real_intersections_absent(coefficients: [&Atom; 9]) -> Atom {
    let real = 1 - any_nonzero(coefficients.map(im));
    let [a, b, c, g, h, d, e, f, j] = coefficients.map(re);
    let quarter = |x: Atom| x / 4;
    let adj11 = &b * &f - quarter(&e * &e);
    let adj22 = &a * &f - quarter(&d * &d);
    let adj33 = &a * &b - quarter(&c * &c);
    let adj12 = quarter(&d * &e) - &c * &f / 2;
    let adj13 = quarter(&c * &e) - &b * &d / 2;
    let adj23 = quarter(&c * &d) - &a * &e / 2;
    let delta = &g * &g * adj11
        + &h * &h * adj22
        + &j * &j * adj33
        + Atom::num(2) * (&g * &h * adj12 + &g * &j * adj13 + &h * &j * adj23);
    real * negative(&grouped(&-delta))
}

/// Nonzero where the line edge `index` (x = 1, y = x, y = 0, as in
/// [tfun_lines]) of a parallel chart takes its product (2 where the product
/// is also probed, see below): the arguments are
/// `index`, the chart's poles, `beta, w, g, h, j`, the line's sums `gh, gj,
/// ghj`, the three edges' products ([line_product]), `a1` and `scale =
/// |g|² + |beta h|²`.
///
/// The products describe the poles' offsets from the line edges' roots
/// where those offsets scale with `a1` (the line nearly parallel to the
/// contour direction, next to a vanishing Gram determinant): there the
/// rounded `y - z` of every edge loses digits to its own subtraction, while
/// the products share the poles' `a1`. Elsewhere the products break the
/// cancellations of the rounded pole that `y - z` keeps, and binary64 loses
/// more with them than without (rows of the verifier's 32,600-row sweep at
/// 1e-12 to 1e-5 where 09b2c74, `y - z` everywhere, kept 1.5e-13 to 7e-7).
/// An edge whose nearer pole's `y - z` cancels by more than 64 (the edges
/// that [r1_line] gives the product) is tested, and the chart takes no
/// product where such an edge has none (its intersection with the line lies
/// on the conic, [line_product], and keeps `y - z`: v1 5635, whose y = 0 edge
/// lies within 4e-16 of its pole). An edge whose pole lies within 8
/// offsets of an end of the edge (`min(|y|², |1 - y|²) < 64 |y - z|²`) keeps
/// `y - z`: there [r1_line] also takes the pole's distance to that end from
/// the rounded pole, which the product's offset no longer matches (v1 19829,
/// 2.4e-7 from the corner where the edges x = 1 and y = x meet, with
/// offsets 5.4e-8 and 5.6e-8; v1 1255, 13854 and g1 6075).
///
/// A kept product whose offset `|y - z|/|z|` is more than 2^10 times
/// `a1/sqrt(scale)` (an offset that `a1` does not account for) gives 2: the
/// chart is then probed ([tfun_selected], [tfun_edge_units]) as if that
/// product disagreed with its pole by at least 128 units, which takes `y - z`
/// everywhere where the products move the chart more under the pole's shift
/// (v1 8793 and 18833, 460 to 5700 times, as 09b2c74). Dropping every
/// product of such a chart instead (59b96e5) lost the chart's digits where
/// the products carry them: reflected, fs 15926 of the second verifier's
/// 36,840-row sweep (one edge's offset 1.3e4 times `a1`, the two others 1.4
/// and 2.6e-4 times) was wrong by 4.9e-6 to 1.4e-5, where every product keeps
/// 2.4e-14.
///
/// An `a1` within 16 units of its rounding (`2^96 |a1|² < scale`, the Gram
/// determinant at the rounding level) has too few digits to compare the
/// offsets with, and a product may replace a rounded `y - z` that vanishes
/// (whose logarithms would be `log 0`, [tfun_lines]): there every product is
/// kept (SymJIT lost v1 12721, whose flat charts have `a1` at 6e-18 to 6e-17
/// of their terms, by its value without them, as at 09b2c74).
fn products_kept(x: [&Atom; 16]) -> Atom {
    let [index, y1, y2, beta, w, g, h, j, gh, gj, ghj] =
        [x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10]];
    let [a1, scale] = [x[14], x[15]];
    let [zero, one] = [Atom::num(0), Atom::num(1)];
    let [minus_w, minus_beta] = [-w, -beta];
    let edges: [[&Atom; 5]; 3] = [
        [beta, &one, h, gj, ghj],
        [&zero, &minus_w, gh, j, ghj],
        [&zero, &minus_beta, g, j, gj],
    ];
    let mut chart = Atom::num(1);
    let mut corners = vec![];
    let mut probes = vec![];
    let mut vanishing = Atom::num(1);
    for ([d, e, b, c, far], product) in edges.into_iter().zip([x[11], x[12], x[13]]) {
        let valid = if_nonzero(e, Atom::num(1)) * if_nonzero(b, Atom::num(1));
        let [p1, p2] = [y1, y2].map(|y| (d + y) / e);
        let z = if_nonzero_else(far, -c / b, Atom::num(1));
        let [n1, n2] = [norm(&(&p1 - &z)), norm(&(&p2 - &z))];
        let first = less(&n1, &n2);
        let nearer = if_nonzero_else(&first, n1, n2);
        let pole = if_nonzero_else(&first, p1.clone(), p2.clone());
        let hzz = norm(&z);
        let cancels = valid
            * less(&(Atom::num(4096) * &nearer), &hzz)
            * if_nonzero(&(&p1 - &p2), Atom::num(1));
        let unscaled = less(&(Atom::num(1u64 << 20) * &hzz * norm(a1)), &(&nearer * scale));
        let fails = &cancels * (1 - if_nonzero(product, Atom::num(1)));
        chart = chart * (1 - fails);
        vanishing = vanishing
            * (1 - &cancels * if_nonzero_else(&nearer, Atom::num(0), Atom::num(1)) * if_nonzero(product, Atom::num(1)));
        let end = min(&norm(&pole), &norm(&(&pole - 1)));
        corners.push(&cancels * less(&end, &(Atom::num(64) * nearer)));
        probes.push(cancels * unscaled);
    }
    let by_index = |values: &[Atom]| {
        if_nonzero_else(
            index,
            if_nonzero_else(&(index - 1), values[2].clone(), values[1].clone()),
            values[0].clone(),
        )
    };
    let [corner, probe] = [by_index(&corners), by_index(&probes)];
    // An `a1` within 16 units of its rounding, or a product where the
    // rounded `y - z` vanishes, keeps every product.
    let unresolved = less(&(Atom::num(symbolica::domains::integer::Integer::from(2).pow(96)) * norm(a1)), scale);
    let all = 1 - (1 - unresolved) * vanishing;
    if_nonzero_else(&all, Atom::num(1), chart * (1 - corner) * (1 + probe))
}

/// The chart of [tfun] for the contour `beta` (with `w = beta - 1` and the
/// chart's `k`), with the line products only where `products` (and the line
/// is nearly parallel to the contour direction; 1: those of [products_kept],
/// 2: all of them) and the poles moved by `shift1`, `shift2`
/// ([tfun_selected]).
fn tfun_lines(
    [
        a,
        b,
        c,
        g,
        h,
        d,
        e,
        f,
        j,
        dpe,
        dpf,
        beta,
        w,
        dd,
        de,
        df,
        dj,
        ddpe,
        ddpf,
        abc,
        lambda_a,
        lambda_abc,
        lambda_bc,
        q10,
        q11,
        gh,
        gj,
        ghj,
        k,
        products,
        shift1,
        shift2,
    ]: [&Atom; 32],
) -> Atom {
    let mut data = beta_data(beta, w, b, c, g, h, d, e, f, j, dpe, dd, de, df, dj, ddpe, Some(k));
    for (y, shift) in data.y.iter_mut().zip([shift1, shift2]) {
        y.value = &y.value + shift;
    }
    let [y1, y2] = &data.y;
    let zero = Atom::num(0);
    let one = Atom::num(1);
    let p1 = call("inside", [beta, &y1.value, &y1.first, &y1.second]);
    let p2 = call("inside", [beta, &y2.value, &y2.first, &y2.second]);
    let sum = coefficient_sum(&[a, b, c], abc);
    let edge = |lambda, p, near, far| Edge {
        lambda,
        p,
        near,
        far,
    };
    // The line edges x = 1, y = x and y = 0 as the `[d, e, b, c, far]` of
    // their s3 terms, and h u + j over a1 at their roots' preimages u.
    //
    // The products replace the poles' offsets from the roots only for a line
    // nearly parallel to the contour direction (1, beta), where a1 = g + beta h
    // cancels by more than 64: there a pole lies next to every edge's root
    // because a1 is small (the vanishing Gram determinant). A pole also lies
    // next to an edge's root where the conic nearly vanishes at the line's
    // intersection with the edge, with an `a1` of the size of its terms;
    // there the product's polynomial cancels at least as much as `y - z` (by
    // 4e5 against 1e4 in a three-mass box), and the line's term cancels
    // against the conic's term on the same edge, whose root lies next to the
    // same pole: both keep `y - z` from the same rounded pole (with the
    // product, three-mass boxes lost 9e-12 to 3.9e-10 where `y - z` kept
    // 2.9e-12). Nor is the product taken on an edge whose intersection with
    // the line lies on the conic ([line_product]). A zero product keeps
    // `y - z` ([r1_line]).
    //
    // An `a1` that rounds to zero (the Gram determinant at the rounding
    // level) is taken at its rounding level, `2^-52 g`: the poles' offsets
    // are proportional to it, and the logarithms of their common factor
    // cancel between the edges, where exact zero offsets left `log 0` in
    // each edge's term (O(1) errors on single routes in a few such boxes).
    let conic_norm = [a, b, c, d, e, f]
        .into_iter()
        .fold(Atom::num(0), |sum, coefficient| sum + norm(coefficient));
    let scale = norm(g) + norm(&(beta * h));
    let parallel = if_nonzero(products, Atom::num(1)) * less(&(Atom::num(4096) * norm(&data.a)), &scale);
    let a1 = if_nonzero_else(&data.a, data.a.clone(), g / Atom::num(1u64 << 52));
    let product = |line: [&Atom; 5], ratio: Atom| {
        if_nonzero_else(
            &parallel,
            line_product(&data, &a1, &conic_norm, line, &ratio, [b, e, f]),
            Atom::num(0),
        )
    };
    let raw = [
        product([beta, &one, h, gj, ghj], Atom::num(-1)),
        product([&zero, &(-w), gh, j, ghj], if_nonzero(gh, j / gh)),
        product([&zero, &(-beta), g, j, gj], if_nonzero(g, j / g)),
    ];
    // With `products` 1, in binary64, the products of a parallel chart only
    // where [products_kept] keeps them; with 2 all of them, for the chart
    // tests of [tfun_chart] (the inspection form, which is the complete
    // expansion, takes them all: the variants are algebraically identical).
    let binary64 = rounds_like_binary64(&grouped(&(1 + &scale)))
        * if_nonzero(&(products - Atom::num(2)), Atom::num(1));
    let kept = |index: usize| {
        let arguments = [
            &Atom::num(index as i64), &y1.value, &y2.value, beta, w, g, h, j, gh, gj, ghj, &raw[0],
            &raw[1], &raw[2], &a1, &scale,
        ];
        Q::rounding_choice(
            if_nonzero_else(
                &(&binary64 * &parallel),
                if_nonzero_else(&call("products_kept", arguments), raw[index].clone(), Atom::num(0)),
                raw[index].clone(),
            ),
            raw[index].clone(),
        )
    };
    let [product_x1, product_diagonal, product_y0] = [kept(0), kept(1), kept(2)];
    // 1 - beta is -w.
    let integral = s3_line_call(y1, y2, beta, &one, h, gj, dj, ghj, &product_x1)
        - s3_line_call(y1, y2, &zero, &(-w), gh, j, dj, ghj, &product_diagonal)
        + s3_line_call(y1, y2, &zero, &(-beta), g, j, dj, gj, &product_y0)
        - s3_call(y1, y2, beta, &one, b, &(c + e), &(a + dpf), de, ddpf, edge(lambda_bc, b, q10, q11))
        + s3_call(y1, y2, &zero, &(-w), &sum, dpe, f, ddpe, df, edge(lambda_abc, abc, f, q11))
        - s3_call(y1, y2, &zero, &(-beta), a, d, f, dd, df, edge(lambda_a, a, f, q10));
    let real = real_intersections_absent([a, b, c, g, h, d, e, f, j]);
    let residue = call(
        "planar",
        [
            &y1.value, &y1.first, &y1.second, &y2.value, &y2.first, &y2.second, &p1, &p2, &data.a,
            &data.b, &data.c, &data.dc, &real,
        ],
    );
    (integral + if_nonzero_else(&less(&Atom::num(0), &im(beta)), -&residue, residue)) / data.k
}

/// The chart of [tfun] for the contour `beta`: [tfun_lines] with the line
/// products, unless, in binary64, a product offset differs from the rounded
/// `y - z` it replaces by more than 8 units of that subtraction's rounding
/// `2^-52 |z|` (the pole itself is then that inaccurate). There the pole is
/// moved by the shift `e (offset - (y - z))` that the product implies, and
/// the variant whose value moves less is taken: the products, or `y - z`
/// everywhere (as 09b2c74 took). Where the rounded pole carries the error,
/// the line and conic terms of an edge share it and cancel it with `y - z`,
/// which the products' exact offsets in the line terms alone break: with
/// the products alone (f2adce4), eight boxes of the verifier's 32,600-row
/// sweep lost 1e-12 to 1.3e-5 on some binary64 route where 09b2c74 kept
/// 1.8e-14 to 1.3e-6 (v1 7832, 8793, 9431, 13198, 13274, 13881, 14824 and
/// 15393; tests/data/contour_pole_errors.txt). The probe restores 09b2c74's
/// level on all of them but 7832 on the hand-written and generated routes
/// (a reflected chart, [zero_external]) and 8793 on the generated route.
/// Where the subtraction alone loses the digits, `y - z` moves with the
/// pole, and the products are kept (the near-threshold boxes of
/// tests/data/near_threshold_boxes.txt). The hand-written kernel records
/// the same disagreement in `r1_with`.
fn tfun_selected(x: [&Atom; 29]) -> Atom {
    let [g, h, beta, w] = [x[3], x[4], x[11], x[12]];
    let zero = Atom::num(0);
    let lines = |products: i64, [shift1, shift2]: [&Atom; 2]| {
        let products = Atom::num(products);
        let mut arguments: Vec<&Atom> = x.to_vec();
        arguments.extend([&products, shift1, shift2]);
        call::<32>("tfun_lines", arguments.try_into().unwrap())
    };
    let value = lines(1, [&zero, &zero]);
    // [tfun_edge_probe] of the line edge `index`: the disagreement in units
    // (`which` 0) or the shift of pole 1 or 2 (`which` 1, 2).
    let edge = |index: i64, which: i64| {
        let [index, which] = [Atom::num(index), Atom::num(which)];
        let mut arguments: Vec<&Atom> = x.to_vec();
        arguments.extend([&index, &which]);
        call::<31>("tfun_edge_probe", arguments.try_into().unwrap())
    };
    // The largest disagreement; the first edge where several are equal.
    let units = [0, 1, 2].map(|index| edge(index, 0));
    let second = less(&units[0], &units[1]);
    let best = if_nonzero_else(&second, units[1].clone(), units[0].clone());
    let third = less(&best, &units[2]);
    let largest = if_nonzero_else(&third, units[2].clone(), best);
    let shift = |which: i64| {
        if_nonzero_else(
            &third,
            edge(2, which),
            if_nonzero_else(&second, edge(1, which), edge(0, which)),
        )
    };
    let [shift1, shift2] = [shift(1), shift(2)];
    let probe = {
        let moved_products = lines(1, [&shift1, &shift2]);
        let direct = lines(0, [&zero, &zero]);
        let moved_direct = lines(0, [&shift1, &shift2]);
        let take_direct = finite(&moved_products)
            * finite(&direct)
            * finite(&moved_direct)
            * le(&norm(&(&moved_direct - &direct)), &norm(&(&moved_products - &value)));
        if_nonzero_else(&take_direct, direct, value.clone())
    };
    // Only a binary64 chart whose line is nearly parallel to the contour
    // direction takes the products ([tfun_lines]) and is probed; the test
    // takes `a1` alone, so that the edges are evaluated only there.
    let [b, c, d, e, f, j, dpe] = [x[1], x[2], x[5], x[6], x[7], x[8], x[9]];
    let [a1, ..] = beta_polynomial(beta, w, b, c, g, h, d, e, f, j, dpe);
    let scale = norm(g) + norm(&(beta * h));
    let parallel = less(&(Atom::num(4096) * norm(&a1)), &scale);
    let binary64 = rounds_like_binary64(&grouped(&(1 + &scale)));
    // The variants are algebraically identical: the inspection form takes the
    // products without the probe (which made the complete expansion of D0
    // exceed its node budget, tests/full_expansion.rs).
    Q::rounding_choice(
        if_nonzero_else(
            &(binary64 * parallel),
            if_nonzero_else(&less(&Atom::num(64), &largest), probe, value.clone()),
            value.clone(),
        ),
        value,
    )
}

/// One line edge of [tfun_selected]'s probe: the `[d, e, b, c, far]` and the
/// ratio of the edge `index` (x = 1, y = x, y = 0, as in [tfun_lines]), its
/// product ([line_product]) and [tfun_edge_units] with the chart's poles.
/// The arguments are [tfun_selected]'s, then `index` and `which`.
fn tfun_edge_probe(x: [&Atom; 31]) -> Atom {
    let [a, b, c, g, h, d, e, f, j, dpe] = [x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]];
    let [beta, w, dd, de, df, dj, ddpe] = [x[11], x[12], x[13], x[14], x[15], x[16], x[17]];
    let [gh, gj, ghj, k] = [x[25], x[26], x[27], x[28]];
    let [index, which] = [x[29], x[30]];
    let data = beta_data(beta, w, b, c, g, h, d, e, f, j, dpe, dd, de, df, dj, ddpe, Some(k));
    let a1 = if_nonzero_else(&data.a, data.a.clone(), g / Atom::num(1u64 << 52));
    let conic_norm = [a, b, c, d, e, f]
        .into_iter()
        .fold(Atom::num(0), |sum, coefficient| sum + norm(coefficient));
    let [zero, one] = [Atom::num(0), Atom::num(1)];
    let pick = |x1: Atom, diagonal: Atom, y0: Atom| {
        if_nonzero_else(index, if_nonzero_else(&(index - 1), y0, diagonal), x1)
    };
    let ed = pick(beta.clone(), zero.clone(), zero);
    let ee = pick(one, -w, -beta);
    let eb = pick(h.clone(), gh.clone(), g.clone());
    let ec = pick(gj.clone(), j.clone(), j.clone());
    let efar = pick(ghj.clone(), ghj.clone(), gj.clone());
    // The three edges' products, of which [products_kept] keeps this edge's.
    let raw = [
        line_product(&data, &a1, &conic_norm, [beta, &Atom::num(1), h, gj, ghj], &Atom::num(-1), [b, e, f]),
        line_product(&data, &a1, &conic_norm, [&Atom::num(0), &(-w), gh, j, ghj], &if_nonzero(gh, j / gh), [b, e, f]),
        line_product(&data, &a1, &conic_norm, [&Atom::num(0), &(-beta), g, j, gj], &if_nonzero(g, j / g), [b, e, f]),
    ];
    let scale = norm(g) + norm(&(beta * h));
    let kept = call(
        "products_kept",
        [
            index, &data.y[0].value, &data.y[1].value, beta, w, g, h, j, gh, gj, ghj, &raw[0], &raw[1],
            &raw[2], &a1, &scale,
        ],
    );
    let product = if_nonzero(&kept, pick(raw[0].clone(), raw[1].clone(), raw[2].clone()));
    let forced = if_nonzero_else(&(&kept - 2), Atom::num(0), Atom::num(1));
    call(
        "tfun_edge_units",
        [&data.y[0].value, &data.y[1].value, &ed, &ee, &eb, &ec, &efar, &product, &forced, which],
    )
}

/// [r1_line]'s offsets of a line edge `[d, e, b, c, far]` with its product,
/// for the poles `y1`, `y2`: where the product is taken, its disagreement with
/// the rounded `y - z` in units of `2^-52 |z|`, squared (`which` 0; zero where
/// no product is taken; at least 128 where `forced`, a product that
/// [products_kept] probes), or the shift `e (offset - (y - z))` of pole 1 or 2
/// (`which` 1 or 2; zero for the other pole).
fn tfun_edge_units([y1, y2, d, e, b, c, far, product, forced, which]: [&Atom; 10]) -> Atom {
    let one = Atom::num(1);
    let rounding = Atom::num(symbolica::domains::integer::Integer::from(2).pow(104));
    let [p1, p2] = [y1, y2].map(|y| (d + y) / e);
    let z = if_nonzero_else(far, -c / b, Atom::num(1));
    let [o1, o2] = [&p1 - &z, &p2 - &z];
    let [n1, n2] = [norm(&o1), norm(&o2)];
    let first = less(&n1, &n2);
    let near = if_nonzero_else(&first, o1.clone(), o2.clone());
    let far_offset = if_nonzero_else(&first, o2, o1);
    let nearer = if_nonzero_else(&first, n1, n2);
    let hzz = norm(&z);
    let taken = less(&(Atom::num(4096) * &nearer), &hzz)
        * if_nonzero(product, one.clone())
        * if_nonzero(&(&p1 - &p2), one.clone())
        * if_nonzero(b, one.clone())
        * if_nonzero(e, one);
    let offset = if_nonzero(&far_offset, product / &far_offset);
    let offset = if_nonzero_else(&im(&near), offset.clone(), re(&offset));
    let difference = &offset - &near;
    let units = &taken
        * if_nonzero(
            &offset,
            rounding * norm(&difference) * &nearer / (&hzz * norm(&offset)),
        );
    let floor = Atom::num(128);
    let units = if_nonzero_else(
        &(forced * taken),
        if_nonzero_else(&less(&units, &floor), floor.clone(), units.clone()),
        units,
    );
    let shift = difference * e;
    if_nonzero_else(
        which,
        if_nonzero_else(&(which - 1), (1 - &first) * &shift, &first * shift),
        units,
    )
}

/// [tfun_selected] with the contour `use_first` selects ([tfun]), except
/// where that contour's a1 keeps fewer than half of binary64's precision
/// (`2^56 |a1|² < |g|² + |beta h|²`, 24 bits, it is flat) and the other's does
/// not and the other has a pole within 2^7: there the other contour, unless
/// the two charts agree to 2^-26 (or the other is not finite).
///
/// In binary64 arithmetic the chart used is also tested where it is the flat
/// selected one without a regular chart that agrees with it, and where its K
/// rounds to zero (its far pole is then placed by [chart_k]; 09b2c74 returned
/// 0/0 there). The flat chart is NaN unless the chart with the line's `g`
/// moved to `g (1 + 2^-52)` (and the line's sums by the same exact amount)
/// keeps 10 bits of it, and the chart with a vanishing K (or the other chart
/// taken, whose K is `b a1` of the flat one) unless the chart with K at the
/// rounding level of its terms, `2^-52 (|b a1| + |h b1|)`, keeps 13 bits.
/// With the Gram determinant at the rounding level a flat chart's poles sit
/// at the line edges' roots with lips that the rounding decides, and a
/// vanishing K's far pole is a rounding residue: such charts were wrong by up
/// to 0.2 (and NaN at 09b2c74) where no other chart is regular, and the tests
/// make them NaN, which the caller detects (NinjaRS rescues it in higher
/// precision). `DoubleFloat` and `Float` keep their value: they resolve what
/// binary64 does not, and the moves are binary64 rounding units. The moved
/// values are grouped, so that every route rounds the same moved line
/// (moving a1 alone breaks the line products' identity `h u + j = ratio a1`
/// ([line_product])).
///
/// The thresholds are fixed powers of two, so that the hand-written,
/// generated, interpreted, SymJIT and `DoubleFloat` routes take the same chart
/// on binary64 inputs: a test that follows the working precision (`1 + x`
/// rounding to 1) is never true in `DoubleFloat` (which keeps `1 + x` exactly
/// for any representable `x`) and so took other charts there. `Float`
/// ([double_double] tells it apart) takes the selected chart: the agreement
/// rule remedies binary64's rounding, and with binary64's thresholds, or with
/// ones that follow its working precision, two 128-bit boxes of a 32,600-row
/// sweep took the other chart and lost 4.2e-7 and 1.7e-6 (28 and 22 tracked
/// bits), where 09b2c74 kept 4.5e-9 and 1.1e-16.
fn tfun_chart(x: [&Atom; 31]) -> Atom {
    let [beta1, w1, beta2, w2, use_first] = [x[26], x[27], x[28], x[29], x[30]];
    let [b, c, g, h, d, e, f, j, dpe] = [x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]];
    let [dd, de, df, dj, ddpe] = [x[11], x[12], x[13], x[14], x[15]];
    let data = |beta: &Atom, w: &Atom| {
        beta_data(
            beta, w, b, c, g, h, d, e, f, j, dpe, dd, de, df, dj, ddpe, None,
        )
    };
    let [data1, data2] = [data(beta1, w1), data(beta2, w2)];
    let scale = |beta: &Atom| norm(g) + norm(&(beta * h));
    // The arithmetic: binary64 and `DoubleFloat` test the charts with
    // binary64's fixed thresholds, `Float` takes the selected one.
    let binary64 = rounds_like_binary64(&grouped(&(1 + scale(beta1))));
    let float = (1 - &binary64) * (1 - double_double(use_first));
    let flat = |data: &BetaData, beta: &Atom| {
        if_nonzero_else(&float, Atom::num(0), below(56, &norm(&data.a), &scale(beta)))
    };
    let near = |data: &BetaData| {
        let within = |y: &Continued| less(&norm(&y.value), &Atom::num(1 << 14));
        1 - (1 - within(&data.y[0])) * (1 - within(&data.y[1]))
    };
    let [flat1, flat2] = [flat(&data1, beta1), flat(&data2, beta2)];
    let check = if_nonzero_else(
        use_first,
        &flat1 * (1 - &flat2) * near(&data2),
        &flat2 * (1 - &flat1) * near(&data1),
    );
    // The chart of the contour `first`: as selected (`Chart::Plain`); with
    // the line's g moved to `g (1 + 2^-52)` and its sums by the same exact
    // `delta`, as grouped values, so that every route rounds the same moved
    // line (`Chart::MovedLine`); or with K at the rounding level of its
    // terms, `2^-52 (|b a1| + |h b1|)` (`Chart::RoundingK`).
    //
    // These charts, which the tests and the agreement compare, take every line
    // product (`products` 2, as before [products_kept]); the chart returned
    // takes the kept ones (the probed selected chart, or the other chart as
    // `Chart::Kept`).
    #[derive(Clone, Copy, PartialEq)]
    enum Chart {
        Plain,
        MovedLine,
        RoundingK,
        Probed,
        Kept,
    }
    let chart = |first: &Atom, kind: Chart| {
        let beta = if_nonzero_else(first, beta1.clone(), beta2.clone());
        let w = if_nonzero_else(first, w1.clone(), w2.clone());
        let [g, gh, gj, ghj] = if kind == Chart::MovedLine {
            let shifted = grouped(&(g * Atom::num((1u64 << 52) + 1) / Atom::num(1u64 << 52)));
            let delta = grouped(&(&shifted - g));
            let sum = |value: &Atom| grouped(&(value + &delta));
            [shifted, sum(x[23]), sum(x[24]), sum(x[25])]
        } else {
            [g, x[23], x[24], x[25]].map(Clone::clone)
        };
        let k = if kind == Chart::RoundingK {
            let [a1, b1, ..] = beta_polynomial(&beta, &w, b, c, &g, h, d, e, f, j, dpe);
            let magnitude = |z: &Atom| Symbol::ABS.call((re(z),)) + Symbol::ABS.call((im(z),));
            (magnitude(&(b * &a1)) + magnitude(&(h * &b1))) / Atom::num(1u64 << 52)
        } else {
            chart_k(&beta, &w, b, c, &g, h, d, e, f, j, dpe)
        };
        let arguments = [
            x[0], b, c, &g, h, d, e, f, j, dpe, x[10], &beta, &w, dd, de, df, dj, ddpe, x[16],
            x[17], x[18], x[19], x[20], x[21], x[22], &gh, &gj, &ghj, &k,
        ];
        if kind == Chart::Probed {
            return call("tfun_selected", arguments);
        }
        let products = Atom::num(if kind == Chart::Kept { 1 } else { 2 });
        let zero = Atom::num(0);
        let mut arguments: Vec<&Atom> = arguments.to_vec();
        arguments.extend([&products, &zero, &zero]);
        call::<32>("tfun_lines", arguments.try_into().unwrap())
    };
    let value = chart(use_first, Chart::Plain);
    // The selected chart's value where the tests keep it ([tfun_selected]).
    let probed = chart(use_first, Chart::Probed);
    // The chart used (`first`, its value `used`), tested in binary64: where
    // `flat`, NaN unless the chart of the moved line keeps 10 bits of it;
    // where `zero` (its K rounds to zero), NaN unless the chart with K at its
    // rounding level keeps 13 bits.
    let tested = |first: &Atom, used: &Atom, result: &Atom, flat: &Atom, zero: &Atom| {
        let keeps = |kind: Chart, bits: u32| {
            let moved = chart(first, kind);
            finite(&moved) * close(bits, used, &moved)
        };
        let line = if_nonzero_else(flat, keeps(Chart::MovedLine, 10), Atom::num(1));
        let k = if_nonzero_else(zero, keeps(Chart::RoundingK, 13), Atom::num(1));
        if_nonzero_else(&(line * k), result.clone(), undetermined(used))
    };
    let zero_k = |data: &BetaData| &binary64 * if_nonzero_else(&data.k, Atom::num(0), Atom::num(1));
    let [zero1, zero2] = [zero_k(&data1), zero_k(&data2)];
    let either = |p: &Atom, q: &Atom| 1 - (1 - p) * (1 - q);
    let none = Atom::num(0);
    // check: the other chart, taken unless the two agree (or it is not
    // finite). The other chart used is tested as one whose K rounds to zero
    // (its K is `b a1` of the flat selected chart, at most 24 bits), a
    // selected chart that it does not confirm (a non-finite other) as below.
    let other_first = 1 - use_first;
    let other = chart(&other_first, Chart::Plain);
    let kept = chart(&other_first, Chart::Kept);
    let disagree = finite(&other) * (1 - finite(&value) * agree(&value, &other));
    let zero_selected = if_nonzero_else(use_first, zero1.clone(), zero2.clone());
    // Where the other chart's K is a rounding residue too (binary64; `K = b
    // a1 - h b1` below 2^-52 of `|b a1|, |h b1|`), neither chart is regular:
    // of two that disagree, the selected one is kept where K at its rounding
    // level moves the other (relatively) more than 2^10 times as much as the
    // selected moves with its line, and it must then keep 20 bits under that
    // move (13 under K at its rounding level where its own K rounds to zero),
    // or the chart is NaN. Taking the other there lost fs 16695 of the second
    // verifier's 36,840-row sweep (exact vanishing Gram determinant; the
    // other's K rounds to zero and moves it by 2e-6, the selected moves by
    // 1.3e-9) to 5.3e-6 in all 24 relabellings, where the selected chart
    // keeps 2.5e-8 to 1.0e-7 (as 09b2c74), and made fs 15250 NaN (6.3e-9 at
    // 09b2c74 and now). Without the margin, fs 14684 (selected moved by 8e-15
    // with its line, other by 1.7e-14 with its K) took a selected chart wrong
    // by 7e-3; with any `K` below 2^-26 of its terms (not only a residue),
    // SymJIT took a selected chart wrong by 0.13 in v1 11845 (its line
    // moves it by 1e-15 at most on the other routes). The boxes of
    // tests/data/contour_gram_rounding.txt move their other chart by 1e-12
    // or less and the selected by 1e-3 or more, and keep the other.
    let residue = Atom::num(symbolica::domains::integer::Integer::from(2).pow(104));
    let weak_k = |data: &BetaData| {
        &binary64
            * less(
                &(&residue * norm(&data.k)),
                &(norm(&(b * &data.a)) + norm(&(h * &data.b))),
            )
    };
    let weak_other = if_nonzero_else(use_first, weak_k(&data2), weak_k(&data1));
    let moved_selected = chart(use_first, Chart::MovedLine);
    let moved_other = chart(&other_first, Chart::RoundingK);
    let moves_more = less(
        &(norm(&(&moved_other - &other)) * norm(&value)),
        &(Atom::num(1u64 << 20) * norm(&(&moved_selected - &value)) * norm(&other)),
    );
    let prefer_selected = &disagree
        * weak_other
        * finite(&value)
        * finite(&moved_selected)
        * (1 - finite(&moved_other) * moves_more);
    let take_other = &disagree * (1 - &prefer_selected);
    let selected_kept = {
        let line = finite(&moved_selected) * close(20, &value, &moved_selected);
        let rounding_k = chart(use_first, Chart::RoundingK);
        let k = if_nonzero_else(
            &zero_selected,
            finite(&rounding_k) * close(13, &value, &rounding_k),
            Atom::num(1),
        );
        if_nonzero_else(&(line * k), probed.clone(), undetermined(&value))
    };
    let checked = if_nonzero_else(
        &take_other,
        if_nonzero_else(
            &binary64,
            tested(&other_first, &other, &kept, &none, &binary64),
            kept.clone(),
        ),
        if_nonzero_else(
            &prefer_selected,
            selected_kept,
            if_nonzero_else(
                &finite(&other),
                probed.clone(),
                tested(use_first, &value, &probed, &binary64, &zero_selected),
            ),
        ),
    );
    // Otherwise a flat selected chart (no regular alternative) and one whose
    // K rounds to zero are tested.
    let flat_selected = &binary64 * if_nonzero_else(use_first, flat1.clone(), flat2.clone());
    if_nonzero_else(
        &check,
        checked,
        if_nonzero_else(
            &either(&flat_selected, &zero_selected),
            tested(use_first, &value, &probed, &flat_selected, &zero_selected),
            probed.clone(),
        ),
    )
}

/// The two contours of [tfun]: the roots `beta1` (of smaller magnitude) and
/// `beta2` of `b beta² + c beta + a` and their `w = beta - 1`.
fn tfun_betas(a: &Atom, b: &Atom, c: &Atom, abc: &Atom) -> [Atom; 4] {
    let [root1, root2] = roots(b, c, a);
    // beta - 1: the root nearer 1 from (beta1 - 1)(beta2 - 1) = (a + b + c)/b.
    let [o1, o2] = [&root1 - 1, &root2 - 1];
    let w1 = if_nonzero_else(&less(&norm(&o1), &norm(&o2)), abc / (b * &o2), o1.clone());
    let w2 = if_nonzero_else(&less(&norm(&o2), &norm(&o1)), abc / (b * &o1), o2);
    let swap = less(&norm(&root2), &norm(&root1));
    let beta1 = if_nonzero_else(&swap, root2.clone(), root1.clone());
    let beta2 = if_nonzero_else(&swap, root1, root2);
    let (w1, w2) = (
        if_nonzero_else(&swap, w2.clone(), w1.clone()),
        if_nonzero_else(&swap, w1, w2),
    );
    [beta1, w1, beta2, w2]
}

/// Nonzero where [tfun] takes its first contour, from tfun's sign-adjusted
/// `[a, b, c, g, h, d, e, f, j, dpe, dd, de, df, dj, ddpe, abc]`: the contour
/// opposite the pole nearest a boundary, or the first where it encloses no
/// pole. A registered function, so that [tfun] and [tfun_vertex_suspicious]
/// share one copy of the poles it compares.
fn tfun_contour_choice(
    [a, b, c, g, h, d, e, f, j, dpe, dd, de, df, dj, ddpe, abc]: [&Atom; 16],
) -> Atom {
    let [beta1, w1, beta2, w2] = tfun_betas(a, b, c, abc);
    let data = |beta: &Atom, w: &Atom| {
        beta_data(
            beta, w, b, c, g, h, d, e, f, j, dpe, dd, de, df, dj, ddpe, None,
        )
    };
    let data1 = data(&beta1, &w1);
    let data2 = data(&beta2, &w2);
    let distance_y =
        |beta: &Atom, y: &Continued| call("distance", [beta, &y.value, &y.first, &y.second]);
    let d1 = min(
        &distance_y(&beta1, &data1.y[0]),
        &distance_y(&beta1, &data1.y[1]),
    );
    let d2 = min(
        &distance_y(&beta2, &data2.y[0]),
        &distance_y(&beta2, &data2.y[1]),
    );
    if_nonzero_else(&(&d1 - 1), less(&d2, &d1), Atom::num(1))
}

/// Nonzero where a pole `y` of the contour [tfun] selects lies next to the
/// vertex `-beta` of its triangle `(0, -beta, 1 - beta)`, the image of the
/// corner (1, 0), with an offset `y + beta` that the rounding resolves only
/// to a factor 16 of itself (`max(|beta|, |y|) > 16 |y + beta|`), or, with the
/// pole's own rounding, to 2^12 of itself ([tfun_vertex_test]). Zero for the
/// `t1` branches (`a` or `b` zero). [zero_external] reflects the chart only
/// there.
fn tfun_vertex_suspicious(
    [
        a,
        b,
        c,
        gin,
        hin,
        d,
        e,
        f,
        jin,
        dpe,
        _dpf,
        dd,
        de,
        df,
        djin,
        ddpe,
        _ddpf,
        abc,
    ]: [&Atom; 18],
) -> Atom {
    let sj = imaginary_sign(jin, &imaginary_sign(djin, &Atom::num(-1)));
    let g = -&sj * gin;
    let h = -&sj * hin;
    let j = -&sj * jin;
    let dj = -&sj * djin;
    let [beta1, w1, beta2, w2] = tfun_betas(a, b, c, abc);
    let use_first = call(
        "tfun_contour_choice",
        [
            a, b, c, &g, &h, d, e, f, &j, dpe, dd, de, df, &dj, ddpe, abc,
        ],
    );
    let test = |beta: &Atom, w: &Atom| {
        call(
            "tfun_vertex_test",
            [beta, w, b, c, &g, &h, d, e, f, &j, dpe],
        )
    };
    let selected = if_nonzero_else(&use_first, test(&beta1, &w1), test(&beta2, &w2));
    if_nonzero_else(b, if_nonzero_else(a, selected, Atom::num(0)), Atom::num(0))
}

/// [tfun_vertex_suspicious]'s test of the contour `(beta, w)`, from
/// `[beta, w, b, c, g, h, d, e, f, j, dpe]`: a pole `y` of `K y² + L y + N`
/// ([beta_data]) with `max(|beta|, |y|) > 16 |y + beta|`, or with
/// `|beta| + dy > 4096 |y + beta|`, where `dy = (Kt |y|² + Lt |y| + Nt)/|2 K y
/// + L|` is the root's sensitivity to the rounding of the products that form
/// `K = b a1 - h b1`, `L = e a1 - h c1 - j b1` and `N = f a1 - j c1` (`Kt`,
/// `Lt`, `Nt` their magnitudes).
fn tfun_vertex_test([beta, w, b, c, g, h, d, e, f, j, dpe]: [&Atom; 11]) -> Atom {
    let [a1, b1, c1, l, n] = beta_polynomial(beta, w, b, c, g, h, d, e, f, j, dpe);
    let k = b * &a1 - h * &b1;
    let size = |z: &Atom| Symbol::ABS.call((z,));
    let either = |p: &Atom, q: &Atom| 1 - (1 - p) * (1 - q);
    let kt = size(b) * size(&a1) + size(h) * size(&b1);
    let lt = size(e) * size(&a1) + size(h) * size(&c1) + size(j) * size(&b1);
    let nt = size(f) * size(&a1) + size(j) * size(&c1);
    let mb = size(beta);
    let pole = |y: &Atom| {
        let offset = size(&(y + beta));
        let my = size(y);
        let larger = if_nonzero_else(&less(&mb, &my), my.clone(), mb.clone());
        let geometric = less(&(Atom::num(16) * &offset), &larger);
        let dy = (&kt * &my * &my + &lt * &my + &nt) / size(&(Atom::num(2) * &k * y + &l));
        either(&geometric, &less(&(Atom::num(4096) * &offset), &(&mb + dy)))
    };
    let [y1, y2] = roots(&k, &l, &n);
    either(&pole(&y1), &pole(&y2))
}

fn tfun(
    [
        a,
        b,
        c,
        gin,
        hin,
        d,
        e,
        f,
        jin,
        dpe,
        dpf,
        dd,
        de,
        df,
        djin,
        ddpe,
        ddpf,
        abc,
        lambda_a,
        lambda_abc,
        lambda_bc,
        q10,
        q11,
        ghin,
        gjin,
        ghjin,
    ]: [&Atom; 26],
) -> Atom {
    let sj = imaginary_sign(jin, &imaginary_sign(djin, &Atom::num(-1)));
    let g = -&sj * gin;
    let h = -&sj * hin;
    let j = -&sj * jin;
    let dj = -&sj * djin;
    // The line's coefficient sums, formed by the caller from the invariants:
    // a vanishing one is an exactly vanishing linear coefficient or an exact
    // root at an end of the edge it restricts to (`gh` = p3 - p5 - p0 in the
    // zero-external chart: for p3 = p5 the rounded g + h was 1e-16 and put the
    // root of one edge at 1e15).
    let gh = -&sj * ghin;
    let gj = -&sj * gjin;
    let ghj = -&sj * ghjin;
    let [beta1, w1, beta2, w2] = tfun_betas(a, b, c, abc);
    let use_first = call(
        "tfun_contour_choice",
        [
            a, b, c, &g, &h, d, e, f, &j, dpe, dd, de, df, &dj, ddpe, abc,
        ],
    );
    // The contour opposite the pole nearest a boundary is used. If the first
    // contour encloses no pole, the first (smaller-magnitude) beta is
    // selected.
    //
    // Next to a vanishing Gram determinant one contour's a1 = g + beta h
    // vanishes (the line is nearly parallel to its direction (1, beta)) and
    // the other's K = b (g + beta' h) with it. With the line products
    // ([line_product]) the first chart keeps its digits for a small a1, but
    // as a1 approaches the rounding level of its terms its poles sit at the
    // line edges' roots with lips decided by the rounding, and it can be
    // wrong by O(1): by 35 times the value in all 24 relabellings of a box
    // with the Gram determinant at the rounding level (row 1 of
    // tests/data/contour_gram_rounding.txt). The other chart has a pole at infinity
    // ([chart_k]) and is regular where its other pole stays within 2^7 of the
    // triangle (with both poles far, K and L vanish together and its terms
    // cancel). Where the selected contour's a1 keeps fewer than 24 bits
    // (`2^56 |a1|² < |g|² + |beta h|²`), the other's does not and the other
    // has a pole within 2^7, [tfun_chart] evaluates both charts and keeps the
    // selected one only if the two agree to 2^-26; in binary64 it returns
    // NaN where the chart used, flat without a regular alternative or with a
    // K that rounds to zero, does not keep 10 bits under a move of its line
    // or of its K within their rounding.
    let general = call(
        "tfun_chart",
        [
            a, b, c, &g, &h, d, e, f, &j, dpe, dpf, dd, de, df, &dj, ddpe, ddpf, abc, lambda_a,
            lambda_abc, lambda_bc, q10, q11, &gh, &gj, &ghj, &beta1, &w1, &beta2, &w2, &use_first,
        ],
    );
    // b = 0: a + c = abc.
    let no_b = call(
        "t1",
        [
            a, c, &g, &h, d, e, f, &j, dpe, dd, de, df, &dj, ddpe, abc, lambda_abc, q11,
            lambda_a, q10, &gh, &gj, &ghj,
        ],
    );
    // a = 0: b + c = abc, and the new a + c is b. The new conic takes the
    // values q11 (its f), q10 and f at the corners (0, 0), (1, 0) and (1, 1),
    // and the new line (-g - h, g, g + h + j) has the sums -h, j and g + j.
    let bc = coefficient_sum(&[b, c], abc);
    let no_a = call(
        "t1",
        [
            &bc,
            &(-c),
            &(-&gh),
            &g,
            &(-dpe - Atom::num(2) * &bc),
            &(d + c),
            q11,
            &ghj,
            &(-e - Atom::num(2) * b - c),
            &(-ddpe),
            dd,
            &(ddpe + df),
            &dj,
            &(-de),
            b,
            lambda_bc,
            q10,
            lambda_abc,
            f,
            &(-&h),
            &j,
            &gj,
        ],
    );
    -sj * if_nonzero_else(b, if_nonzero_else(a, general, no_a), no_b)
}

fn no_opposite_chart(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    let dm = mass_directions(m);
    let a = p[4] - p[0];
    let c = p[3] - p[4] - p[2];
    let g = p[1];
    let h = p[5] - p[1] - p[2];
    let d = in_order(&[m[2] - m[3], -p[2]]);
    let e = in_order(&[m[0] - m[2], p[2].clone(), -p[3]]);
    let f = m[3];
    let k = in_order(&[m[1] - m[2], -p[5], p[2].clone()]);
    let dpe = in_order(&[m[0] - m[3], -p[3]]);
    let dpk = in_order(&[m[1] - m[3], -p[5]]);
    let dpf = m[2] - p[2];
    // The edge quadratics' values at x = 1 (m0 and m1 for p0 = p2 = 0).
    let far_ae = in_order(&[m[0].clone(), -p[0], -p[2]]);
    let far_gk = m[1] - p[2];
    call(
        "t13",
        [
            &a,
            &c,
            g,
            &h,
            &d,
            &e,
            f,
            &k,
            &dpe,
            &dpk,
            &dpf,
            &(&dm[2] - &dm[3]),
            &(&dm[0] - &dm[2]),
            &dm[3],
            &(&dm[1] - &dm[2]),
            &(&dm[0] - &dm[3]),
            &(&dm[1] - &dm[3]),
            &dm[2],
            &kallen_call(&a, &dpf, &far_ae),
            &kallen_call(g, &dpf, &far_gk),
            &kallen_call(&(g + &h), f, &far_gk),
            &kallen_call(&(&a + &c), f, &far_ae),
            &far_ae,
            &far_gk,
            &(g + &h),
            &(&a + &c),
        ],
    )
}

fn no_opposite_chart_regular(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    let a = p[4] - p[0];
    let c = p[3] - p[4] - p[2];
    let g = p[1];
    let h = p[5] - p[1] - p[2];
    let d = in_order(&[m[2] - m[3], -p[2]]);
    let e = in_order(&[m[0] - m[2], p[2].clone(), -p[3]]);
    let j = in_order(&[m[1] - m[2], -p[5], p[2].clone()]);
    let k = &h * &a - &c * g;
    let l = &a * &d + &h * &e - &d * g - &c * &j;
    let n = &d * (&e - &j) + (&h - &c) * m[3];
    let endpoints = &n * (&k + &l + &n);
    if_nonzero_else(
        &any_nonzero([p[0].clone(), p[2].clone()]),
        Atom::num(0),
        if_nonzero(&endpoints, Atom::num(1)),
    )
}

fn no_opposite_externals(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    // The separate S3 terms diverge when an auxiliary pole is exactly an
    // integration endpoint, although their t13 sum can be regular. A cyclic
    // change of Feynman parameters is exact and often removes this chart
    // singularity. Only use charts preserving the two opposite zero external
    // invariants, and choose by exact polynomial endpoint tests, not a cutoff.
    const CYCLES: [[usize; 6]; 4] = [
        [0, 1, 2, 3, 4, 5],
        [1, 2, 3, 0, 5, 4],
        [2, 3, 0, 1, 4, 5],
        [3, 0, 1, 2, 5, 4],
    ];
    let arguments = CYCLES.map(|cycle| {
        core::array::from_fn::<_, 10, _>(|idx| {
            if idx < 6 {
                p[cycle[idx]].clone()
            } else {
                m[cycle[idx - 6]].clone()
            }
        })
    });
    let mut selected = Atom::num(1);
    for rotation in (0..4).rev() {
        let regular = call("opposite_regular", arguments[rotation].each_ref());
        selected = if_nonzero_else(&regular, Atom::num((rotation + 1) as i64), selected);
    }
    let selected: [Atom; 10] =
        core::array::from_fn(|idx| pick(&selected, arguments.iter().map(|row| row[idx].clone())));
    call("opposite_chart", selected.each_ref())
}

/// [zero_external_chart] in the orientation that puts a massless line's
/// corner of the Feynman-parameter triangle at its origin when the line
/// passes next to it.
///
/// The reflection `(p1, p4, p3, p2, s23, s12; m1, m0, m3, m2)` keeps the
/// vanishing `p1²` first and exchanges the corners (0, 0) and (1, 0): the
/// conic's values there (`f = m3`, `q10 = m2`) and the line's (`k` and
/// `h + k`). With a massless `m2` and the line nearer (1, 0) than (0, 0), a
/// pole lies next to that corner, at an `x` next to 1 whose offset binary64
/// resolves only to eps; reflected, it lies next to x = 0 with its relative
/// precision. The box `4 1.0 A 2.25 -0.3298 0 A A' (A - 1.3458e-4 i) 0 2.25
/// (A - 1.3458e-4 i)`, `A = 0.0154585...` and `A'` one ulp below, was wrong by
/// 44 to 128 times its value (and NaN with the line products), and with `A'`
/// 1e-8 or 1e-4 below `A` by 0.4 and 6e-6; reflected, it keeps 2e-12.
///
/// The line must pass nearer (1, 0) by more than 1/128 (`64 |at_10|² < 63
/// |at_00|²`): where the corners are equally near, the reflection moves no
/// pole away from a corner, and of the verifier's 32,600 boxes two took a
/// reflected chart that lost 1.2e-12 (v1 777, both corners equally near
/// within 0.6 %) and 0.3 to 1.2 times the value (v1 15934, equally near
/// within rounding), where the original orientation keeps 2.6e-14 and
/// 5.4e-14. The nearest box that the reflection improves has the line 2 %
/// nearer (1, 0) (v1 12400: 1.4e-6 unreflected, 1.2e-13 reflected).
///
/// The line's distance alone does not tell whether a pole lies next to the
/// corner: the reflection is taken only where one of the two [tfun] selects
/// a contour with a pole that the rounding resolves poorly next to the
/// corner's vertex `-beta` ([zero_external_suspicious],
/// [tfun_vertex_suspicious]). Elsewhere the reflected chart has no such pole
/// to gain and can lose to its own cancellations: a far pole of a nearly
/// vanishing `b` (v1 1804 of the verifier's sweep), a line product (v1 7832)
/// or a `K` that cancels by 2e3 and 400 (v1 16613, g1 10856). Reflected, these
/// four lost 1.3e-12 to 1.0e-11 on their worst binary64 route; the original
/// orientation keeps at most 8.2e-13. Of the 1,505 boxes where the line test
/// holds (the verifier's v1, CT, MS, NR and g1 sets, the NinjaRS streams'
/// sweep, the tests' rows in their 24 relabellings and the 22 benchmark
/// streams), 739 keep the original orientation, and none that the
/// reflection improves by more than 10x (above 1e-12) on any binary64 route;
/// the nearest have a pole 20 times nearer the vertex than its size (SW 29392
/// of the streams' sweep, 4.1e-12 unreflected) or one whose own rounding the
/// offset from the vertex magnifies 1e4 times (v1 17278, 1.8e-11 unreflected).
fn zero_external(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    let at_10 = in_order(&[m[0] - m[1], p[1].clone(), -p[4]]);
    let at_00 = in_order(&[m[0] - m[1], p[5].clone(), -p[3]]);
    // Not into a chart whose conic has a double contour direction: the
    // discriminant `c² - 4 a b` of `b beta² + c beta + a` (a Källén function
    // of three of the invariants) within 2^-26 of its terms (binary64's half
    // precision, on every route; [below]). There the two
    // contours nearly coincide and the chart is singular: reflected, a box
    // with the Gram determinant at the rounding level kept 2.4e-6 where it
    // keeps 2.3e-9 unreflected (row 2 of tests/data/contour_gram_rounding.txt),
    // and another 1.9e-7 where it keeps 2.7e-9.
    let distinct = |a: &Atom, b: &Atom, c: Atom| {
        let [cc, ab] = [&c * &c, Atom::num(4) * a * b];
        1 - below(52, &norm(&(&cc - &ab)), &(norm(&cc) + norm(&ab)))
    };
    let reflect = if_nonzero_else(
        m[2],
        Atom::num(0),
        if_nonzero_else(
            m[3],
            if_nonzero_else(
                p[1],
                less(
                    &(Atom::num(64) * norm(&at_10)),
                    &(Atom::num(63) * norm(&at_00)),
                )
                    * distinct(p[2], p[3], p[4] - p[3] - p[2])
                    * distinct(p[2], &(p[5] - p[0]), p[1] - p[5] - p[2]),
                Atom::num(0),
            ),
            Atom::num(0),
        ),
    );
    let reflected = [p[0], p[3], p[2], p[1], p[5], p[4], m[1], m[0], m[3], m[2]];
    let original = [p[0], p[1], p[2], p[3], p[4], p[5], m[0], m[1], m[2], m[3]];
    // Only where the chart itself has a pole next to the corner (1, 0) that
    // the rounding resolves poorly ([zero_external_suspicious]). The two
    // orientations are algebraically identical: the inspection form keeps
    // the test above alone.
    let suspicious = call("zero_external_suspicious", original);
    let reflect = Q::rounding_choice(if_nonzero_else(&reflect, suspicious, Atom::num(0)), reflect);
    let selected: [Atom; 10] = core::array::from_fn(|index| {
        if_nonzero_else(&reflect, reflected[index].clone(), original[index].clone())
    });
    call("zero_external_chart", selected.each_ref())
}

fn zero_external_chart(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    let [first, second] = zero_external_tfun_arguments(p, m);
    call("tfun", first.each_ref()) - call("tfun", second.each_ref())
}

/// Nonzero where one of the two [tfun] of [zero_external_chart] selects a
/// contour with a pole next to the corner (1, 0) that the rounding resolves
/// poorly ([tfun_vertex_suspicious]).
fn zero_external_suspicious(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    let [first, second] = zero_external_tfun_arguments(p, m);
    let flag =
        |x: &[Atom; 26]| call::<18>("tfun_vertex_suspicious", core::array::from_fn(|i| &x[i]));
    1 - (1 - flag(&first)) * (1 - flag(&second))
}

/// The arguments of the two [tfun] of [zero_external_chart].
fn zero_external_tfun_arguments(p: [&Atom; 6], m: [&Atom; 4]) -> [[Atom; 26]; 2] {
    let dm = mass_directions(m);
    let a = p[2];
    let b = p[1];
    let c = p[5] - p[1] - p[2];
    let h = p[3] - p[4] - p[5] + p[1];
    let j = p[4] - p[0] - p[1];
    let d = in_order(&[m[2] - m[3], -p[2]]);
    let e = in_order(&[m[1] - m[2], -p[5], p[2].clone()]);
    let k = in_order(&[m[0] - m[1], p[5].clone(), -p[3]]);
    let f = m[3];
    let cph = p[3] - p[4] - p[2];
    let dpe = in_order(&[m[1] - m[3], -p[5]]);
    let epk = in_order(&[m[0] - m[2], p[2].clone(), -p[3]]);
    let dek = in_order(&[m[0] - m[3], -p[3]]);
    let dpf = m[2] - p[2];
    let bj = b + &j;
    let abc = a + b + &c;
    let [lambda_a, lambda_abc, lambda_bc] = edge_kallens(a, b, &abc, f, m[2], m[1]);
    let abj = a + &bj + &cph;
    let corner = m[0] - p[0];
    let [_, lambda_abj, lambda_bj] = edge_kallens(a, &bj, &abj, f, m[2], &corner);
    // The line (h, j, k): h + j, h + k and h + j + k.
    let line = [
        in_order(&[p[3].clone(), -p[5], -p[0]]),
        in_order(&[m[0] - m[1], p[1].clone(), -p[4]]),
        in_order(&[m[0] - m[1], -p[0]]),
    ];
    [
        [
            a,
            b,
            &c,
            &h,
            &j,
            &d,
            &e,
            f,
            &k,
            &dpe,
            &dpf,
            &(&dm[2] - &dm[3]),
            &(&dm[1] - &dm[2]),
            &dm[3],
            &(&dm[0] - &dm[1]),
            &(&dm[1] - &dm[3]),
            &dm[2],
            &abc,
            &lambda_a,
            &lambda_abc,
            &lambda_bc,
            m[2],
            m[1],
            &line[0],
            &line[1],
            &line[2],
        ],
        [
            a,
            &bj,
            &cph,
            &h,
            &j,
            &d,
            &epk,
            f,
            &k,
            &dek,
            &dpf,
            &(&dm[2] - &dm[3]),
            &(&dm[0] - &dm[2]),
            &dm[3],
            &(&dm[0] - &dm[1]),
            &(&dm[0] - &dm[3]),
            &dm[2],
            &abj,
            &lambda_a,
            &lambda_abj,
            &lambda_bj,
            m[2],
            &corner,
            &line[0],
            &line[1],
            &line[2],
        ],
    ]
    .map(|arguments| arguments.map(Clone::clone))
}

fn ordered(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    let dm = mass_directions(m);
    let a = p[2];
    let b = p[1];
    let g = p[0];
    let c = p[5] - p[1] - p[2];
    let h = p[3] - p[4] - p[5] + p[1];
    let j = p[4] - p[0] - p[1];
    let d = in_order(&[m[2] - m[3], -p[2]]);
    let e = in_order(&[m[1] - m[2], -p[5], p[2].clone()]);
    let k = in_order(&[m[0] - m[1], p[5].clone(), -p[3]]);
    let f = m[3];
    let abc = p[5];
    let bgj = p[4];
    let jph = p[3] - p[0] - p[5];
    let cph = p[3] - p[4] - p[2];
    let dpe = in_order(&[m[1] - m[3], -p[5]]);
    let epk = in_order(&[m[0] - m[2], p[2].clone(), -p[3]]);
    let dek = in_order(&[m[0] - m[3], -p[3]]);
    let dpf = m[2] - p[2];
    let def = m[1] - p[5];
    let [x1, x2] = roots(g, &j, b);
    let alpha = if_nonzero_else(&less(&norm(&re(&x1)), &norm(&re(&x2))), x2, x1);
    let o1 = 1 - &alpha;
    let j1 = &j + Atom::num(2) * g * &alpha;
    let e1 = &e + &k * &alpha;
    let de1 = &dm[1] - &dm[2] + (&dm[0] - &dm[1]) * &alpha;
    // The six Källén functions λ(p_ij, m_i, m_j) of the box.
    let first = edge_kallens(abc, g, &(abc + g + &jph), f, m[1], m[0]);
    let second = edge_kallens(a, bgj, &(a + bgj + &cph), f, m[2], m[0]);
    let third = edge_kallens(a, b, &(a + b + &c), f, m[2], m[1]);
    // The lines and their coefficient sums (g + h, g + j, g + h + j).
    let sums = |[g, h, j]: &[Atom; 3]| {
        [
            in_order(&[g.clone(), h.clone()]),
            in_order(&[g.clone(), j.clone()]),
            in_order(&[g.clone(), h.clone(), j.clone()]),
        ]
    };
    let line1 = [&c + Atom::num(2) * b + (&h + &j) * &alpha, j1.clone(), e1.clone()];
    let line2 = [&c + &h * &alpha, &o1 * &j1, e1.clone()];
    let line3 = [&c + &h * &alpha, -&j1 * &alpha, e1.clone()];
    let [sums1, sums2, sums3] = [&line1, &line2, &line3].map(sums);
    -call(
        "tfun",
        [
            abc,
            g,
            &jph,
            &line1[0],
            &line1[1],
            &dpe,
            &k,
            f,
            &e1,
            &dek,
            &def,
            &(&dm[1] - &dm[3]),
            &(&dm[0] - &dm[1]),
            &dm[3],
            &de1,
            &(&dm[0] - &dm[3]),
            &dm[1],
            &(abc + g + &jph),
            &first[0],
            &first[1],
            &first[2],
            m[1],
            m[0],
            &sums1[0],
            &sums1[1],
            &sums1[2],
        ],
    ) + &o1
        * call(
            "tfun",
            [
                a,
                bgj,
                &cph,
                &line2[0],
                &line2[1],
                &d,
                &epk,
                f,
                &e1,
                &dek,
                &dpf,
                &(&dm[2] - &dm[3]),
                &(&dm[0] - &dm[2]),
                &dm[3],
                &de1,
                &(&dm[0] - &dm[3]),
                &dm[2],
                &(a + bgj + &cph),
                &second[0],
                &second[1],
                &second[2],
                m[2],
                m[0],
                &sums2[0],
                &sums2[1],
                &sums2[2],
            ],
        )
        + &alpha
            * call(
                "tfun",
                [
                    a,
                    b,
                    &c,
                    &line3[0],
                    &line3[1],
                    &d,
                    &e,
                    f,
                    &e1,
                    &dpe,
                    &dpf,
                    &(&dm[2] - &dm[3]),
                    &(&dm[1] - &dm[2]),
                    &dm[3],
                    &de1,
                    &(&dm[1] - &dm[3]),
                    &dm[2],
                    &(a + b + &c),
                    &third[0],
                    &third[1],
                    &third[2],
                    m[2],
                    m[1],
                    &sums3[0],
                    &sums3[1],
                    &sums3[2],
                ],
            )
}

const ORD: [[usize; 10]; 24] = [
    [1, 2, 3, 4, 1, 2, 3, 4, 5, 6],
    [4, 2, 3, 1, 6, 2, 5, 4, 3, 1],
    [2, 4, 3, 1, 6, 3, 5, 1, 2, 4],
    [1, 4, 3, 2, 4, 3, 2, 1, 5, 6],
    [4, 1, 3, 2, 4, 5, 2, 6, 3, 1],
    [2, 1, 3, 4, 1, 5, 3, 6, 2, 4],
    [2, 3, 4, 1, 2, 3, 4, 1, 6, 5],
    [1, 3, 4, 2, 5, 3, 6, 1, 4, 2],
    [3, 1, 4, 2, 5, 4, 6, 2, 3, 1],
    [2, 1, 4, 3, 1, 4, 3, 2, 6, 5],
    [1, 2, 4, 3, 1, 6, 3, 5, 4, 2],
    [3, 2, 4, 1, 2, 6, 4, 5, 3, 1],
    [3, 4, 1, 2, 3, 4, 1, 2, 5, 6],
    [2, 4, 1, 3, 6, 4, 5, 2, 1, 3],
    [4, 2, 1, 3, 6, 1, 5, 3, 4, 2],
    [3, 2, 1, 4, 2, 1, 4, 3, 5, 6],
    [2, 3, 1, 4, 2, 5, 4, 6, 1, 3],
    [4, 3, 1, 2, 3, 5, 1, 6, 4, 2],
    [4, 1, 2, 3, 4, 1, 2, 3, 6, 5],
    [3, 1, 2, 4, 5, 1, 6, 3, 2, 4],
    [1, 3, 2, 4, 5, 2, 6, 4, 1, 3],
    [4, 3, 2, 1, 3, 2, 1, 4, 6, 5],
    [3, 4, 2, 1, 3, 6, 1, 5, 2, 4],
    [1, 4, 2, 3, 4, 6, 2, 5, 1, 3],
];

const PRM_A: [i64; 24] = [
    17, 7, 22, 18, 8, 21, 23, 13, 4, 24, 14, 3, 5, 19, 10, 6, 20, 9, 11, 1, 16, 12, 2, 15,
];
const PRM_B: [i64; 24] = [
    5, 4, 1, 6, 3, 2, 11, 10, 7, 12, 9, 8, 17, 16, 13, 18, 15, 14, 23, 22, 19, 24, 21, 20,
];

fn sort_four([a, b, c, d]: [&Atom; 4]) -> Atom {
    let outputs = [
        10, 9, 8, 4, 11, 12, 7, 3, 24, 23, 22, 2, 14, 13, 18, 5, 15, 16, 17, 6, 19, 20, 21, 1,
    ];
    let leaf = |index: usize| Atom::num(outputs[index - 1]);
    let choose = |x: &Atom, y: &Atom, yes: Atom, no: Atom| if_nonzero_else(&le(x, y), yes, no);
    let branch1 = choose(
        c,
        d,
        leaf(1),
        choose(b, d, leaf(2), choose(a, d, leaf(3), leaf(4))),
    );
    let branch2 = choose(
        b,
        d,
        leaf(5),
        choose(c, d, leaf(6), choose(a, d, leaf(7), leaf(8))),
    );
    let branch3 = choose(
        b,
        d,
        leaf(9),
        choose(a, d, leaf(10), choose(c, d, leaf(11), leaf(12))),
    );
    let branch4 = choose(
        c,
        d,
        leaf(13),
        choose(a, d, leaf(14), choose(b, d, leaf(15), leaf(16))),
    );
    let branch5 = choose(
        a,
        d,
        leaf(17),
        choose(c, d, leaf(18), choose(b, d, leaf(19), leaf(20))),
    );
    let branch6 = choose(
        a,
        d,
        leaf(21),
        choose(b, d, leaf(22), choose(c, d, leaf(23), leaf(24))),
    );
    choose(
        a,
        b,
        choose(b, c, branch1, choose(a, c, branch2, branch3)),
        choose(a, c, branch4, choose(b, c, branch5, branch6)),
    )
}

fn pick(index: &Atom, values: impl IntoIterator<Item = Atom>) -> Atom {
    let values: Vec<_> = values.into_iter().collect();
    let mut result = values[values.len() - 1].clone();
    for (slot, value) in values[..values.len() - 1].iter().enumerate().rev() {
        result = if_nonzero_else(
            &(index - Atom::num((slot + 1) as i64)),
            result,
            value.clone(),
        );
    }
    result
}

fn lambdas(p: [&Atom; 6]) -> [Atom; 4] {
    [(0, 1, 4), (1, 2, 5), (2, 3, 4), (3, 0, 5)]
        .map(|(a, b, c)| re(&((p[c] - p[a] - p[b]).pow(2) - Atom::num(4) * p[a] * p[b])))
}

fn general(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    let lambda = lambdas(p);
    let sorted = call("sort", lambda.each_ref());
    let all_positive = lambda
        .iter()
        .fold(Atom::num(1), |pred, l| pred * (1 - negative(l)));
    let permutation = if_nonzero_else(
        &all_positive,
        pick(&sorted, PRM_B.map(Atom::num)),
        pick(&sorted, PRM_A.map(Atom::num)),
    );
    let selected_m: [Atom; 4] = core::array::from_fn(|slot| {
        pick(&permutation, ORD.iter().map(|row| m[row[slot] - 1].clone()))
    });
    let selected_p: [Atom; 6] = core::array::from_fn(|slot| {
        pick(
            &permutation,
            ORD.iter().map(|row| p[row[slot + 4] - 1].clone()),
        )
    });
    call(
        "ordered",
        [
            &selected_p[0],
            &selected_p[1],
            &selected_p[2],
            &selected_p[3],
            &selected_p[4],
            &selected_p[5],
            &selected_m[0],
            &selected_m[1],
            &selected_m[2],
            &selected_m[3],
        ],
    )
}

fn body(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    // The same canonical external-zero sectors used by boxc's permtable.
    const PERM: [[usize; 6]; 16] = [
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
    const CASE: [usize; 16] = [0, 1, 1, 2, 1, 5, 2, 3, 1, 2, 5, 3, 2, 3, 3, 4];
    let mask = p[..4]
        .iter()
        .enumerate()
        .fold(Atom::num(0), |s, (idx, value)| {
            s + if_nonzero(value, Atom::num(1 << (3 - idx)))
        });
    let values = PERM.iter().enumerate().map(|(sector, order)| {
        let p = order.map(|idx| p[idx]);
        let m: [&Atom; 4] = core::array::from_fn(|idx| m[order[idx]]);
        let args = [p[0], p[1], p[2], p[3], p[4], p[5], m[0], m[1], m[2], m[3]];
        call(
            if sector == 15 {
                "general"
            } else if matches!(CASE[sector], 0 | 1 | 5) {
                "opposite"
            } else {
                "zero_external"
            },
            args,
        )
    });
    pick(&(mask + 1), values)
}

fn define(map: &mut FunctionMap, name: &str, arity: usize, build: impl FnOnce(&[Atom]) -> Atom) {
    let variables = (0..arity)
        .map(|idx| symbolica::symbol!(format!("__olo_boxc_arg_{idx}")))
        .collect::<Vec<_>>();
    let args = variables
        .iter()
        .map(|symbol| symbol.to_atom())
        .collect::<Vec<_>>();
    map.add_function_with_options(
        function(name),
        variables,
        build(&args),
        native_function_options(),
    )
    .unwrap();
}

/// Registers [norm]. Inlined, so that the evaluators bind each operand once
/// (as the written-out norm did); the definitions keep one copy of it per use.
fn register_norm(map: &mut FunctionMap) {
    let argument = symbolica::symbol!("__olo_boxc_arg_0");
    map.add_function_with_options(
        function("norm"),
        vec![argument],
        norm_body(&argument.to_atom()),
        symbolica::evaluate::FunctionRegistrationOptions::new()
            .inlining(symbolica::evaluate::InliningPolicy::Always),
    )
    .unwrap();
}

pub(crate) fn register(map: &mut FunctionMap) {
    register_norm(map);
    define(map, "root_1", 3, |x| {
        root_bodies(&x[0], &x[1], &x[2])[0].clone()
    });
    define(map, "root_2", 3, |x| {
        root_bodies(&x[0], &x[1], &x[2])[1].clone()
    });
    define(map, "kallen", 3, |x| kallen(core::array::from_fn(|idx| &x[idx])));
    define(map, "r0", 6, |x| r0(core::array::from_fn(|idx| &x[idx])));
    define(map, "r0_line", 7, |x| {
        r0_line(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "r1", 10, |x| r1(core::array::from_fn(|idx| &x[idx])));
    define(map, "r1_line", 11, |x| {
        r1_line(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "s3", 17, |x| s3(core::array::from_fn(|idx| &x[idx])));
    define(map, "s3_line", 13, |x| {
        s3_line(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "s3_double_real", 4, |x| {
        s3_double_real(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "planar", 13, |x| {
        planar(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "inside", 4, |x| {
        inside(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "distance", 4, |x| {
        distance(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "t1", 22, |x| t1(core::array::from_fn(|idx| &x[idx])));
    define(map, "t13", 26, |x| t13(core::array::from_fn(|idx| &x[idx])));
    define(map, "tfun_lines", 32, |x| {
        tfun_lines(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "products_kept", 16, |x| {
        products_kept(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "tfun_edge_units", 10, |x| {
        tfun_edge_units(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "tfun_edge_probe", 31, |x| {
        tfun_edge_probe(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "tfun_selected", 29, |x| {
        tfun_selected(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "tfun_chart", 31, |x| {
        tfun_chart(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "tfun", 26, |x| {
        tfun(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "tfun_contour_choice", 16, |x| {
        tfun_contour_choice(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "tfun_vertex_suspicious", 18, |x| {
        tfun_vertex_suspicious(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "tfun_vertex_test", 11, |x| {
        tfun_vertex_test(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "sort", 4, |x| {
        sort_four(core::array::from_fn(|idx| &x[idx]))
    });
    for (name, build) in [
        (
            "opposite_regular",
            no_opposite_chart_regular as fn([&Atom; 6], [&Atom; 4]) -> Atom,
        ),
        ("opposite_chart", no_opposite_chart),
        ("opposite", no_opposite_externals),
        ("zero_external", zero_external),
        ("zero_external_chart", zero_external_chart),
        ("zero_external_suspicious", zero_external_suspicious),
        ("ordered", ordered),
        ("general", general),
        ("box", body),
    ] {
        define(map, name, 10, |x| {
            build(
                core::array::from_fn(|idx| &x[idx]),
                core::array::from_fn(|idx| &x[idx + 6]),
            )
        });
    }
}

/// Calls the finite complex-mass representation; its native definitions must
/// be registered in the same map as `sheet_exact`.
pub(crate) fn expression(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    call(
        "box",
        [p[0], p[1], p[2], p[3], p[4], p[5], m[0], m[1], m[2], m[3]],
    )
}

/// The nonzero-external construction requires at least one positive Källén
/// discriminant. Zero-external cases use the separate source reductions.
/// Vanishing channel invariants need analytic limits: the numerical source
/// perturbs them by a finite epsilon, which this exact representation avoids.
pub(crate) fn supported_kinematics(p: [&Atom; 6]) -> Atom {
    let all_nonzero = p[..4].iter().fold(Atom::num(1), |pred, value| {
        pred * if_nonzero(value, Atom::num(1))
    });
    let positive_lambda = any_nonzero(lambdas(p).map(|l| negative(&(-l))));
    if_nonzero(p[4], Atom::num(1))
        * if_nonzero(p[5], Atom::num(1))
        * if_nonzero_else(&all_nonzero, positive_lambda, Atom::num(1))
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::prelude::{Complex, Float, RealLike};

    #[test]
    fn composite_norms_keep_real_comparisons() {
        use symbolica::domains::float::{FloatLike, SingleFloat};
        use symbolica::domains::rational::Rational;

        let x = ["boxc_norm_z", "boxc_norm_w"].map(|name| symbolica::symbol!(name).to_atom());
        // These are the composite distances compared by the contour/root
        // choices. Each is real for arbitrary complex inputs, before any
        // model parameter or numerical tolerance is considered.
        let first = norm(&(Atom::num(1) - &x[0]));
        let second = norm(&(&x[1] - &x[0]));
        let difference = &first - &second;
        let outputs = [
            grouped(&re(&grouped(&x[0]))),
            grouped(&im(&grouped(&x[0]))),
            norm(&x[0]),
            first,
            second,
            difference.clone(),
            sign_nonnegative(&difference),
            sign_nonnegative(&(-&difference)),
            sign_nonnegative(&Atom::num(0)),
        ];
        let mut map = FunctionMap::new();
        sheet_exact::register(&mut map);
        // `norm` is a registered function of this module.
        register_norm(&mut map);
        let exact = Atom::evaluator_multiple(&outputs.each_ref().map(Atom::as_view), &x)
            .function_map(map.into())
            .direct_translation(true)
            .build()
            .unwrap();
        for bits in [128, 266] {
            let prototype = Float::new(bits);
            let mut evaluator = exact.clone().map_coeff_with_prec(
                &|value| {
                    Complex::new(
                        prototype.from_rational(&value.re),
                        prototype.from_rational(&value.im),
                    )
                },
                bits,
            );
            let mut cases = Vec::new();
            for n in 1..=40 {
                let ratio = |n, d| prototype.from_rational(&Rational::from((n, d)));
                let z = Complex::new(ratio(n, 123), ratio(41 - n, 287));
                let w = Complex::new(ratio(n + 3, 11), -ratio(n, 13));
                let mut mixed_z = z.clone();
                mixed_z.re.set_precision(bits - 3);
                let mut mixed_w = w.clone();
                mixed_w.im.set_precision(bits - 2);
                cases.extend([
                    [z.clone(), w],
                    [z.clone(), z],
                    [mixed_z.clone(), mixed_w],
                    [mixed_z.clone(), mixed_z],
                ]);
            }
            // Actual bound operand from complex parity row233: 128/126 bits.
            // The second run also checks the same norm at 266/264 bits.
            let captured = Complex::new(
                Float::parse("1.0292061724814945285656539261050293211619088743365883613260552074468938284134112", Some(bits)).unwrap(),
                Float::parse("-4.4458254896654443057813999812257865754273155939406320738943731102215691718465756e-4", Some(bits - 2)).unwrap(),
            );
            cases.push([captured.clone(), captured]);
            for (case, inputs) in cases.into_iter().enumerate() {
                let mut values = core::array::from_fn::<_, 9, _>(|_| {
                    Complex::new(prototype.zero(), prototype.zero())
                });
                evaluator.evaluate(&inputs, &mut values);
                for value in &values[..6] {
                    assert!(value.im.is_zero(), "bits={bits} case={case}: {value:?}");
                }
                assert_eq!(values[0].re, inputs[0].re);
                assert_eq!(values[1].re, inputs[0].im);
                assert_eq!(values[2].re, inputs[0].norm_squared());
                let forward = if values[5].re < prototype.zero() {
                    -1
                } else {
                    1
                };
                let backward = if values[5].re > prototype.zero() {
                    -1
                } else {
                    1
                };
                for (value, expected) in values[6..].iter().zip([forward, backward, 1]) {
                    assert_eq!(
                        *value,
                        Complex::new(prototype.from_i64(expected), prototype.zero())
                    );
                }
            }
        }
    }

    #[test]
    fn principal_sqrt_complex_axis_derivatives() {
        crate::analytic_continuation_tests::check_complex_axis_germs(principal_sqrt, false, false);
    }

    #[test]
    fn contour_quadratic_roots_preserve_exact_axes() {
        let x = ["boxc_root_a", "boxc_root_b", "boxc_root_c"]
            .map(|name| symbolica::symbol!(name).to_atom());
        let mut map = FunctionMap::new();
        sheet_exact::register(&mut map);
        register_norm(&mut map);
        define(&mut map, "root_1", 3, |x| {
            root_bodies(&x[0], &x[1], &x[2])[0].clone()
        });
        define(&mut map, "root_2", 3, |x| {
            root_bodies(&x[0], &x[1], &x[2])[1].clone()
        });
        let values = roots(&x[0], &x[1], &x[2]);
        let mut evaluator = Atom::evaluator_multiple(&values.each_ref().map(Atom::as_view), &x)
            .function_map(map.into())
            .direct_translation(true)
            .build()
            .unwrap()
            .map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64()));
        for (input, expected) in [
            ([1., 0., 1.], [(0., 1.), (0., -1.)]),
            ([-1., 0., -1.], [(0., -1.), (0., 1.)]),
            ([1., 0., -1.], [(1., 0.), (-1., 0.)]),
            ([-1., 0., 1.], [(-1., 0.), (1., 0.)]),
            ([1., 2., 1.], [(-1., 0.), (-1., 0.)]),
        ] {
            let input = input.map(|v| Complex::new(v, 0.));
            let mut actual = [Complex::new(0., 0.); 2];
            evaluator.evaluate(&input, &mut actual);
            assert_eq!(
                actual,
                expected.map(|(re, im)| Complex::new(re, im)),
                "{input:?}"
            );
        }
    }

    #[test]
    fn complex_box_finite_continuation_regression() {
        super::super::tests::large_stack(|| {
            let x: [Atom; 10] = core::array::from_fn(|index| {
                symbolica::symbol!(format!("boxc_regression_{index}")).to_atom()
            });
            let integral = expression(
                core::array::from_fn(|index| &x[index]),
                core::array::from_fn(|index| &x[index + 6]),
            );
            let mut map = FunctionMap::new();
            sheet_exact::register(&mut map);
            register(&mut map);
            let exact = Atom::evaluator_multiple(&[integral.as_view()], &x)
                .function_map(map.into())
                .direct_translation(true)
                .build()
                .unwrap();
            let fixture = include_str!("../tests/data/scalar_audit.txt")
                .lines()
                .skip_while(|line| *line != "# 293 mass_masks_widths")
                .nth(1)
                .expect("audited complex-mass box fixture");
            let values = fixture
                .split_whitespace()
                .map(|v| v.parse::<f64>().unwrap())
                .collect::<Vec<_>>();
            assert_eq!(values.len(), 22);
            assert_eq!(values[0], 4.0);
            let base_input: [Complex<f64>; 10] = core::array::from_fn(|index| {
                if index < 6 {
                    Complex::new(values[index + 2], 0.0)
                } else {
                    Complex::new(values[8 + 2 * (index - 6)], values[9 + 2 * (index - 6)])
                }
            });
            let expected = Complex::new(values[16], values[17]);
            let mut double = exact
                .clone()
                .map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64()));
            let mut precise = exact.map_coeff_with_prec(
                &|v| Complex::new(v.re.to_multi_prec_float(128), v.im.to_multi_prec_float(128)),
                128,
            );
            for rotation in 0..4 {
                let input = core::array::from_fn::<_, 10, _>(|index| {
                    let source = if index < 4 {
                        (index + rotation) % 4
                    } else if index < 6 {
                        4 + (index - 4 + rotation) % 2
                    } else {
                        6 + (index - 6 + rotation) % 4
                    };
                    base_input[source]
                });
                let mut out = [Complex::new(0.0, 0.0)];
                double.evaluate(&input, &mut out);
                let difference = out[0] - expected;
                assert!(
                    difference.re.hypot(difference.im) < 3e-10,
                    "rotation={rotation}, f64 actual={:?}, expected={expected:?}",
                    out[0]
                );
                let precise_input = input
                    .map(|v| Complex::new(Float::with_val(128, v.re), Float::with_val(128, v.im)));
                let mut precise_out = [Complex::new(Float::new(128), Float::new(128))];
                precise.evaluate(&precise_input, &mut precise_out);
                let difference =
                    Complex::new(precise_out[0].re.to_f64(), precise_out[0].im.to_f64()) - expected;
                assert!(
                    difference.re.hypot(difference.im) < 3e-10,
                    "rotation={rotation}, 128-bit actual={:?}, expected={expected:?}",
                    precise_out[0]
                );
            }
        });
    }
}
