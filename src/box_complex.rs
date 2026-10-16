//! Native complex-mass box representation from OneLOop's `avh_olo_boxc`.
//!
//! The logarithms, dilogarithms, contour residues and permutation choices are
//! transparent Symbolica definitions. Infinitesimal root lips are carried by
//! `SheetAtom`; no finite regulator or numerical integral callback is used.
#![allow(clippy::too_many_arguments)] // Match the source's named polynomial coefficients.
use super::*;
use sheet_exact::{
    SheetAtom as Q, divided_difference as dd, i, im, negative, re, sign_nonnegative,
};

fn function(name: &str) -> Symbol {
    symbolica::symbol!(format!("__olo_boxc_{name}"))
}

fn call<const N: usize>(name: &str, args: [&Atom; N]) -> Atom {
    function(name).call(args.map(Clone::clone).as_slice())
}

fn norm(z: &Atom) -> Atom {
    z * z.conj()
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
    // Insert y + eps*y' + eps^2*y'' into K*y^2+L*y+N=0. The source
    // has mass widths at order EPSN and its explicit -IEPS*abs(Re f)
    // at order EPSN^2. Retaining both coefficients fixes cancellations
    // between equal mass widths without ever choosing a finite epsilon.
    let derivative = Atom::num(2) * k * &y + l;
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

fn r1_call(z: &Continued, y1: &Continued, y2: &Continued, fy: &Atom) -> Atom {
    call(
        "r1",
        [
            &z.value, &z.first, &z.second, &y1.value, &y1.first, &y1.second, &y2.value, &y2.first,
            &y2.second, fy,
        ],
    )
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
) -> Atom {
    call(
        "s3",
        [
            &y1.value, &y1.first, &y1.second, &y2.value, &y2.first, &y2.second, d, e, a, b, c, db,
            dc,
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

fn r0(x: [&Atom; 6]) -> Atom {
    let y1 = Continued::from_args([x[0], x[1], x[2]]);
    let y2 = Continued::from_args([x[3], x[4], x[5]]);
    (&y2.complement(0).sheet() / &y1.complement(0).sheet()).log_over_one_minus() / &y1.value
        + (&y2.complement(1).sheet() / &y1.complement(1).sheet()).log_over_one_minus()
            / (1 - &y1.value)
}

fn r1(x: [&Atom; 10]) -> Atom {
    let zz = Continued::from_args([x[0], x[1], x[2]]);
    let yy1 = Continued::from_args([x[3], x[4], x[5]]);
    let yy2 = Continued::from_args([x[6], x[7], x[8]]);
    let (z, y1, y2, fy) = (&zz.value, &yy1.value, &yy2.value, x[9]);
    let oz: Atom = 1 - z;
    let yz1 = y1 - z;
    let yz2 = y2 - z;
    let q1 = yy1.subtract(&zz).sheet();
    let q2 = yy2.subtract(&zz).sheet();
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
    let near_y2 = fy * q1.log() - r0_call(&yy2, &zz) * &ratio;
    let near_y1 = fy * q2.log() - r0_call(&yy1, &zz) * (&q2 / &q1).log_over_one_minus();
    let ordinary = (if_nonzero(&yz1, &yz1 * q1.log() * r0_call(&yy1, &zz))
        - if_nonzero(&yz2, &yz2 * q2.log() * r0_call(&yy2, &zz)))
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
    let dilog_z = if_nonzero_else(
        &zsmall,
        if_nonzero(z, dd(&(&qz / &q1), &(&qz / &q2)) * z / (&yz1 * &yz2)),
        dd(&(&q1 / &qz), &(&q2 / &qz)) / z,
    );
    // The first matching branch in r1fun determines which inversion is used.
    let dilog_o = if_nonzero_else(
        &((1 - &zsmall) * &osmall),
        if_nonzero(&oz, dd(&(&qo / &q1), &(&qo / &q2)) * &oz / (&yz1 * &yz2)),
        dd(&(&q1 / &qo), &(&q2 / &qo)) / oz,
    );
    logarithms + dilog_z + dilog_o
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

fn s3(x: [&Atom; 13]) -> Atom {
    let mut y1 = Continued::from_args([x[0], x[1], x[2]]);
    let mut y2 = Continued::from_args([x[3], x[4], x[5]]);
    let [d, e, a, b, c, db, dc] = [x[6], x[7], x[8], x[9], x[10], x[11], x[12]];
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
    let [root1, root2] = roots(&re(a), b, c);
    let rea = sign_nonnegative(&re(a));
    let mut z1 = continued_root(root1, a, b, db, dc, &Atom::num(0));
    let mut z2 = continued_root(root2, a, b, db, dc, &Atom::num(0));
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
    let z0 = -b / (Atom::num(2) * a);
    let real_coefficients = 1 - any_nonzero([im(a), im(b), im(c)]);
    let positive_a = negative(&(-re(a)));
    let interior_z = negative(&(-re(&z0))) * negative(&(re(&z0) - 1));
    let regular_poles = outside_closed_interval(&y1.value) * outside_closed_interval(&y2.value);
    let double_root = if_nonzero_else(
        &(b * b - Atom::num(4) * a * c),
        Atom::num(0),
        real_coefficients * positive_a * interior_z * regular_poles,
    );
    let quadratic = if_nonzero_else(
        &double_root,
        call("s3_double_real", [&y1.value, &y2.value, &z0, a]),
        ordinary_quadratic,
    );
    let mut z = continued_root(-c / b, &Atom::num(0), b, db, dc, &Atom::num(0));
    let known = any_nonzero([im(&z.value), im(&z.first), im(&z.second)]);
    z.second = if_nonzero_else(
        &known,
        z.second.clone(),
        -i() * &simc * sign_nonnegative(&re(b)) * Symbol::ABS.call((re(&z.value),)),
    );
    let sz = z.sign();
    let sb = imaginary_sign(db, &simc);
    let linear = &fy
        * (Q::with_sign(b.clone(), sb.clone()).log()
            + eta(b, &sb, &(-&z.value), &(-&sz), c, &simc))
        + r1_call(&z, &y1, &y2, &fy);
    let constant = if_nonzero(c, Q::with_sign(c.clone(), simc).log() * &fy);
    if_nonzero(
        e,
        if_nonzero_else(a, quadratic, if_nonzero_else(b, linear, constant)) / e,
    )
}

fn planar(x: [&Atom; 12]) -> Atom {
    let y1 = Continued::from_args([x[0], x[1], x[2]]);
    let y2 = Continued::from_args([x[3], x[4], x[5]]);
    let [p1, p2, a, b, c, dc] = [x[6], x[7], x[8], x[9], x[10], x[11]];
    let x1 = b * &y1.value + c;
    let x2 = b * &y2.value + c;
    let sheet = |y: &Continued, value: &Atom| {
        let first = b * &y.first + dc;
        let second = b * &y.second;
        Continued {
            value: a / value,
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

fn t1([a, c, g, h, d, e, f, j, dpe, dd, de, df, dj, ddpe]: [&Atom; 14]) -> Atom {
    let k = h * a - c * g;
    let l = h * d - c * j - e * g;
    let n = h * f - e * j;
    let dl = h * dd - c * dj - de * g;
    let dn = h * df - de * j - e * dj;
    let dn2 = -de * dj - i() * Symbol::ABS.call((re(f),)) * h;
    let [y1, y2] = roots(&k, &l, &n).map(|y| continued_root(y, &k, &l, &dl, &dn, &dn2));
    let zero = Atom::num(0);
    let one = Atom::num(1);
    (-s3_call(&y1, &y2, &zero, &one, &(a + c), dpe, f, ddpe, df)
        + s3_call(&y1, &y2, &zero, &one, &zero, &(g + h), j, &zero, dj)
        - s3_call(&y1, &y2, &zero, &one, &zero, g, j, &zero, dj)
        + s3_call(&y1, &y2, &zero, &one, a, d, f, dd, df))
        / k
}

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
    ]: [&Atom; 18],
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
    (-s3_call(&y1, &y2, &zero, &one, a, &(e + c), dpf, de, ddpf)
        + s3_call(&y1, &y2, &zero, &one, g, &(j + h), dpf, dj, ddpf)
        - s3_call(&y1, &y2, &zero, &one, &(g + h), dpj, f, ddpj, df)
        + s3_call(&y1, &y2, &zero, &one, &(a + c), dpe, f, ddpe, df))
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

fn beta_data(
    beta: &Atom,
    b: &Atom,
    c: &Atom,
    g: &Atom,
    h: &Atom,
    d: &Atom,
    e: &Atom,
    f: &Atom,
    j: &Atom,
    dd: &Atom,
    de: &Atom,
    df: &Atom,
    dj: &Atom,
) -> BetaData {
    let a1 = g + beta * h;
    let b1 = c + Atom::num(2) * beta * b;
    let c1 = d + beta * e;
    let k = b * &a1 - h * &b1;
    let l = e * &a1 - h * &c1 - j * &b1;
    let n = f * &a1 - j * &c1;
    let dc1 = dd + beta * de;
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

fn tfun_selected(
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
        dd,
        de,
        df,
        dj,
        ddpe,
        ddpf,
    ]: [&Atom; 18],
) -> Atom {
    let data = beta_data(beta, b, c, g, h, d, e, f, j, dd, de, df, dj);
    let [y1, y2] = &data.y;
    let zero = Atom::num(0);
    let one = Atom::num(1);
    let p1 = call("inside", [beta, &y1.value, &y1.first, &y1.second]);
    let p2 = call("inside", [beta, &y2.value, &y2.first, &y2.second]);
    let integral = s3_call(y1, y2, beta, &one, &zero, h, &(g + j), &zero, dj)
        - s3_call(y1, y2, &zero, &(1 - beta), &zero, &(g + h), j, &zero, dj)
        + s3_call(y1, y2, &zero, &(-beta), &zero, g, j, &zero, dj)
        - s3_call(y1, y2, beta, &one, b, &(c + e), &(a + dpf), de, ddpf)
        + s3_call(y1, y2, &zero, &(1 - beta), &(a + b + c), dpe, f, ddpe, df)
        - s3_call(y1, y2, &zero, &(-beta), a, d, f, dd, df);
    let residue = call(
        "planar",
        [
            &y1.value, &y1.first, &y1.second, &y2.value, &y2.first, &y2.second, &p1, &p2, &data.a,
            &data.b, &data.c, &data.dc,
        ],
    );
    (integral + if_nonzero_else(&less(&Atom::num(0), &im(beta)), -&residue, residue)) / data.k
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
    ]: [&Atom; 17],
) -> Atom {
    let sj = imaginary_sign(jin, &imaginary_sign(djin, &Atom::num(-1)));
    let g = -&sj * gin;
    let h = -&sj * hin;
    let j = -&sj * jin;
    let dj = -&sj * djin;
    let [root1, root2] = roots(b, c, a);
    let swap = less(&norm(&root2), &norm(&root1));
    let beta1 = if_nonzero_else(&swap, root2.clone(), root1.clone());
    let beta2 = if_nonzero_else(&swap, root1, root2);
    let data1 = beta_data(&beta1, b, c, &g, &h, d, e, f, &j, dd, de, df, &dj);
    let data2 = beta_data(&beta2, b, c, &g, &h, d, e, f, &j, dd, de, df, &dj);
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
    // The contour opposite the pole nearest a boundary is used. If neither
    // contour encloses a pole, the first (smaller-magnitude) beta is selected.
    let use_first = if_nonzero_else(&(&d1 - 1), less(&d2, &d1), Atom::num(1));
    let beta = if_nonzero_else(&use_first, beta1, beta2);
    let general = call(
        "tfun_selected",
        [
            a, b, c, &g, &h, d, e, f, &j, dpe, dpf, &beta, dd, de, df, &dj, ddpe, ddpf,
        ],
    );
    let no_b = call(
        "t1",
        [a, c, &g, &h, d, e, f, &j, dpe, dd, de, df, &dj, ddpe],
    );
    let no_a = call(
        "t1",
        [
            &(b + c),
            &(-c),
            &(-&g - &h),
            &g,
            &(-dpe - Atom::num(2) * (b + c)),
            &(d + c),
            &(dpe + b + c + f),
            &(&g + &h + &j),
            &(-e - Atom::num(2) * b - c),
            &(-ddpe),
            dd,
            &(ddpe + df),
            &dj,
            &(-de),
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
    let d = m[2] - m[3] - p[2];
    let e = m[0] - m[2] + p[2] - p[3];
    let f = m[3];
    let k = m[1] - m[2] - p[5] + p[2];
    let dpe = m[0] - m[3] - p[3];
    let dpk = m[1] - m[3] - p[5];
    let dpf = m[2] - p[2];
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
        ],
    )
}

fn no_opposite_chart_regular(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    let a = p[4] - p[0];
    let c = p[3] - p[4] - p[2];
    let g = p[1];
    let h = p[5] - p[1] - p[2];
    let d = m[2] - m[3] - p[2];
    let e = m[0] - m[2] + p[2] - p[3];
    let j = m[1] - m[2] - p[5] + p[2];
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

fn zero_external(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    let dm = mass_directions(m);
    let a = p[2];
    let b = p[1];
    let c = p[5] - p[1] - p[2];
    let h = p[3] - p[4] - p[5] + p[1];
    let j = p[4] - p[0] - p[1];
    let d = m[2] - m[3] - p[2];
    let e = m[1] - m[2] - p[5] + p[2];
    let k = m[0] - m[1] + p[5] - p[3];
    let f = m[3];
    let cph = p[3] - p[4] - p[2];
    let dpe = m[1] - m[3] - p[5];
    let epk = m[0] - m[2] + p[2] - p[3];
    let dek = m[0] - m[3] - p[3];
    let dpf = m[2] - p[2];
    call(
        "tfun",
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
        ],
    ) - call(
        "tfun",
        [
            a,
            &(b + &j),
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
        ],
    )
}

fn ordered(p: [&Atom; 6], m: [&Atom; 4]) -> Atom {
    let dm = mass_directions(m);
    let a = p[2];
    let b = p[1];
    let g = p[0];
    let c = p[5] - p[1] - p[2];
    let h = p[3] - p[4] - p[5] + p[1];
    let j = p[4] - p[0] - p[1];
    let d = m[2] - m[3] - p[2];
    let e = m[1] - m[2] - p[5] + p[2];
    let k = m[0] - m[1] + p[5] - p[3];
    let f = m[3];
    let abc = p[5];
    let bgj = p[4];
    let jph = p[3] - p[0] - p[5];
    let cph = p[3] - p[4] - p[2];
    let dpe = m[1] - m[3] - p[5];
    let epk = m[0] - m[2] + p[2] - p[3];
    let dek = m[0] - m[3] - p[3];
    let dpf = m[2] - p[2];
    let def = m[1] - p[5];
    let [x1, x2] = roots(g, &j, b);
    let alpha = if_nonzero_else(&less(&norm(&re(&x1)), &norm(&re(&x2))), x2, x1);
    let o1 = 1 - &alpha;
    let j1 = &j + Atom::num(2) * g * &alpha;
    let e1 = &e + &k * &alpha;
    let de1 = &dm[1] - &dm[2] + (&dm[0] - &dm[1]) * &alpha;
    -call(
        "tfun",
        [
            abc,
            g,
            &jph,
            &(&c + Atom::num(2) * b + (&h + &j) * &alpha),
            &j1,
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
        ],
    ) + &o1
        * call(
            "tfun",
            [
                a,
                bgj,
                &cph,
                &(&c + &h * &alpha),
                &(&o1 * &j1),
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
            ],
        )
        + &alpha
            * call(
                "tfun",
                [
                    a,
                    b,
                    &c,
                    &(&c + &h * &alpha),
                    &(-&j1 * &alpha),
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

pub(crate) fn register(map: &mut FunctionMap) {
    define(map, "root_1", 3, |x| {
        root_bodies(&x[0], &x[1], &x[2])[0].clone()
    });
    define(map, "root_2", 3, |x| {
        root_bodies(&x[0], &x[1], &x[2])[1].clone()
    });
    define(map, "r0", 6, |x| r0(core::array::from_fn(|idx| &x[idx])));
    define(map, "r1", 10, |x| r1(core::array::from_fn(|idx| &x[idx])));
    define(map, "s3", 13, |x| s3(core::array::from_fn(|idx| &x[idx])));
    define(map, "s3_double_real", 4, |x| {
        s3_double_real(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "planar", 12, |x| {
        planar(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "inside", 4, |x| {
        inside(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "distance", 4, |x| {
        distance(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "t1", 14, |x| t1(core::array::from_fn(|idx| &x[idx])));
    define(map, "t13", 18, |x| t13(core::array::from_fn(|idx| &x[idx])));
    define(map, "tfun_selected", 18, |x| {
        tfun_selected(core::array::from_fn(|idx| &x[idx]))
    });
    define(map, "tfun", 17, |x| {
        tfun(core::array::from_fn(|idx| &x[idx]))
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
    fn principal_sqrt_complex_axis_derivatives() {
        crate::analytic_continuation_tests::check_complex_axis_germs(principal_sqrt, false, false);
    }

    #[test]
    fn contour_quadratic_roots_preserve_exact_axes() {
        let x = ["boxc_root_a", "boxc_root_b", "boxc_root_c"]
            .map(|name| symbolica::symbol!(name).to_atom());
        let mut map = FunctionMap::new();
        sheet_exact::register(&mut map);
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
