//! Binary64 port of the contour box `box_complex.rs` (OneLOop `boxc`:
//! Dao Thi Nhung & Le Duc Ninh, 't Hooft & Veltman), branch for branch.
//!
//! Function names follow the registered `__olo_boxc_*` definitions. Every
//! exact zero test, comparison, infinitesimal direction ("Continued" first and
//! second coefficients) and sheet choice is the expression's; only the
//! evaluation order differs (inactive branches are never evaluated).
#![allow(clippy::too_many_arguments)] // The source's named polynomial coefficients.
use super::super::primitives::{C, NativeFloat, inverse};
use super::Constants;
use super::sheet::{Sheet, dd};

pub(super) type Z = C<f64>;

const ZERO: Z = Z { re: 0.0, im: 0.0 };
const ONE: Z = Z { re: 1.0, im: 0.0 };

#[inline(always)]
pub(super) fn r(x: f64) -> Z {
    Z { re: x, im: 0.0 }
}

#[inline(always)]
pub(super) fn dv(a: Z, b: Z) -> Z {
    a * inverse(&b)
}

#[inline(always)]
pub(super) fn nz(a: Z) -> bool {
    a.re != 0.0 || a.im != 0.0
}

#[inline(always)]
fn norm(a: Z) -> f64 {
    a.re * a.re + a.im * a.im
}

/// `coefficient_sum`: the left-to-right sum of `terms`, or the invariant
/// `exact` it equals when the sum cancels by more than 2^16.
#[inline]
fn coefficient_sum(terms: &[Z], exact: Z) -> Z {
    let numeric = terms.iter().fold(ZERO, |sum, &term| sum + term);
    let scale: f64 = terms.iter().map(|&term| norm(term)).sum();
    if negative(4294967296.0 * norm(exact) - scale) { exact } else { numeric }
}

/// `negative` of a real value: `x - |x|` is nonzero (x < 0, inf or NaN).
#[inline(always)]
pub(super) fn negative(x: f64) -> bool {
    x - x.abs() != 0.0
}

/// `sign_nonnegative` of a real value.
#[inline(always)]
fn sgn(x: f64) -> i32 {
    if negative(x) { -1 } else { 1 }
}

/// `imaginary_sign(z, fallback)`.
#[inline(always)]
fn isign(z: Z, fallback: i32) -> i32 {
    if z.im != 0.0 { sgn(z.im) } else { fallback }
}

/// A value and the first two coefficients of its infinitesimal continuation.
#[derive(Clone, Copy, Debug)]
struct Cont {
    v: Z,
    f: Z,
    s: Z,
}

impl Cont {
    #[inline]
    fn sign(&self) -> i32 {
        isign(self.v, isign(self.f, isign(self.s, 1)))
    }

    #[inline]
    fn sheet(&self) -> Sheet<f64> {
        Sheet::with_sign(self.v, self.sign())
    }

    #[inline]
    fn complement(&self, constant: f64) -> Self {
        Self {
            v: r(constant) - self.v,
            f: -self.f,
            s: -self.s,
        }
    }

    #[inline]
    fn minus(&self, other: &Self) -> Self {
        Self {
            v: self.v - other.v,
            f: self.f - other.f,
            s: self.s - other.s,
        }
    }
}

/// `mass_directions`: -i*|Re m| for exactly real masses.
#[inline]
fn mass_directions(m: [Z; 4]) -> [Z; 4] {
    m.map(|mass| {
        if mass.im != 0.0 {
            ZERO
        } else {
            Z {
                re: 0.0,
                im: -mass.re.abs(),
            }
        }
    })
}

#[inline]
fn continued_root(y: Z, k: Z, l: Z, dl: Z, dn: Z, dn2: Z) -> Cont {
    let derivative = r(2.0) * k * y + l;
    let inv = inverse(&derivative);
    let first = -((dl * y + dn) * inv);
    let second = -((k * first * first + dl * first + dn2) * inv);
    Cont {
        v: y,
        f: first,
        s: second,
    }
}

/// `principal_sqrt`: exact axes and rectangular components off the axis.
#[inline]
fn principal_sqrt(v: Z) -> Z {
    if v.im == 0.0 {
        if negative(v.re) {
            Z {
                re: 0.0,
                im: (-v.re).sqrt(),
            }
        } else {
            Z {
                re: v.re.sqrt(),
                im: v.im,
            }
        }
    } else {
        let absolute = <f64 as NativeFloat>::hypot(&v.re, &v.im);
        if negative(v.re) {
            let w = sgn(v.im) as f64 * ((absolute - v.re) / 2.0).sqrt();
            Z {
                re: v.im / (2.0 * w),
                im: w,
            }
        } else {
            let u = ((absolute + v.re) / 2.0).sqrt();
            Z {
                re: u,
                im: v.im / (2.0 * u),
            }
        }
    }
}

/// `root_1`, `root_2` (solabc labels and product relation).
#[inline]
fn roots(a: Z, b: Z, c: Z) -> [Z; 2] {
    if !nz(c) {
        return [if nz(a) { -dv(b, a) } else { ZERO }, ZERO];
    }
    if !nz(a) {
        let linear = if nz(b) { -dv(c, b) } else { ZERO };
        return [linear, linear];
    }
    let delta = principal_sqrt(b * b - r(4.0) * a * c);
    let plus = -b + delta;
    let minus = -b - delta;
    if negative(norm(plus) - norm(minus)) {
        [dv(r(2.0) * c, minus), dv(minus, r(2.0) * a)]
    } else {
        [dv(plus, r(2.0) * a), dv(r(2.0) * c, plus)]
    }
}

/// `eta(a, sa, b, sb, c, sc)` = i*pi*sc*(1+sa*sb)*(1-sa*sc)/2.
#[inline]
fn eta(a: Z, sa: i32, b: Z, sb: i32, c: Z, sc: i32, k: &Constants<f64>) -> Z {
    let sa = isign(a, sa);
    let sb = isign(b, sb);
    let sc = isign(c, sc);
    let value = sc * (1 + sa * sb) * (1 - sa * sc) / 2;
    if value == 0 {
        ZERO
    } else {
        Z {
            re: 0.0,
            im: k.pi * value as f64,
        }
    }
}

#[inline]
fn lom(q: &Sheet<f64>, k: &Constants<f64>) -> Z {
    q.log_over_one_minus(k)
}

fn r0(y1: &Cont, y2: &Cont, k: &Constants<f64>) -> Z {
    // `r0`: the windings of the two quotients are taken out and added as
    // i pi (n0 - n1)/(y1 - y2), with the exact pole difference.
    let split = |q: Sheet<f64>| {
        let winding = if q.odd { 0 } else { q.phase };
        (
            Sheet {
                phase: q.phase - winding,
                ..q
            },
            winding,
        )
    };
    let (q0, n0) = split(y2.complement(0.0).sheet().div(&y1.complement(0.0).sheet()));
    let (q1, n1) = split(y2.complement(1.0).sheet().div(&y1.complement(1.0).sheet()));
    let value = dv(lom(&q0, k), y1.v) + dv(lom(&q1, k), ONE - y1.v);
    let winding = n0 - n1;
    if winding != 0 {
        let pi_n = k.pi * winding as f64;
        value + dv(Z { re: 0.0, im: pi_n }, y1.v - y2.v)
    } else {
        value
    }
}

fn r1(zz: &Cont, yy1: &Cont, yy2: &Cont, fy: Z, k: &Constants<f64>) -> Z {
    let (z, y1, y2) = (zz.v, yy1.v, yy2.v);
    let oz = ONE - z;
    let yz1 = y1 - z;
    let yz2 = y2 - z;
    let q1 = yy1.minus(zz).sheet();
    let q2 = yy2.minus(zz).sheet();
    let qz = zz.complement(0.0).sheet();
    let qo = zz.complement(1.0).sheet();
    let h12 = norm(y1 - y2);
    let hz1 = norm(yz1);
    let hz2 = norm(yz2);
    let hzz = norm(z);
    let hoz = norm(oz);
    let zsmall = negative(hzz - hz1) && negative(hzz - hz2) && negative(hzz - hoz);
    let osmall = negative(hoz - hz1) && negative(hoz - hz2);
    let logarithms = if zsmall || osmall {
        let product_log = q1.log(k) + q2.log(k);
        let y2sq = y2 * y2;
        let inv_y2sq = inverse(&y2sq);
        let yratio = Cont {
            v: dv(y2 - ONE, y2),
            f: yy2.f * inv_y2sq,
            s: yy2.s * inv_y2sq - yy2.f * yy2.f * inverse(&(y2sq * y2)),
        };
        let ylog = yratio.sheet().log(k);
        let ratio = lom(&q1.div(&q2), k);
        let half = product_log * 0.5;
        let bracket = if zsmall {
            half + ylog - qo.log(k)
        } else {
            -half + ylog + qz.log(k)
        };
        fy * q1.log(k) - dv(bracket * ratio, yz2)
    } else if !negative(hz2 - h12) && !negative(hz1 - hz2) {
        fy * q1.log(k) - r0(yy2, zz, k) * lom(&q1.div(&q2), k)
    } else if !negative(hz1 - h12) && !negative(hz2 - hz1) {
        fy * q2.log(k) - r0(yy1, zz, k) * lom(&q2.div(&q1), k)
    } else {
        let first = if nz(yz1) {
            yz1 * q1.log(k) * r0(yy1, zz, k)
        } else {
            ZERO
        };
        let second = if nz(yz2) {
            yz2 * q2.log(k) * r0(yy2, zz, k)
        } else {
            ZERO
        };
        dv(first - second, y1 - y2)
    };
    let dilog_z = if zsmall {
        if nz(z) {
            dv(dd(&qz.div(&q1), &qz.div(&q2), k) * z, yz1 * yz2)
        } else {
            ZERO
        }
    } else {
        dv(dd(&q1.div(&qz), &q2.div(&qz), k), z)
    };
    let dilog_o = if !zsmall && osmall {
        if nz(oz) {
            dv(dd(&qo.div(&q1), &qo.div(&q2), k) * oz, yz1 * yz2)
        } else {
            ZERO
        }
    } else {
        dv(dd(&q1.div(&qo), &q2.div(&qo), k), oz)
    };
    logarithms + dilog_z + dilog_o
}

fn s3_double_real(y1: Z, y2: Z, z: Z, a: Z, k: &Constants<f64>) -> Z {
    let a1 = z - y1;
    let a2 = z - y2;
    let oz = ONE - z;
    let q01 = Sheet::upper(-y1).div(&Sheet::upper(a1));
    let q02 = Sheet::upper(-y2).div(&Sheet::upper(a2));
    let q11 = Sheet::upper(ONE - y1).div(&Sheet::upper(a1));
    let q12 = Sheet::upper(ONE - y2).div(&Sheet::upper(a2));
    let log = |x: Z| super::super::primitives::log(&x);
    let jdd = dv(
        z * dd(&q01, &q02, k) + oz * dd(&q11, &q12, k)
            - dv(z * log(z) * lom(&q01.div(&q02), k), q02.value())
            - dv(oz * log(oz) * lom(&q11.div(&q12), k), q12.value()),
        a1 * a2,
    );
    let pole = |y: Z| Cont {
        v: y,
        f: ZERO,
        s: ZERO,
    };
    r0(&pole(y1), &pole(y2), k) * log(a) + r(2.0) * jdd
}

#[inline]
fn outside_closed_interval(y: Z) -> bool {
    y.im != 0.0 || negative(y.re) || negative(1.0 - y.re)
}

/// `s3`: the y1, y2 poles are mapped by y -> (d + y)/e.
fn s3(y1: &Cont, y2: &Cont, d: Z, e: Z, a: Z, b: Z, c: Z, db: Z, dc: Z, k: &Constants<f64>) -> Z {
    if !nz(e) {
        return ZERO;
    }
    let simc = isign(c, isign(dc, isign(b, isign(db, -1))));
    let inv_e = inverse(&e);
    let map = |y: &Cont| Cont {
        v: (d + y.v) * inv_e,
        f: y.f * inv_e,
        s: y.s * inv_e,
    };
    let y1 = map(y1);
    let y2 = map(y2);
    let fy = r0(&y1, &y2, k);
    let value = if nz(a) {
        let z0 = -dv(b, r(2.0) * a);
        // `s3`'s discriminant, zero when it does not change 1024 |b|².
        let reference = 1024.0 * norm(b);
        let discriminant = b * b - r(4.0) * a * c;
        let double = !nz(Z {
            re: (reference + discriminant.re) - reference,
            im: discriminant.im,
        });
        let real_coefficients = a.im == 0.0 && b.im == 0.0 && c.im == 0.0;
        let double_root = double
            && real_coefficients
            && negative(-a.re)
            && negative(-z0.re)
            && negative(z0.re - 1.0)
            && outside_closed_interval(y1.v)
            && outside_closed_interval(y2.v);
        if double_root {
            s3_double_real(y1.v, y2.v, z0, a, k)
        } else {
            let [root1, root2] = roots(r(a.re), b, c);
            let rea = sgn(a.re);
            let mut z1 = continued_root(root1, a, b, db, dc, ZERO);
            let mut z2 = continued_root(root2, a, b, db, dc, ZERO);
            if double {
                // An exact double root splits like sqrt(eps), one root up and
                // one down, dominating every first-order shift of y1 and y2.
                let size = |z: Z| z.re.hypot(z.im);
                let dominant = 2.0 * (1.0 + size(y1.f) + size(y1.s) + size(y2.f) + size(y2.s));
                let centre = z0;
                z1 = Cont { v: centre, f: Z { re: 0.0, im: dominant }, s: ZERO };
                z2 = Cont { v: centre, f: Z { re: 0.0, im: -dominant }, s: ZERO };
            }
            let fallback1 = simc * rea * sgn(z2.v.re);
            let fallback2 = simc * rea * sgn(z1.v.re);
            for (z, fallback) in [(&mut z1, fallback1), (&mut z2, fallback2)] {
                let known = z.v.im != 0.0 || z.f.im != 0.0 || z.s.im != 0.0;
                if !known {
                    z.s = Z {
                        re: 0.0,
                        im: fallback as f64 * z.v.re.abs(),
                    };
                }
            }
            let s1 = z1.sign();
            let s2 = z2.sign();
            fy * (Sheet::with_sign(a, simc).log(k)
                + eta(-z1.v, -s1, -z2.v, -s2, ZERO, simc * rea, k))
                + r1(&z1, &y1, &y2, fy, k)
                + r1(&z2, &y1, &y2, fy, k)
        }
    } else if nz(b) {
        let mut z = continued_root(-dv(c, b), ZERO, b, db, dc, ZERO);
        let known = z.v.im != 0.0 || z.f.im != 0.0 || z.s.im != 0.0;
        if !known {
            z.s = Z {
                re: 0.0,
                im: -(simc * sgn(b.re)) as f64 * z.v.re.abs(),
            };
        }
        let sz = z.sign();
        let sb = isign(db, simc);
        fy * (Sheet::with_sign(b, sb).log(k) + eta(b, sb, -z.v, -sz, c, simc, k))
            + r1(&z, &y1, &y2, fy, k)
    } else if nz(c) {
        Sheet::with_sign(c, simc).log(k) * fy
    } else {
        ZERO
    };
    value * inv_e
}

/// `planar`; `real` marks exactly real logarithm arguments
/// (`real_intersections_absent`).
fn planar(y1: &Cont, y2: &Cont, p1: bool, p2: bool, a: Z, b: Z, c: Z, dc: Z, real: bool, k: &Constants<f64>) -> Z {
    if !p1 && !p2 {
        return ZERO;
    }
    let x1 = b * y1.v + c;
    let x2 = b * y2.v + c;
    let sheet = |y: &Cont, value: Z| {
        let first = b * y.f + dc;
        let second = b * y.s;
        let v2 = value * value;
        let ratio = dv(a, value);
        Cont {
            v: if real { r(ratio.re) } else { ratio },
            f: -dv(a * first, v2),
            s: dv(a * first * first, v2 * value) - dv(a * second, v2),
        }
        .sheet()
    };
    let value = if p1 && p2 {
        let q1 = sheet(y1, x1);
        let q2 = sheet(y2, x2);
        dv(lom(&q2.div(&q1), k) * b, x2)
    } else if p1 {
        dv(sheet(y1, x1).log(k), y1.v - y2.v)
    } else {
        dv(sheet(y2, x2).log(k), y2.v - y1.v)
    };
    // 2*i*pi*value
    Z {
        re: -2.0 * k.pi * value.im,
        im: 2.0 * k.pi * value.re,
    }
}

fn t1(
    [a, c, g, h, d, e, f, j, dpe, dd_, de, df, dj, ddpe, ac]: [Z; 15],
    k: &Constants<f64>,
) -> Z {
    let kk = h * a - c * g;
    let l = h * d - c * j - e * g;
    let n = h * f - e * j;
    let dl = h * dd_ - c * dj - de * g;
    let dn = h * df - de * j - e * dj;
    let dn2 = -de * dj - Z { re: 0.0, im: f.re.abs() } * h;
    let [r1_, r2_] = roots(kk, l, n);
    let y1 = continued_root(r1_, kk, l, dl, dn, dn2);
    let y2 = continued_root(r2_, kk, l, dl, dn, dn2);
    dv(
        -s3(&y1, &y2, ZERO, ONE, coefficient_sum(&[a, c], ac), dpe, f, ddpe, df, k)
            + s3(&y1, &y2, ZERO, ONE, ZERO, g + h, j, ZERO, dj, k)
            - s3(&y1, &y2, ZERO, ONE, ZERO, g, j, ZERO, dj, k)
            + s3(&y1, &y2, ZERO, ONE, a, d, f, dd_, df, k),
        kk,
    )
}

fn t13(
    [a, c, g, h, d, e, f, j, dpe, dpj, dpf, dd_, de, df, dj, ddpe, ddpj, ddpf]: [Z; 18],
    k: &Constants<f64>,
) -> Z {
    let kk = h * a - c * g;
    let l = a * d + h * e - d * g - c * j;
    let n = d * (e - j) + (h - c) * f;
    let dl = (a - g) * dd_ + h * de - c * dj;
    let dn = dd_ * (e - j) + d * (de - dj) + (h - c) * df;
    let dn2 = dd_ * (de - dj) - Z { re: 0.0, im: f.re.abs() } * (h - c);
    let [r1_, r2_] = roots(kk, l, n);
    let y1 = continued_root(r1_, kk, l, dl, dn, dn2);
    let y2 = continued_root(r2_, kk, l, dl, dn, dn2);
    dv(
        -s3(&y1, &y2, ZERO, ONE, a, e + c, dpf, de, ddpf, k)
            + s3(&y1, &y2, ZERO, ONE, g, j + h, dpf, dj, ddpf, k)
            - s3(&y1, &y2, ZERO, ONE, g + h, dpj, f, ddpj, df, k)
            + s3(&y1, &y2, ZERO, ONE, a + c, dpe, f, ddpe, df, k),
        kk,
    )
}

/// Barycentric coordinates in the triangle (0, -beta, 1-beta).
#[inline]
fn coordinates(beta: Z, y: Z) -> [f64; 2] {
    let sum = -y.im / beta.im;
    let x2 = y.re + beta.re * sum;
    let x1 = sum - x2;
    [x1, x2]
}

#[inline]
fn continued_nonnegative(value: f64, first: f64, second: f64) -> bool {
    if value != 0.0 {
        !negative(value)
    } else if first != 0.0 {
        !negative(first)
    } else {
        !negative(second)
    }
}

fn inside(beta: Z, y: &Cont) -> bool {
    if beta.im == 0.0 {
        return false;
    }
    let [x1, x2] = coordinates(beta, y.v);
    let [d1, d2] = coordinates(beta, y.f);
    let [dd1, dd2] = coordinates(beta, y.s);
    continued_nonnegative(x1, d1, dd1)
        && continued_nonnegative(x2, d2, dd2)
        && continued_nonnegative(1.0 - x1 - x2, -d1 - d2, -dd1 - dd2)
}

fn distance(beta: Z, y: &Cont) -> f64 {
    if inside(beta, y) {
        let [x1, x2] = coordinates(beta, y.v);
        if negative(x1 - x2) { x1 } else { x2 }
    } else {
        1.0
    }
}

struct BetaData {
    a: Z,
    b: Z,
    c: Z,
    dc: Z,
    k: Z,
    y: [Cont; 2],
}

fn beta_data(beta: Z, b: Z, c: Z, g: Z, h: Z, d: Z, e: Z, f: Z, j: Z, dd_: Z, de: Z, df: Z, dj: Z) -> BetaData {
    let a1 = g + beta * h;
    let b1 = c + r(2.0) * beta * b;
    let c1 = d + beta * e;
    let kk = b * a1 - h * b1;
    let l = e * a1 - h * c1 - j * b1;
    let n = f * a1 - j * c1;
    let dc1 = dd_ + beta * de;
    let dl = de * a1 - h * dc1 - dj * b1;
    let dn = df * a1 - dj * c1 - j * dc1;
    let dn2 = -dj * dc1 - Z { re: 0.0, im: f.re.abs() } * a1;
    let [r1_, r2_] = roots(kk, l, n);
    BetaData {
        a: a1,
        b: b1,
        c: c1,
        dc: dc1,
        k: kk,
        y: [
            continued_root(r1_, kk, l, dl, dn, dn2),
            continued_root(r2_, kk, l, dl, dn, dn2),
        ],
    }
}

/// `real_intersections_absent`: real coefficients whose line and conic meet at
/// complex conjugate points (`n^T adj(M) n > 0`).
fn real_intersections_absent(coefficients: [Z; 9]) -> bool {
    if coefficients.iter().any(|z| z.im != 0.0) {
        return false;
    }
    let [a, b, c, g, h, d, e, f, j] = coefficients.map(|z| z.re);
    let adj11 = b * f - e * e / 4.0;
    let adj22 = a * f - d * d / 4.0;
    let adj33 = a * b - c * c / 4.0;
    let adj12 = d * e / 4.0 - c * f / 2.0;
    let adj13 = c * e / 4.0 - b * d / 2.0;
    let adj23 = c * d / 4.0 - a * e / 2.0;
    let delta = g * g * adj11
        + h * h * adj22
        + j * j * adj33
        + 2.0 * (g * h * adj12 + g * j * adj13 + h * j * adj23);
    negative(-delta)
}

fn tfun_selected(
    a: Z, b: Z, c: Z, g: Z, h: Z, d: Z, e: Z, f: Z, j: Z, dpe: Z, dpf: Z, beta: Z,
    data: &BetaData, dd_: Z, de: Z, df: Z, dj: Z, ddpe: Z, ddpf: Z, abc: Z, k: &Constants<f64>,
) -> Z {
    let [y1, y2] = &data.y;
    let p1 = inside(beta, y1);
    let p2 = inside(beta, y2);
    let integral = s3(y1, y2, beta, ONE, ZERO, h, g + j, ZERO, dj, k)
        - s3(y1, y2, ZERO, ONE - beta, ZERO, g + h, j, ZERO, dj, k)
        + s3(y1, y2, ZERO, -beta, ZERO, g, j, ZERO, dj, k)
        - s3(y1, y2, beta, ONE, b, c + e, a + dpf, de, ddpf, k)
        + s3(y1, y2, ZERO, ONE - beta, coefficient_sum(&[a, b, c], abc), dpe, f, ddpe, df, k)
        - s3(y1, y2, ZERO, -beta, a, d, f, dd_, df, k);
    let real = real_intersections_absent([a, b, c, g, h, d, e, f, j]);
    let residue = planar(y1, y2, p1, p2, data.a, data.b, data.c, data.dc, real, k);
    let residue = if negative(0.0 - beta.im) { -residue } else { residue };
    dv(integral + residue, data.k)
}

fn tfun(
    [a, b, c, gin, hin, d, e, f, jin, dpe, dpf, dd_, de, df, djin, ddpe, ddpf, abc]: [Z; 18],
    k: &Constants<f64>,
) -> Z {
    let sj = isign(jin, isign(djin, -1));
    let s = r(-(sj as f64));
    let g = s * gin;
    let h = s * hin;
    let j = s * jin;
    let dj = s * djin;
    let value = if nz(b) {
        if nz(a) {
            let [root1, root2] = roots(b, c, a);
            let swap = negative(norm(root2) - norm(root1));
            let (beta1, beta2) = if swap { (root2, root1) } else { (root1, root2) };
            let data1 = beta_data(beta1, b, c, g, h, d, e, f, j, dd_, de, df, dj);
            let data2 = beta_data(beta2, b, c, g, h, d, e, f, j, dd_, de, df, dj);
            let min = |x: f64, y: f64| if negative(x - y) { x } else { y };
            let d1 = min(distance(beta1, &data1.y[0]), distance(beta1, &data1.y[1]));
            let d2 = min(distance(beta2, &data2.y[0]), distance(beta2, &data2.y[1]));
            let use_first = if d1 - 1.0 != 0.0 { negative(d2 - d1) } else { true };
            let (beta, data) = if use_first { (beta1, data1) } else { (beta2, data2) };
            tfun_selected(a, b, c, g, h, d, e, f, j, dpe, dpf, beta, &data, dd_, de, df, dj, ddpe, ddpf, abc, k)
        } else {
            // a = 0: b + c = abc, and the new a + c is b.
            let bc = coefficient_sum(&[b, c], abc);
            t1(
                [
                    bc,
                    -c,
                    -g - h,
                    g,
                    -dpe - r(2.0) * bc,
                    d + c,
                    dpe + bc + f,
                    g + h + j,
                    -e - r(2.0) * b - c,
                    -ddpe,
                    dd_,
                    ddpe + df,
                    dj,
                    -de,
                    b,
                ],
                k,
            )
        }
    } else {
        // b = 0: a + c = abc.
        t1([a, c, g, h, d, e, f, j, dpe, dd_, de, df, dj, ddpe, abc], k)
    };
    -(r(sj as f64) * value)
}

/// `opposite_regular`: the auxiliary poles avoid the integration endpoints.
fn opposite_regular(p: [Z; 6], m: [Z; 4]) -> bool {
    if nz(p[0]) || nz(p[2]) {
        return false;
    }
    let a = p[4] - p[0];
    let c = p[3] - p[4] - p[2];
    let g = p[1];
    let h = p[5] - p[1] - p[2];
    let d = m[2] - m[3] - p[2];
    let e = m[0] - m[2] + p[2] - p[3];
    let j = m[1] - m[2] - p[5] + p[2];
    let kk = h * a - c * g;
    let l = a * d + h * e - d * g - c * j;
    let n = d * (e - j) + (h - c) * m[3];
    nz(n * (kk + l + n))
}

fn opposite_chart(p: [Z; 6], m: [Z; 4], k: &Constants<f64>) -> Z {
    let dm = mass_directions(m);
    let a = p[4] - p[0];
    let c = p[3] - p[4] - p[2];
    let g = p[1];
    let h = p[5] - p[1] - p[2];
    let d = m[2] - m[3] - p[2];
    let e = m[0] - m[2] + p[2] - p[3];
    let f = m[3];
    let kk = m[1] - m[2] - p[5] + p[2];
    let dpe = m[0] - m[3] - p[3];
    let dpk = m[1] - m[3] - p[5];
    let dpf = m[2] - p[2];
    t13(
        [
            a,
            c,
            g,
            h,
            d,
            e,
            f,
            kk,
            dpe,
            dpk,
            dpf,
            dm[2] - dm[3],
            dm[0] - dm[2],
            dm[3],
            dm[1] - dm[2],
            dm[0] - dm[3],
            dm[1] - dm[3],
            dm[2],
        ],
        k,
    )
}

fn opposite(p: [Z; 6], m: [Z; 4], k: &Constants<f64>) -> Z {
    const CYCLES: [[usize; 6]; 4] = [
        [0, 1, 2, 3, 4, 5],
        [1, 2, 3, 0, 5, 4],
        [2, 3, 0, 1, 4, 5],
        [3, 0, 1, 2, 5, 4],
    ];
    let chart = |rotation: usize| {
        let cycle = CYCLES[rotation];
        (
            cycle.map(|i| p[i]),
            core::array::from_fn::<_, 4, _>(|i| m[cycle[i]]),
        )
    };
    let selected = (0..4)
        .find(|&rotation| {
            let (p, m) = chart(rotation);
            opposite_regular(p, m)
        })
        .unwrap_or(0);
    let (p, m) = chart(selected);
    opposite_chart(p, m, k)
}

fn zero_external(p: [Z; 6], m: [Z; 4], k: &Constants<f64>) -> Z {
    let dm = mass_directions(m);
    let a = p[2];
    let b = p[1];
    let c = p[5] - p[1] - p[2];
    let h = p[3] - p[4] - p[5] + p[1];
    let j = p[4] - p[0] - p[1];
    let d = m[2] - m[3] - p[2];
    let e = m[1] - m[2] - p[5] + p[2];
    let kk = m[0] - m[1] + p[5] - p[3];
    let f = m[3];
    let cph = p[3] - p[4] - p[2];
    let dpe = m[1] - m[3] - p[5];
    let epk = m[0] - m[2] + p[2] - p[3];
    let dek = m[0] - m[3] - p[3];
    let dpf = m[2] - p[2];
    tfun(
        [
            a,
            b,
            c,
            h,
            j,
            d,
            e,
            f,
            kk,
            dpe,
            dpf,
            dm[2] - dm[3],
            dm[1] - dm[2],
            dm[3],
            dm[0] - dm[1],
            dm[1] - dm[3],
            dm[2],
            // a + b + c, simplified as in the expression.
            p[5],
        ],
        k,
    ) - tfun(
        [
            a,
            b + j,
            cph,
            h,
            j,
            d,
            epk,
            f,
            kk,
            dek,
            dpf,
            dm[2] - dm[3],
            dm[0] - dm[2],
            dm[3],
            dm[0] - dm[1],
            dm[0] - dm[3],
            dm[2],
            p[3] - p[0],
        ],
        k,
    )
}

fn ordered(p: [Z; 6], m: [Z; 4], k: &Constants<f64>) -> Z {
    let dm = mass_directions(m);
    let a = p[2];
    let b = p[1];
    let g = p[0];
    let c = p[5] - p[1] - p[2];
    let h = p[3] - p[4] - p[5] + p[1];
    let j = p[4] - p[0] - p[1];
    let d = m[2] - m[3] - p[2];
    let e = m[1] - m[2] - p[5] + p[2];
    let kk = m[0] - m[1] + p[5] - p[3];
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
    let [x1, x2] = roots(g, j, b);
    let alpha = if negative(x1.re * x1.re - x2.re * x2.re) { x2 } else { x1 };
    let o1 = ONE - alpha;
    let j1 = j + r(2.0) * g * alpha;
    let e1 = e + kk * alpha;
    let de1 = dm[1] - dm[2] + (dm[0] - dm[1]) * alpha;
    -tfun(
        [
            abc,
            g,
            jph,
            c + r(2.0) * b + (h + j) * alpha,
            j1,
            dpe,
            kk,
            f,
            e1,
            dek,
            def,
            dm[1] - dm[3],
            dm[0] - dm[1],
            dm[3],
            de1,
            dm[0] - dm[3],
            dm[1],
            p[3],
        ],
        k,
    ) + o1
        * tfun(
            [
                a,
                bgj,
                cph,
                c + h * alpha,
                o1 * j1,
                d,
                epk,
                f,
                e1,
                dek,
                dpf,
                dm[2] - dm[3],
                dm[0] - dm[2],
                dm[3],
                de1,
                dm[0] - dm[3],
                dm[2],
                p[3],
            ],
            k,
        )
        + alpha
            * tfun(
                [
                    a,
                    b,
                    c,
                    c + h * alpha,
                    -j1 * alpha,
                    d,
                    e,
                    f,
                    e1,
                    dpe,
                    dpf,
                    dm[2] - dm[3],
                    dm[1] - dm[2],
                    dm[3],
                    de1,
                    dm[1] - dm[3],
                    dm[2],
                    p[5],
                ],
                k,
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

const PRM_A: [usize; 24] = [
    17, 7, 22, 18, 8, 21, 23, 13, 4, 24, 14, 3, 5, 19, 10, 6, 20, 9, 11, 1, 16, 12, 2, 15,
];
const PRM_B: [usize; 24] = [
    5, 4, 1, 6, 3, 2, 11, 10, 7, 12, 9, 8, 17, 16, 13, 18, 15, 14, 23, 22, 19, 24, 21, 20,
];

/// `sort_four`: the leaf of the expression's comparison tree (1-based).
fn sort_four(a: f64, b: f64, c: f64, d: f64) -> usize {
    const OUTPUTS: [usize; 24] = [
        10, 9, 8, 4, 11, 12, 7, 3, 24, 23, 22, 2, 14, 13, 18, 5, 15, 16, 17, 6, 19, 20, 21, 1,
    ];
    let leaf = |index: usize| OUTPUTS[index - 1];
    // le(x, y) = 1 - less(y, x)
    let le = |x: f64, y: f64| !negative(y - x);
    let choose = |x: f64, y: f64, yes: usize, no: usize| if le(x, y) { yes } else { no };
    let branch1 = || choose(c, d, leaf(1), choose(b, d, leaf(2), choose(a, d, leaf(3), leaf(4))));
    let branch2 = || choose(b, d, leaf(5), choose(c, d, leaf(6), choose(a, d, leaf(7), leaf(8))));
    let branch3 = || choose(b, d, leaf(9), choose(a, d, leaf(10), choose(c, d, leaf(11), leaf(12))));
    let branch4 = || choose(c, d, leaf(13), choose(a, d, leaf(14), choose(b, d, leaf(15), leaf(16))));
    let branch5 = || choose(a, d, leaf(17), choose(c, d, leaf(18), choose(b, d, leaf(19), leaf(20))));
    let branch6 = || choose(a, d, leaf(21), choose(b, d, leaf(22), choose(c, d, leaf(23), leaf(24))));
    if le(a, b) {
        if le(b, c) {
            branch1()
        } else if le(a, c) {
            branch2()
        } else {
            branch3()
        }
    } else if le(a, c) {
        branch4()
    } else if le(b, c) {
        branch5()
    } else {
        branch6()
    }
}

/// Re of the four Källén functions of the three-point sub-kinematics.
#[inline]
fn lambdas(p: [Z; 6]) -> [f64; 4] {
    [(0, 1, 4), (1, 2, 5), (2, 3, 4), (3, 0, 5)].map(|(a, b, c)| {
        let x = p[c] - p[a] - p[b];
        (x * x - r(4.0) * p[a] * p[b]).re
    })
}

fn general(p: [Z; 6], m: [Z; 4], k: &Constants<f64>) -> Z {
    let lambda = lambdas(p);
    let sorted = sort_four(lambda[0], lambda[1], lambda[2], lambda[3]);
    let all_positive = lambda.iter().all(|l| !negative(*l));
    let permutation = if all_positive {
        PRM_B[sorted - 1]
    } else {
        PRM_A[sorted - 1]
    };
    let row = ORD[permutation - 1];
    let selected_m = core::array::from_fn(|slot| m[row[slot] - 1]);
    let selected_p = core::array::from_fn(|slot| p[row[slot + 4] - 1]);
    ordered(selected_p, selected_m, k)
}

/// `__olo_boxc_box`: the finite contour box for one mass/momentum labeling.
pub(super) fn boxc(p: [Z; 6], m: [Z; 4], k: &Constants<f64>) -> Z {
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
    let mask = (0..4).fold(0, |s, i| s + if nz(p[i]) { 1 << (3 - i) } else { 0 });
    let order = PERM[mask];
    let pp = order.map(|i| p[i]);
    let mm = core::array::from_fn(|i| m[order[i]]);
    if mask == 15 {
        general(pp, mm, k)
    } else if matches!(CASE[mask], 0 | 1 | 5) {
        opposite(pp, mm, k)
    } else {
        zero_external(pp, mm, k)
    }
}

/// `supported_kinematics` of the contour representation.
pub(super) fn supported_kinematics(p: [Z; 6]) -> bool {
    if !nz(p[4]) || !nz(p[5]) {
        return false;
    }
    if p[..4].iter().all(|v| nz(*v)) {
        lambdas(p).iter().any(|l| negative(-*l))
    } else {
        true
    }
}

/// Each ported definition against the registered `__olo_boxc_*` expression,
/// evaluated by Symbolica in binary64 on random (often degenerate) arguments.
#[cfg(test)]
mod expression_tests {
    use super::*;
    use crate::definitions::FunctionMap;
    use symbolica::atom::{Atom, AtomCore};
    use symbolica::prelude::RealLike;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn uniform(&mut self) -> f64 {
            (self.next() >> 11) as f64 / (1u64 << 53) as f64
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
        fn real(&mut self) -> f64 {
            match self.below(5) {
                0 => 0.0,
                1 => (self.below(9) as f64 - 4.0) / 2.0,
                _ => 4.0 * self.uniform() - 2.0,
            }
        }
        fn complex(&mut self) -> Z {
            let re = self.real();
            let im = if self.below(2) == 0 { 0.0 } else { self.real() };
            Z { re, im }
        }
    }

    fn evaluator(name: &str, arity: usize) -> impl FnMut(&[Z]) -> Z {
        let vars: Vec<Atom> = (0..arity)
            .map(|i| symbolica::symbol!(format!("handwritten_boxc_arg_{i}")).to_atom())
            .collect();
        let call = symbolica::symbol!(format!("__olo_boxc_{name}")).call(&vars);
        let mut map = FunctionMap::new();
        crate::sheet_exact::register(&mut map);
        crate::box_complex::register(&mut map);
        let mut evaluator = Atom::evaluator_multiple(&[call.as_view()], &vars)
            .function_map(map.into())
            .direct_translation(true)
            .build()
            .unwrap()
            .map_coeff(&|v| Z::new(v.re.to_f64(), v.im.to_f64()));
        move |input: &[Z]| {
            let mut out = [ZERO];
            evaluator.evaluate(input, &mut out);
            out[0]
        }
    }

    /// The registered definition at 256 bits (rounded to binary64).
    fn precise(name: &str, x: &[Z]) -> Z {
        use symbolica::domains::float::Float;
        let vars: Vec<Atom> = (0..x.len())
            .map(|i| symbolica::symbol!(format!("handwritten_boxc_arg_{i}")).to_atom())
            .collect();
        let call = symbolica::symbol!(format!("__olo_boxc_{name}")).call(&vars);
        let mut map = FunctionMap::new();
        crate::sheet_exact::register(&mut map);
        crate::box_complex::register(&mut map);
        let mut evaluator = Atom::evaluator_multiple(&[call.as_view()], &vars)
            .function_map(map.into())
            .direct_translation(true)
            .build()
            .unwrap()
            .map_coeff_with_prec(
                &|v| symbolica::domains::float::Complex::new(v.re.to_multi_prec_float(256), v.im.to_multi_prec_float(256)),
                256,
            );
        let input: Vec<_> = x
            .iter()
            .map(|z| symbolica::domains::float::Complex::new(Float::with_val(256, z.re), Float::with_val(256, z.im)))
            .collect();
        let mut out = [symbolica::domains::float::Complex::new(Float::new(256), Float::new(256))];
        evaluator.evaluate(&input, &mut out);
        Z::new(out[0].re.to_f64(), out[0].im.to_f64())
    }

    fn cont(x: &[Z]) -> Cont {
        Cont { v: x[0], f: x[1], s: x[2] }
    }

    /// `allowed`: tolerated disagreements with the 256-bit expression, for a
    /// composition whose binary64 sheet tests can be decided by roundoff (the
    /// generated binary64 code shares them; tests/native_handwritten.rs checks
    /// the D0 accuracy against it).
    fn compare(name: &str, arity: usize, samples: usize, allowed: usize, port: impl Fn(&[Z]) -> Z, generate: &dyn Fn(&mut Rng) -> Vec<Z>) {
        let mut exact = evaluator(name, arity);
        let mut rng = Rng(0x9e3779b97f4a7c15 ^ arity as u64);
        let mut failures = vec![];
        let mut compared = 0;
        let mut rescued = 0;
        for _ in 0..samples {
            let x = generate(&mut rng);
            assert_eq!(x.len(), arity);
            let want = exact(&x);
            if !(want.re.is_finite() && want.im.is_finite()) {
                continue;
            }
            compared += 1;
            let got = port(&x);
            let abs = |z: Z| z.re.hypot(z.im);
            let error = abs(got - want) / abs(want).max(1e-300);
            if !(error < 1e-8) && abs(got - want) > 1e-10 {
                // Decide with the expression at 256 bits: near removable
                // singularities binary64 rounding can select another branch.
                let truth = precise(name, &x);
                if abs(got - truth) <= 1e-8 * abs(truth).max(1e-2) {
                    rescued += 1;
                } else {
                    failures.push(format!("{error:.2e}: {x:?}\n   port={got:?}\n   expr={want:?}\n   256-bit expr={truth:?}"));
                }
            }
        }
        eprintln!(
            "{name}: {compared} compared, {} failures; {rescued} binary64 expression values off the 256-bit expression where the port agrees with it",
            failures.len()
        );
        assert!(
            failures.len() <= allowed,
            "{name}: {} failures\n{}",
            failures.len(),
            failures[..failures.len().min(8)].join("\n")
        );
    }

    fn k() -> Constants<f64> {
        Constants::new(&0.0)
    }

    fn any(n: usize) -> impl Fn(&mut Rng) -> Vec<Z> {
        move |rng: &mut Rng| (0..n).map(|_| rng.complex()).collect()
    }

    /// Physical box arguments: real momenta, masses with nonpositive widths.
    fn kinematics(rng: &mut Rng) -> Vec<Z> {
        let mut x: Vec<Z> = (0..6).map(|_| r(rng.real() * 2.0)).collect();
        // The contour box is used with three or four massive propagators.
        let massless = if rng.below(3) == 0 { rng.below(4) } else { 4 };
        for slot in 0..4 {
            let re = match rng.below(3) {
                _ if slot == massless => 0.0,
                0 => 0.5,
                _ => 0.05 + rng.uniform(),
            };
            let im = if rng.below(3) == 0 { -0.05 * rng.uniform() } else { 0.0 };
            x.push(Z { re, im });
        }
        if rng.below(3) == 0 {
            let (i, j) = (rng.below(10) as usize, rng.below(10) as usize);
            if (i < 6) == (j < 6) {
                x[i] = x[j];
            }
        }
        x
    }

    fn split(x: &[Z]) -> ([Z; 6], [Z; 4]) {
        (core::array::from_fn(|i| x[i]), core::array::from_fn(|i| x[6 + i]))
    }

    #[test]
    fn contour_box_definitions_match_expressions() {
        crate::tests::large_stack(|| {
            let c = k();
            compare("root_1", 3, 3000, 0, |x| roots(x[0], x[1], x[2])[0], &any(3));
            compare("root_2", 3, 3000, 0, |x| roots(x[0], x[1], x[2])[1], &any(3));
            compare("r0", 6, 3000, 0, |x| r0(&cont(&x[0..3]), &cont(&x[3..6]), &c), &any(6));
            compare("r1", 10, 3000, 0, |x| r1(&cont(&x[0..3]), &cont(&x[3..6]), &cont(&x[6..9]), x[9], &c), &any(10));
            compare(
                "s3",
                13,
                3000,
                0,
                |x| s3(&cont(&x[0..3]), &cont(&x[3..6]), x[6], x[7], x[8], x[9], x[10], x[11], x[12], &c),
                &any(13),
            );
            compare(
                "planar",
                13,
                2000,
                0,
                |x| planar(&cont(&x[0..3]), &cont(&x[3..6]), nz(x[6]), nz(x[7]), x[8], x[9], x[10], x[11], nz(x[12]), &c),
                &|rng: &mut Rng| {
                    let mut x = any(13)(rng);
                    x[6] = r(rng.below(2) as f64);
                    x[7] = r(rng.below(2) as f64);
                    x[12] = r(rng.below(2) as f64);
                    x
                },
            );
            // The last argument is `a + c`, which the callers form exactly; the
            // first 14 keep the random stream of the 14-argument comparison.
            let stream = std::cell::RefCell::new(Rng(0x9e3779b97f4a7c15 ^ 14));
            compare("t1", 15, 1500, 0, |x| t1(core::array::from_fn(|i| x[i]), &c), &|_| {
                let mut x = any(14)(&mut stream.borrow_mut());
                x.push(x[0] + x[1]);
                x
            });
            compare("t13", 18, 1500, 0, |x| t13(core::array::from_fn(|i| x[i]), &c), &any(18));
            // Physical structure: real momentum combinations, masses with
            // nonpositive widths, imaginary infinitesimal directions.
            let _structured = |rng: &mut Rng| {
                let mut x: Vec<Z> = (0..5).map(|_| r(rng.real())).collect();
                for _ in 5..11 {
                    let mut z = r(rng.real());
                    if rng.below(3) == 0 {
                        z.im = -0.2 * rng.uniform();
                    }
                    x.push(z);
                }
                for _ in 11..17 {
                    x.push(Z { re: 0.0, im: rng.real() });
                }
                x
            };
            // tfun alone is not checked on random arguments: with a real
            // contour parameter its two candidate betas need not agree for
            // inconsistent random widths, and a tie |beta1| = |beta2| is broken
            // by rounding. The physical box compositions below are checked.
            // The charts are composed by `box` (each is only used inside its
            // own external-zero sector); compare on supported kinematics.
            compare(
                "box",
                10,
                3000,
                3,
                |x| {
                    let (p, m) = split(x);
                    boxc(p, m, &c)
                },
                &|rng: &mut Rng| loop {
                    let x = kinematics(rng);
                    if supported_kinematics(split(&x).0) {
                        break x;
                    }
                },
            );
        });
    }
}
