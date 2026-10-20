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
use super::sheet::{Sheet, dd, dd_given};

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

/// `box_complex::below`: `2^bits x < y` for squared norms (binary64's half
/// precision for 52 and 56, as on every expression route).
#[inline(always)]
fn below(bits: i32, x: f64, y: f64) -> bool {
    negative(2f64.powi(bits) * x - y)
}

/// `box_complex::agree`: `|v - w| < 2^-26 |w|`.
#[inline(always)]
fn agree(v: Z, w: Z) -> bool {
    below(52, norm(v - w), norm(w))
}

/// `box_complex::close`: `|v - w| < 2^-bits |w|`.
#[inline(always)]
fn close(bits: i32, v: Z, w: Z) -> bool {
    below(2 * bits, norm(v - w), norm(w))
}

/// Finite real and imaginary parts.
#[inline(always)]
fn finite(z: Z) -> bool {
    z.re.is_finite() && z.im.is_finite()
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
    continued_root_with_slope(y, r(2.0) * k * y + l, k, dl, dn, dn2)
}

/// `continued_root` given the slope `2 k y + l` at the root.
#[inline]
fn continued_root_with_slope(y: Z, derivative: Z, k: Z, dl: Z, dn: Z, dn2: Z) -> Cont {
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

/// `kallen`: λ(p, mi, mj) = (p - (√mi + √mj)²)(p - (√mi - √mj)²), each factor
/// formed from the invariants (exactly zero at a threshold).
#[inline]
fn kallen(p: Z, mi: Z, mj: Z) -> Z {
    if p.im == 0.0 && mi.im == 0.0 && mj.im == 0.0 && !negative(mi.re) && !negative(mj.re) {
        // Real masses and invariant: the same operations in real arithmetic.
        let (p, mi, mj) = (p.re, mi.re, mj.re);
        let plus = (mi + mj) + 2.0 * (mi * mj).sqrt();
        let minus = if mj == 0.0 {
            mi
        } else if mi == 0.0 {
            mj
        } else {
            (mi - mj) * (mi - mj) / plus
        };
        return r((p - plus) * (p - minus));
    }
    let plus = (mi + mj) + r(2.0) * principal_sqrt(mi * mj);
    let minus = if !nz(mj) {
        mi
    } else if !nz(mi) {
        mj
    } else {
        let difference = mi - mj;
        difference * difference * inverse(&plus)
    };
    (p - plus) * (p - minus)
}

/// `edge_kallens`: [λ_a, λ_abc, λ_bc] of a `tfun` conic with corner values
/// `f`, `q10`, `q11`.
#[inline]
fn edge_kallens(a: Z, b: Z, abc: Z, f: Z, q10: Z, q11: Z) -> [Z; 3] {
    [kallen(a, q10, f), kallen(abc, q11, f), kallen(b, q11, q10)]
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
    let (q0, n0, q1, n1) = r0_quotients(y1, y2);
    let value = dv(lom(&q0, k), y1.v) + dv(lom(&q1, k), ONE - y1.v);
    let winding = n0 - n1;
    if winding != 0 {
        let pi_n = k.pi * winding as f64;
        value + dv(Z { re: 0.0, im: pi_n }, y1.v - y2.v)
    } else {
        value
    }
}

/// The unwound quotients `y2/y1` and `(1 - y2)/(1 - y1)` of `r0` and their
/// windings.
#[inline(always)]
fn r0_quotients(y1: &Cont, y2: &Cont) -> (Sheet<f64>, i32, Sheet<f64>, i32) {
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
    (q0, n0, q1, n1)
}

/// `r0_line`: `r0` given `delta = y1 - y2`, as `(log(q0) - log(q1) + i pi
/// (n0 - n1))/delta` with the logarithms from the offsets `-delta/y1` and
/// `delta/(1 - y1)`; `r0` itself at `delta = 0`.
fn r0_line(y1: &Cont, y2: &Cont, delta: Z, k: &Constants<f64>) -> Z {
    if !nz(delta) {
        return r0(y1, y2, k);
    }
    let (q0, n0, q1, n1) = r0_quotients(y1, y2);
    // -delta/y1 and delta/(1 - y1) from one reciprocal.
    let oy1 = ONE - y1.v;
    let scaled = delta * inverse(&(y1.v * oy1));
    let mut logarithms = q0.log_given(&-(scaled * oy1), k) - q1.log_given(&(scaled * y1.v), k);
    logarithms.im += k.pi * (n0 - n1) as f64;
    dv(logarithms, delta)
}

/// `far_root_series`: `r1` for a far root, fy log(-z) - sum_{n<=8} M_n/(n zⁿ).
fn far_root_series(zz: &Cont, yy1: &Cont, yy2: &Cont, fy: Z, k: &Constants<f64>) -> Z {
    let l2 = yy2.complement(1.0).sheet().log(k) - yy2.complement(0.0).sheet().log(k);
    let sum_y = yy1.v + yy2.v;
    let product_y = yy1.v * yy2.v;
    let inv = inverse(&zz.v);
    let mut m = [ZERO; 9];
    m[0] = fy;
    m[1] = yy1.v * fy + l2;
    for n in 2..=8 {
        m[n] = r(1.0 / (n - 1) as f64) + sum_y * m[n - 1] - product_y * m[n - 2];
    }
    let mut sum = m[8] * r(1.0 / 8.0);
    for n in (1..8).rev() {
        sum = m[n] * r(1.0 / n as f64) + inv * sum;
    }
    fy * zz.complement(0.0).sheet().log(k) - inv * sum
}

fn r1(zz: &Cont, yy1: &Cont, yy2: &Cont, fy: Z, k: &Constants<f64>) -> Z {
    r1_with(zz, yy1, yy2, fy, None, k)
}

/// `r1` (no line) or `r1_line`: where the nearer pole's `y - z` cancels by
/// more than 64 (`4096 |y - z|² < |z|²`), its offset is the line's product
/// `(y1 - z)(y2 - z)` (`LineProduct`) over the farther pole's.
#[inline(always)]
fn r1_with(zz: &Cont, yy1: &Cont, yy2: &Cont, fy: Z, line: Option<&LineProduct>, k: &Constants<f64>) -> Z {
    let (z, y1, y2) = (zz.v, yy1.v, yy2.v);
    let oz = ONE - z;
    let (mut yz1, mut yz2) = (y1 - z, y2 - z);
    let (mut hz1, mut hz2) = (norm(yz1), norm(yz2));
    let hzz = norm(z);
    let mut given = [false; 2];
    if let Some(line) = line {
        let nearer_first = negative(hz1 - hz2);
        let nearer = if nearer_first { hz1 } else { hz2 };
        if negative(4096.0 * nearer - hzz) {
            // An exactly real `y - z` stays real (its lip follows the
            // directions); a vanishing product keeps `y - z`.
            let over = |near: Z, far: Z| {
                line.over(far).map(|offset| if near.im != 0.0 { offset } else { r(offset.re) })
            };
            if nearer_first {
                if let Some(offset) = over(yz1, yz2) {
                    yz1 = offset;
                    hz1 = norm(yz1);
                    given[0] = true;
                }
            } else if let Some(offset) = over(yz2, yz1) {
                yz2 = offset;
                hz2 = norm(yz2);
                given[1] = true;
            }
        }
    }
    let q1 = Cont { v: yz1, ..yy1.minus(zz) }.sheet();
    let q2 = Cont { v: yz2, ..yy2.minus(zz) }.sheet();
    let offsets = [yz1, yz2];
    let r0 = |y: &Cont, z: &Cont, index: usize, k: &Constants<f64>| {
        if given[index] { r0_line(y, z, offsets[index], k) } else { r0(y, z, k) }
    };
    let qz = zz.complement(0.0).sheet();
    let qo = zz.complement(1.0).sheet();
    let h12 = norm(y1 - y2);
    let hoz = norm(oz);
    // `far_root`; below |z|² = 2^14 its bound is at least 2^-52.
    if hzz >= 16384.0 {
        let larger = |a: f64, b: f64| if negative(a - b) { b } else { a };
        let scale = larger(larger(1.0, norm(y1)), norm(y2));
        let ratio = scale / hzz;
        let squared = ratio * ratio;
        let bound = 16.0 * squared * squared;
        if (1.0 + bound) - 1.0 == 0.0 {
            return far_root_series(zz, yy1, yy2, fy, k);
        }
    }
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
        fy * q1.log(k) - r0(yy2, zz, 1, k) * lom(&q1.div(&q2), k)
    } else if !negative(hz1 - h12) && !negative(hz2 - hz1) {
        fy * q2.log(k) - r0(yy1, zz, 0, k) * lom(&q2.div(&q1), k)
    } else {
        let first = if nz(yz1) {
            yz1 * q1.log(k) * r0(yy1, zz, 0, k)
        } else {
            ZERO
        };
        let second = if nz(yz2) {
            yz2 * q2.log(k) * r0(yy2, zz, 1, k)
        } else {
            ZERO
        };
        dv(first - second, y1 - y2)
    };
    // The complements 1 - x of the dilogarithms' arguments from the poles
    // (`r1`): a pole next to an end of the interval puts one next to one.
    // They are formed only where the dilogarithm takes them (`dd_given`).
    let dilog_z = if zsmall {
        if nz(z) {
            dv(dd_given(&qz.div(&q1), || dv(y1, yz1), &qz.div(&q2), || dv(y2, yz2), k) * z, yz1 * yz2)
        } else {
            ZERO
        }
    } else {
        dv(dd_given(&q1.div(&qz), || dv(y1, z), &q2.div(&qz), || dv(y2, z), k), z)
    };
    let dilog_o = if !zsmall && osmall {
        if nz(oz) {
            dv(dd_given(&qo.div(&q1), || dv(y1 - ONE, yz1), &qo.div(&q2), || dv(y2 - ONE, yz2), k) * oz, yz1 * yz2)
        } else {
            ZERO
        }
    } else {
        dv(dd_given(&q1.div(&qo), || dv(ONE - y1, oz), &q2.div(&qo), || dv(ONE - y2, oz), k), oz)
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

/// `s3`: the y1, y2 poles are mapped by y -> (d + y)/e; a quadratic passes
/// its two-point subgraph `[lambda, p, near, far]` (`Edge`), and a linear
/// b x + c passes `far` = b + c.
fn s3(
    y1: &Cont, y2: &Cont, d: Z, e: Z, a: Z, b: Z, c: Z, db: Z, dc: Z, edge: [Z; 4],
    k: &Constants<f64>,
) -> Z {
    s3_with(y1, y2, d, e, a, b, c, db, dc, edge, None, k)
}

/// `s3` (no line) or, for `a = db = 0`, `s3_line`: the linear term with the
/// line's product of the mapped poles' offsets from its root.
#[inline(always)]
fn s3_with(
    y1: &Cont, y2: &Cont, d: Z, e: Z, a: Z, b: Z, c: Z, db: Z, dc: Z, [lambda, p, near, far]: [Z; 4],
    line: Option<&LineProduct>,
    k: &Constants<f64>,
) -> Z {
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
        // A massless end (far = 0 or near = 0) is a root.
        let z0 = if !nz(far) {
            ONE
        } else if !nz(near) {
            ZERO
        } else {
            -dv(b, r(2.0) * a)
        };
        let double = !nz(lambda);
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
            let rea = sgn(a.re);
            let (mut z1, mut z2) = if double {
                // An exact double root splits like sqrt(eps), one root up and
                // one down, dominating every first-order shift of y1 and y2.
                let size = |z: Z| z.re.hypot(z.im);
                let dominant = 2.0 * (1.0 + size(y1.f) + size(y1.s) + size(y2.f) + size(y2.s));
                let centre = z0;
                (
                    Cont { v: centre, f: Z { re: 0.0, im: dominant }, s: ZERO },
                    Cont { v: centre, f: Z { re: 0.0, im: -dominant }, s: ZERO },
                )
            } else if !nz(far) || !nz(near) {
                // A massless end is a root and the other lies slope/p from it,
                // with the slope from the invariants; both continue from the
                // end value's direction.
                let (end, slope, direction) = if !nz(far) { (ONE, near - p, db + dc) } else { (ZERO, p - far, dc) };
                let inv_slope = inverse(&slope);
                let end_first = direction * inv_slope;
                let near_first = -dv(db, a) - direction * inv_slope;
                (
                    Cont {
                        v: end,
                        f: end_first,
                        s: (a * end_first * end_first + db * end_first) * inv_slope,
                    },
                    Cont {
                        v: end + dv(slope, r(p.re)),
                        f: near_first,
                        s: dv(near_first * a * direction, slope * slope),
                    },
                )
            } else {
                // Roots with slopes 2 a z + b = +sqrt(lambda), -sqrt(lambda).
                let ra = r(a.re);
                let delta = principal_sqrt(lambda);
                let plus = -b + delta;
                let minus = -b - delta;
                let (root1, root2) = if negative(norm(plus) - norm(minus)) {
                    (dv(r(2.0) * c, minus), dv(minus, r(2.0) * ra))
                } else {
                    (dv(plus, r(2.0) * ra), dv(r(2.0) * c, plus))
                };
                (
                    continued_root_with_slope(root1, delta, a, db, dc, ZERO),
                    continued_root_with_slope(root2, -delta, a, db, dc, ZERO),
                )
            };
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
        let root = if nz(far) { -dv(c, b) } else { ONE };
        let mut z = continued_root(root, ZERO, b, db, dc, ZERO);
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
            + r1_with(&z, &y1, &y2, fy, line, k)
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

/// `s3` of a linear `b x + c`, with `far` = b + c.
#[inline]
fn s3_linear(y1: &Cont, y2: &Cont, d: Z, e: Z, b: Z, c: Z, dc: Z, far: Z, k: &Constants<f64>) -> Z {
    s3(y1, y2, d, e, ZERO, b, c, ZERO, dc, [ZERO, ZERO, c, far], k)
}

/// `s3_line`: a line edge of `tfun_selected` (its `[d, e, b, c, far]` in
/// the `LineProduct`), whose product is given only for a line nearly parallel
/// to the contour direction (`parallel`).
fn s3_line(y1: &Cont, y2: &Cont, dc: Z, line: &LineProduct, parallel: bool, k: &Constants<f64>) -> Z {
    let [d, e, b, c, far] = line.line;
    s3_with(y1, y2, d, e, ZERO, b, c, ZERO, dc, [ZERO, ZERO, c, far], parallel.then_some(line), k)
}

/// `line_product`: the product `(Y1 - z)(Y2 - z)` of the mapped poles'
/// offsets from the root `z` of a line edge `b x + c` with `[d, e, b, c,
/// far]`, `a1 (Q(u) - ratio B(u))/(K e²)` at the preimage `u = e z - d`,
/// formed by `r1_with` only where an offset cancels. With the scale `s = b`
/// (1 for `far = 0`), `U = s u` and `ratio = num/den`, the nearer offset is
/// `a1 (den (Q_b U² + Q_e U s + Q_f s²) - num s (b1 U + c1 s))/(K s² den e²)`
/// over the farther one, with one division. `a1` is the poles' own, or its
/// rounding level `2^-52 g` where it vanishes (`tfun_selected`); `conic`
/// holds `b, e, f` and then `a, c, d`.
struct LineProduct<'a> {
    line: [Z; 5],
    ratio: [Z; 2],
    conic: &'a [Z; 6],
    data: &'a BetaData,
    a1: Z,
}

impl LineProduct<'_> {
    /// The product over the farther pole's offset `far` (zero for `far = 0`),
    /// or `None` where the product vanishes or the line's intersection with
    /// the edge lies on the conic to within the square root of the rounding
    /// (`line_product`; the expression's `r1_line` then keeps `y - z`).
    fn over(&self, far_offset: Z) -> Option<Z> {
        let [d, e, b, c, far] = self.line;
        let [qb, qe, qf, ..] = *self.conic;
        let BetaData { b: b1, c: c1, k: kk, .. } = *self.data;
        let a1 = self.a1;
        let (scale, scaled_u) = if nz(far) { (b, -(e * c + d * b)) } else { (ONE, e - d) };
        // ratio = 0 for a vanishing denominator, as the expression's if_nonzero.
        let [num, den] = if nz(self.ratio[1]) { self.ratio } else { [ZERO, ONE] };
        let quadratic = (qb * scaled_u + qe * scale) * scaled_u + qf * scale * scale;
        let linear = b1 * scaled_u + c1 * scale;
        let value = den * quadratic - num * scale * linear;
        // The conic's value against its coefficients, with the value's scale
        // den s² (`line_product`).
        let conic_norm = self.conic.iter().fold(0.0, |sum, &q| sum + norm(q));
        if !negative(conic_norm * norm(den * scale * scale) - 4503599627370496.0 * norm(value)) {
            return None;
        }
        let numerator = a1 * value;
        if !nz(numerator) {
            return None;
        }
        let denominator = kk * scale * scale * den * e * e * far_offset;
        Some(if nz(far_offset) && nz(e) { dv(numerator, denominator) } else { ZERO })
    }
}

/// `t1`; `gh`, `gj`, `ghj` are the line's g + h, g + j, g + h + j.
fn t1(
    [a, c, g, h, d, e, f, j, dpe, dd_, de, df, dj, ddpe, ac, lambda_ac, far_ac, lambda_a, far_a, gh, gj, ghj]: [Z; 22],
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
    let sum = coefficient_sum(&[a, c], ac);
    dv(
        -s3(&y1, &y2, ZERO, ONE, sum, dpe, f, ddpe, df, [lambda_ac, ac, f, far_ac], k)
            + s3_linear(&y1, &y2, ZERO, ONE, gh, j, dj, ghj, k)
            - s3_linear(&y1, &y2, ZERO, ONE, g, j, dj, gj, k)
            + s3(&y1, &y2, ZERO, ONE, a, d, f, dd_, df, [lambda_a, a, f, far_a], k),
        kk,
    )
}

fn t13(
    [
        a, c, g, h, d, e, f, j, dpe, dpj, dpf, dd_, de, df, dj, ddpe, ddpj, ddpf, lambda_ae, lambda_gj,
        lambda_gh, lambda_ac, far_ae, far_gj, gh, ac,
    ]: [Z; 26],
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
        -s3(&y1, &y2, ZERO, ONE, a, e + c, dpf, de, ddpf, [lambda_ae, a, dpf, far_ae], k)
            + s3(&y1, &y2, ZERO, ONE, g, j + h, dpf, dj, ddpf, [lambda_gj, g, dpf, far_gj], k)
            - s3(&y1, &y2, ZERO, ONE, g + h, dpj, f, ddpj, df, [lambda_gh, gh, f, far_gj], k)
            + s3(&y1, &y2, ZERO, ONE, a + c, dpe, f, ddpe, df, [lambda_ac, ac, f, far_ae], k),
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

/// `beta_data`; `w` is beta - 1, and nearer 1 than 0 the coefficients are
/// expanded about beta = 1 (`dpe` = d + e, `ddpe` = dd + de).
#[inline(always)]
fn beta_data(
    beta: Z, w: Z, b: Z, c: Z, g: Z, h: Z, d: Z, e: Z, f: Z, j: Z, dpe: Z, dd_: Z, de: Z, df: Z,
    dj: Z, ddpe: Z, far_pole: bool,
) -> BetaData {
    let (a1, b1, c1, dc1) = if negative(norm(w) - norm(beta)) {
        ((g + h) + w * h, (c + r(2.0) * b) + r(2.0) * w * b, dpe + w * e, ddpe + w * de)
    } else {
        (g + beta * h, c + r(2.0) * beta * b, d + beta * e, dd_ + beta * de)
    };
    let kk = b * a1 - h * b1;
    let l = e * a1 - h * c1 - j * b1;
    let n = f * a1 - j * c1;
    // In the chart (`far_pole`), a K that rounds to zero puts its far pole
    // 2^52 (1 + |N/L|²) away (`box_complex::beta_data`).
    let kk = if nz(kk) || !far_pole {
        kk
    } else if nz(l) {
        -dv(l, r(4503599627370496.0 * (1.0 + norm(dv(n, l)))))
    } else {
        ZERO
    };
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

/// `beta_data` with `K` replaced by `kk` (`box_complex::tfun_chart`'s chart
/// with K at its rounding level).
#[cold]
fn beta_data_with_k(
    beta: Z, w: Z, b: Z, c: Z, g: Z, h: Z, d: Z, e: Z, f: Z, j: Z, dpe: Z, dd_: Z, de: Z, df: Z,
    dj: Z, ddpe: Z, kk: Z,
) -> BetaData {
    let (a1, b1, c1, dc1) = if negative(norm(w) - norm(beta)) {
        ((g + h) + w * h, (c + r(2.0) * b) + r(2.0) * w * b, dpe + w * e, ddpe + w * de)
    } else {
        (g + beta * h, c + r(2.0) * beta * b, d + beta * e, dd_ + beta * de)
    };
    let l = e * a1 - h * c1 - j * b1;
    let n = f * a1 - j * c1;
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

/// `tfun_selected`; `w` is beta - 1 and `line` the line's g + h, g + j and
/// g + h + j.
fn tfun_selected(
    a: Z, b: Z, c: Z, g: Z, h: Z, d: Z, e: Z, f: Z, j: Z, dpe: Z, dpf: Z, beta: Z, w: Z,
    data: &BetaData, dd_: Z, de: Z, df: Z, dj: Z, ddpe: Z, ddpf: Z, abc: Z,
    [lambda_a, lambda_abc, lambda_bc, q10, q11]: [Z; 5], [gh, gj, ghj]: [Z; 3], k: &Constants<f64>,
) -> Z {
    let [y1, y2] = &data.y;
    let p1 = inside(beta, y1);
    let p2 = inside(beta, y2);
    let sum = coefficient_sum(&[a, b, c], abc);
    // The line edges x = 1, y = x and y = 0 as the `[d, e, b, c, far]` of
    // their s3 terms, and h u + j over a1 at their roots' preimages u.
    let conic = [b, e, f, a, c, d];
    // The products only for a line nearly parallel to the contour direction
    // (a1 = g + beta h cancels by more than 64), with an a1 that rounds to
    // zero at its rounding level 2^-52 g (`box_complex::tfun_selected`).
    let parallel = negative(4096.0 * norm(data.a) - (norm(g) + norm(beta * h)));
    let a1 = if nz(data.a) { data.a } else { g * r(1.0 / 4503599627370496.0) };
    let line = |line: [Z; 5], ratio: [Z; 2]| LineProduct { line, ratio, conic: &conic, data, a1 };
    let integral = s3_line(y1, y2, dj, &line([beta, ONE, h, gj, ghj], [r(-1.0), ONE]), parallel, k)
        - s3_line(y1, y2, dj, &line([ZERO, -w, gh, j, ghj], [j, gh]), parallel, k)
        + s3_line(y1, y2, dj, &line([ZERO, -beta, g, j, gj], [j, g]), parallel, k)
        - s3(y1, y2, beta, ONE, b, c + e, a + dpf, de, ddpf, [lambda_bc, b, q10, q11], k)
        + s3(y1, y2, ZERO, -w, sum, dpe, f, ddpe, df, [lambda_abc, abc, f, q11], k)
        - s3(y1, y2, ZERO, -beta, a, d, f, dd_, df, [lambda_a, a, f, q10], k);
    let real = real_intersections_absent([a, b, c, g, h, d, e, f, j]);
    let residue = planar(y1, y2, p1, p2, data.a, data.b, data.c, data.dc, real, k);
    let residue = if negative(0.0 - beta.im) { -residue } else { residue };
    dv(integral + residue, data.k)
}

/// The rare charts of `tfun` (`box_complex::tfun_chart`): the selected
/// contour (`[beta, w]`, `data`) with a far pole where its K rounds to zero;
/// where `check`, the other contour's chart too, taken unless the two agree
/// to 2^-26 (or it is not finite). The chart used is tested where it is the
/// selected `flat` one without a regular chart that agrees with it (NaN
/// unless the chart with the line's g moved to `g (1 + 2^-52)`, and its sums
/// by the same exact amount, keeps 10 bits of it) and where its K rounds to
/// zero or it is the other chart (whose K is `b a1` of the flat selected
/// one): NaN unless the chart with K at its rounding level keeps 13 bits.
#[cold]
#[inline(never)]
fn tfun_checked(
    [a, b, c, g, h, d, e, f, j, dpe, dpf, dd_, de, df, dj, ddpe, ddpf, abc]: [Z; 18],
    [beta, w, other_beta, other_w]: [Z; 4],
    [data, other]: [&BetaData; 2],
    check: bool,
    flat: bool,
    edges: [Z; 5],
    sums: [Z; 3],
    k: &Constants<f64>,
) -> Z {
    let chart = |beta: Z, w: Z, data: &BetaData| {
        let far = (!nz(data.k)).then(|| beta_data(beta, w, b, c, g, h, d, e, f, j, dpe, dd_, de, df, dj, ddpe, true));
        tfun_selected(a, b, c, g, h, d, e, f, j, dpe, dpf, beta, w, far.as_ref().unwrap_or(data), dd_, de, df, dj, ddpe, ddpf, abc, edges, sums, k)
    };
    let tested = |beta: Z, w: Z, data: &BetaData, used: Z, flat: bool, zero: bool| {
        let keeps = |moved: Z, bits: i32| finite(moved) && close(bits, used, moved);
        if flat {
            let shifted = g * r(1.0 + f64::EPSILON);
            let delta = shifted - g;
            let moved = beta_data(beta, w, b, c, shifted, h, d, e, f, j, dpe, dd_, de, df, dj, ddpe, true);
            let moved = tfun_selected(a, b, c, shifted, h, d, e, f, j, dpe, dpf, beta, w, &moved, dd_, de, df, dj, ddpe, ddpf, abc, edges, sums.map(|sum| sum + delta), k);
            if !keeps(moved, 10) {
                return Z { re: f64::NAN, im: f64::NAN };
            }
        }
        if zero {
            let magnitude = |z: Z| z.re.abs() + z.im.abs();
            let rounding = r((magnitude(b * data.a) + magnitude(h * data.b)) / 4503599627370496.0);
            let probe = beta_data_with_k(beta, w, b, c, g, h, d, e, f, j, dpe, dd_, de, df, dj, ddpe, rounding);
            let moved = tfun_selected(a, b, c, g, h, d, e, f, j, dpe, dpf, beta, w, &probe, dd_, de, df, dj, ddpe, ddpf, abc, edges, sums, k);
            if !keeps(moved, 13) {
                return Z { re: f64::NAN, im: f64::NAN };
            }
        }
        used
    };
    let value = chart(beta, w, data);
    if check {
        let regular = chart(other_beta, other_w, other);
        if finite(regular) && !(finite(value) && agree(value, regular)) {
            return tested(other_beta, other_w, other, regular, false, true);
        }
        return if finite(regular) { value } else { tested(beta, w, data, value, true, !nz(data.k)) };
    }
    if flat || !nz(data.k) { tested(beta, w, data, value, flat, !nz(data.k)) } else { value }
}

/// `tfun`; `edges` holds `[λ_a, λ_abc, λ_bc, q10, q11]` (`edge_kallens` and
/// the conic's values at the corners (1, 0) and (1, 1)), `line_sums` the
/// line's g + h, g + j and g + h + j from the invariants.
fn tfun(
    [a, b, c, gin, hin, d, e, f, jin, dpe, dpf, dd_, de, df, djin, ddpe, ddpf, abc]: [Z; 18],
    edges: [Z; 5],
    line_sums: [Z; 3],
    k: &Constants<f64>,
) -> Z {
    let [lambda_a, lambda_abc, lambda_bc, q10, q11] = edges;
    let sj = isign(jin, isign(djin, -1));
    let s = r(-(sj as f64));
    let g = s * gin;
    let h = s * hin;
    let j = s * jin;
    let dj = s * djin;
    let [gh, gj, ghj] = line_sums.map(|v| s * v);
    let value = if nz(b) {
        if nz(a) {
            let [root1, root2] = roots(b, c, a);
            // beta - 1, the root nearer 1 from (beta1 - 1)(beta2 - 1) = abc/b.
            let [o1, o2] = [root1 - ONE, root2 - ONE];
            let w1 = if negative(norm(o1) - norm(o2)) { dv(abc, b * o2) } else { o1 };
            let w2 = if negative(norm(o2) - norm(o1)) { dv(abc, b * o1) } else { o2 };
            let swap = negative(norm(root2) - norm(root1));
            let ((beta1, w1), (beta2, w2)) = if swap { ((root2, w2), (root1, w1)) } else { ((root1, w1), (root2, w2)) };
            let data1 = beta_data(beta1, w1, b, c, g, h, d, e, f, j, dpe, dd_, de, df, dj, ddpe, false);
            let data2 = beta_data(beta2, w2, b, c, g, h, d, e, f, j, dpe, dd_, de, df, dj, ddpe, false);
            let min = |x: f64, y: f64| if negative(x - y) { x } else { y };
            let d1 = min(distance(beta1, &data1.y[0]), distance(beta1, &data1.y[1]));
            let d2 = min(distance(beta2, &data2.y[0]), distance(beta2, &data2.y[1]));
            // By default the first (`tfun`). Where the selected contour's a1
            // keeps fewer than 24 bits (`2^56 |a1|² < |g|² + |beta h|²`): if
            // the other's does not and the other has a pole within 2^7, the
            // other chart is computed too and taken unless the two agree to
            // 2^-26; otherwise the chart is tested against a1 moved by one
            // rounding unit (`tfun_chart`).
            let use_first = if d1 - 1.0 != 0.0 { negative(d2 - d1) } else { true };
            let flat = |data: &BetaData, beta: Z| below(56, norm(data.a), norm(g) + norm(beta * h));
            let near = |data: &BetaData| negative(norm(data.y[0].v) - 16384.0) || negative(norm(data.y[1].v) - 16384.0);
            let ((beta, data), (other_beta, other)) = if use_first { ((beta1, &data1), (beta2, &data2)) } else { ((beta2, &data2), (beta1, &data1)) };
            let flat_selected = flat(data, beta);
            let check = flat_selected && !flat(other, other_beta) && near(other);
            if flat_selected || !nz(data.k) {
                let [w, other_w] = if use_first { [w1, w2] } else { [w2, w1] };
                let coefficients = [a, b, c, g, h, d, e, f, j, dpe, dpf, dd_, de, df, dj, ddpe, ddpf, abc];
                return -(r(sj as f64) * tfun_checked(coefficients, [beta, w, other_beta, other_w], [data, other], check, flat_selected, edges, [gh, gj, ghj], k));
            }
            let w = if use_first { w1 } else { w2 };
            tfun_selected(a, b, c, g, h, d, e, f, j, dpe, dpf, beta, w, data, dd_, de, df, dj, ddpe, ddpf, abc, edges, [gh, gj, ghj], k)
        } else {
            // a = 0: b + c = abc, and the new a + c is b. The new conic takes
            // the values q11, q10 and f at the corners (0, 0), (1, 0), (1, 1),
            // and the new line (-g - h, g, g + h + j) the sums -h, j, g + j.
            let bc = coefficient_sum(&[b, c], abc);
            t1(
                [
                    bc,
                    -c,
                    -gh,
                    g,
                    -dpe - r(2.0) * bc,
                    d + c,
                    q11,
                    ghj,
                    -e - r(2.0) * b - c,
                    -ddpe,
                    dd_,
                    ddpe + df,
                    dj,
                    -de,
                    b,
                    lambda_bc,
                    q10,
                    lambda_abc,
                    f,
                    -h,
                    j,
                    gj,
                ],
                k,
            )
        }
    } else {
        // b = 0: a + c = abc.
        t1(
            [a, c, g, h, d, e, f, j, dpe, dd_, de, df, dj, ddpe, abc, lambda_abc, q11, lambda_a, q10, gh, gj, ghj],
            k,
        )
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
    let far_ae = m[0] - p[0] - p[2];
    let far_gk = m[1] - p[2];
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
            kallen(a, dpf, far_ae),
            kallen(g, dpf, far_gk),
            // g + h and a + c, simplified as in the expression.
            kallen(p[5] - p[2], f, far_gk),
            kallen(p[3] - p[0] - p[2], f, far_ae),
            far_ae,
            far_gk,
            p[5] - p[2],
            p[3] - p[0] - p[2],
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

/// `zero_external`: `zero_external_chart` reflected where `m2` is massless,
/// the line passes nearer the corner (1, 0) than (0, 0), and neither of the
/// reflected chart's conics has a double contour direction (a discriminant
/// `c² - 4 a b` within 2^-26 of its terms, `below`).
fn zero_external(p: [Z; 6], m: [Z; 4], k: &Constants<f64>) -> Z {
    if !nz(m[2]) && nz(m[3]) && nz(p[1]) {
        let at_10 = m[0] - m[1] + p[1] - p[4];
        let at_00 = m[0] - m[1] + p[5] - p[3];
        if negative(norm(at_10) - norm(at_00)) {
            let reflected = [p[0], p[3], p[2], p[1], p[5], p[4]];
            let distinct = |a: Z, b: Z, c: Z| {
                let [cc, ab] = [c * c, r(4.0) * a * b];
                !below(52, norm(cc - ab), norm(cc) + norm(ab))
            };
            let [p0, p1, p2, p3, p4, p5] = reflected;
            if distinct(p2, p1, p5 - p1 - p2) && distinct(p2, p4 - p0, p3 - p4 - p2) {
                return zero_external_chart(reflected, [m[1], m[0], m[3], m[2]], k);
            }
        }
    }
    zero_external_chart(p, m, k)
}

fn zero_external_chart(p: [Z; 6], m: [Z; 4], k: &Constants<f64>) -> Z {
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
    let [lambda_a, lambda_abc, lambda_bc] = edge_kallens(a, b, p[5], f, m[2], m[1]);
    // b + j, a + b + j + cph and the conic at (1, 1), simplified as in the
    // expression.
    let bj = p[4] - p[0];
    let corner = m[0] - p[0];
    let lambda_abj = kallen(p[3] - p[0], corner, f);
    let lambda_bj = kallen(bj, corner, m[2]);
    // The line (h, j, kk): h + j, h + kk and h + j + kk.
    let line = [p[3] - p[5] - p[0], m[0] - m[1] + p[1] - p[4], m[0] - m[1] - p[0]];
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
        [lambda_a, lambda_abc, lambda_bc, m[2], m[1]],
        line,
        k,
    ) - tfun(
        [
            a,
            bj,
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
        [lambda_a, lambda_abj, lambda_bj, m[2], corner],
        line,
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
    // The six Källén functions of the box (`edge_kallens` of the three tfun).
    let l01 = kallen(p[0], m[0], m[1]);
    let l12 = kallen(p[1], m[1], m[2]);
    let l23 = kallen(p[2], m[2], m[3]);
    let l03 = kallen(p[3], m[0], m[3]);
    let l02 = kallen(p[4], m[0], m[2]);
    let l13 = kallen(p[5], m[1], m[3]);
    // The lines and their coefficient sums (g + h, g + j, g + h + j).
    let sums = |[g, h, j]: [Z; 3]| [g + h, g + j, (g + h) + j];
    let line1 = [c + r(2.0) * b + (h + j) * alpha, j1, e1];
    let line2 = [c + h * alpha, o1 * j1, e1];
    let line3 = [c + h * alpha, -j1 * alpha, e1];
    -tfun(
        [
            abc,
            g,
            jph,
            line1[0],
            line1[1],
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
        [l13, l03, l01, m[1], m[0]],
        sums(line1),
        k,
    ) + o1
        * tfun(
            [
                a,
                bgj,
                cph,
                line2[0],
                line2[1],
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
            [l23, l03, l02, m[2], m[0]],
            sums(line2),
            k,
        )
        + alpha
            * tfun(
                [
                    a,
                    b,
                    c,
                    line3[0],
                    line3[1],
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
                [l23, l13, l12, m[2], m[1]],
                sums(line3),
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

    /// `r1` of a far root (|z| = 1e2 to 1e17, real or complex, poles in the
    /// strip around the interval): the hand-written port and the binary64
    /// expression against the expression at 256 bits, within 1e-12 (below the
    /// series threshold, |z| ~ 140, the general form keeps about eps |z|). Its
    /// dilogarithm terms, taken next to the cut on different sheets, cancel by
    /// eps |z| (1e-6 at |z| = 1e9, O(1) at 1e15) unless `far_root_series` is
    /// used.
    #[test]
    fn far_root_r1_matches_high_precision() {
        crate::tests::large_stack(|| {
            let c = k();
            let mut exact = evaluator("r1", 10);
            let mut rng = Rng(0x5eed);
            let mut failures = vec![];
            for sample in 0..32 {
                let y = |rng: &mut Rng| Cont {
                    v: Z { re: 2.0 * rng.uniform() - 0.5, im: rng.uniform() - 0.5 },
                    f: Z { re: rng.uniform() - 0.5, im: rng.uniform() - 0.5 },
                    s: Z { re: rng.uniform() - 0.5, im: rng.uniform() - 0.5 },
                };
                let (y1, y2) = (y(&mut rng), y(&mut rng));
                let exponent = 2 + (sample % 16) as i32;
                let sign = if rng.below(2) == 0 { 1.0 } else { -1.0 };
                let big = sign * 10f64.powi(exponent) * (1.0 + rng.uniform());
                let im = if sample % 3 == 0 { big * (rng.uniform() - 0.5) } else { 0.0 };
                let zz = Cont { v: Z { re: big, im }, f: Z { re: 0.0, im: -big }, s: ZERO };
                let fy = r0(&y1, &y2, &c);
                let args = [zz.v, zz.f, zz.s, y1.v, y1.f, y1.s, y2.v, y2.f, y2.s, fy];
                let truth = precise("r1", &args);
                let abs = |z: Z| z.re.hypot(z.im);
                for (route, value) in [("port", r1(&zz, &y1, &y2, fy, &c)), ("expression", exact(&args))] {
                    let error = abs(value - truth) / abs(truth);
                    if !(error <= 1e-12) {
                        failures.push(format!("{route} |z| = {:.1e}: {error:.1e} ({value:?} vs {truth:?})", abs(zz.v)));
                    }
                }
            }
            assert!(failures.is_empty(), "{failures:#?}");
        });
    }

    #[test]
    fn contour_box_definitions_match_expressions() {
        crate::tests::large_stack(|| {
            let c = k();
            compare("root_1", 3, 3000, 0, |x| roots(x[0], x[1], x[2])[0], &any(3));
            compare("root_2", 3, 3000, 0, |x| roots(x[0], x[1], x[2])[1], &any(3));
            compare("r0", 6, 3000, 0, |x| r0(&cont(&x[0..3]), &cont(&x[3..6]), &c), &any(6));
            compare("r1", 10, 3000, 0, |x| r1(&cont(&x[0..3]), &cont(&x[3..6]), &cont(&x[6..9]), x[9], &c), &any(10));
            // `r0` given `delta = y1 - y2`.
            compare("r0_line", 7, 3000, 0, |x| r0_line(&cont(&x[0..3]), &cont(&x[3..6]), x[6], &c), &|rng: &mut Rng| {
                let mut x = any(6)(rng);
                x.push(x[0] - x[3]);
                x
            });
            // The line variants take the product (y1 - z)(y2 - z) of the
            // poles' offsets; a third of the samples put a pole within 1e-6
            // of the root, where the product replaces its offset.
            let conic = [ZERO, ZERO, ONE, ZERO, ZERO, ZERO];
            let data = |product: Z| BetaData { a: product, b: ZERO, c: ZERO, dc: ZERO, k: ONE, y: [Cont { v: ZERO, f: ZERO, s: ZERO }; 2] };
            fn given<'a>(data: &'a BetaData, conic: &'a [Z; 6]) -> LineProduct<'a> {
                LineProduct { line: [ZERO, ONE, ONE, ZERO, ONE], ratio: [ZERO, ONE], conic, data, a1: data.a }
            }
            let near = |rng: &mut Rng, root: Z| root + Z { re: 1e-6 * rng.real(), im: 1e-6 * rng.real() };
            compare(
                "r1_line",
                11,
                3000,
                0,
                |x| r1_with(&cont(&x[0..3]), &cont(&x[3..6]), &cont(&x[6..9]), x[9], Some(&given(&data(x[10]), &conic)), &c),
                &|rng: &mut Rng| {
                    let mut x = any(10)(rng);
                    if rng.below(3) == 0 {
                        x[3] = near(rng, x[0]);
                    }
                    x.push((x[3] - x[0]) * (x[6] - x[0]));
                    x
                },
            );
            compare(
                "s3_line",
                13,
                3000,
                0,
                |x| {
                    // `over` maps the product by 1/e², the expression takes it mapped.
                    let data = data(x[12] * x[7] * x[7]);
                    let line = LineProduct { line: [x[6], x[7], x[8], x[9], x[11]], ..given(&data, &conic) };
                    s3_line(&cont(&x[0..3]), &cont(&x[3..6]), x[10], &line, true, &c)
                },
                &|rng: &mut Rng| {
                    let mut x = any(11)(rng);
                    let (d, e, b) = (x[6], x[7], x[8]);
                    // far = b + c, or a massless end (c = -b, root 1).
                    let massless = rng.below(4) == 0;
                    if massless {
                        x[9] = -b;
                    }
                    let far = if massless { ZERO } else { b + x[9] };
                    let root = if nz(far) { -dv(x[9], b) } else { ONE };
                    if rng.below(3) == 0 {
                        // A pole whose image (d + y)/e lies next to the root.
                        x[0] = near(rng, e * root - d);
                    }
                    let mapped = |y: Z| dv(d + y, e) - root;
                    x.extend([far, mapped(x[0]) * mapped(x[3])]);
                    x
                },
            );
            compare("kallen", 3, 3000, 0, |x| kallen(x[0], x[1], x[2]), &any(3));
            // A quadratic's two-point subgraph [lambda, p, near, far]; often a
            // massless end. `s3` takes the roots of re(a) x² + b x + c (a is real
            // for every edge of the box), so lambda is their discriminant.
            compare(
                "s3",
                17,
                3000,
                0,
                |x| {
                    let edge = [x[13], x[14], x[15], x[16]];
                    s3(&cont(&x[0..3]), &cont(&x[3..6]), x[6], x[7], x[8], x[9], x[10], x[11], x[12], edge, &c)
                },
                &|rng: &mut Rng| {
                    let mut x = any(13)(rng);
                    let (a, b) = (x[8], x[9]);
                    let end = rng.below(4);
                    if end == 0 {
                        x[10] = -(a + b);
                    } else if end == 1 {
                        x[10] = ZERO;
                    }
                    let c = x[10];
                    x.extend([b * b - r(4.0) * r(a.re) * c, a, c, if end == 0 { ZERO } else { a + b + c }]);
                    x
                },
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
            // `a + c` and the discriminants and values at x = 1 of the edge
            // quadratics, which the callers form exactly; the first 14 keep the
            // random stream of the 14-argument comparison.
            let stream = std::cell::RefCell::new(Rng(0x9e3779b97f4a7c15 ^ 14));
            compare("t1", 22, 1500, 0, |x| t1(core::array::from_fn(|i| x[i]), &c), &|_| {
                let mut x = any(14)(&mut stream.borrow_mut());
                let [a, c, g, h, d, f, j, dpe] = [x[0], x[1], x[2], x[3], x[4], x[6], x[7], x[8]];
                let ac = a + c;
                x.extend([ac, dpe * dpe - r(4.0) * r(ac.re) * f, ac + dpe + f, d * d - r(4.0) * r(a.re) * f, a + d + f]);
                x.extend([g + h, g + j, g + h + j]);
                x
            });
            compare("t13", 26, 1500, 0, |x| t13(core::array::from_fn(|i| x[i]), &c), &|rng: &mut Rng| {
                let mut x = any(18)(rng);
                let [a, c, g, h, e, f, j, dpf] = [x[0], x[1], x[2], x[3], x[5], x[6], x[7], x[10]];
                let far_ae = a + e + c + dpf;
                let far_gj = g + j + h + dpf;
                // dpe and dpj consistent with the conic's corner values.
                let dpe = far_ae - a - c - f;
                let dpj = far_gj - g - h - f;
                x[8] = dpe;
                x[9] = dpj;
                x.extend([
                    (e + c) * (e + c) - r(4.0) * r(a.re) * dpf,
                    (j + h) * (j + h) - r(4.0) * r(g.re) * dpf,
                    dpj * dpj - r(4.0) * r((g + h).re) * f,
                    dpe * dpe - r(4.0) * r((a + c).re) * f,
                    far_ae,
                    far_gj,
                    g + h,
                    a + c,
                ]);
                x
            });
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
