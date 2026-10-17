//! Binary64 D0 dispatch for the hand-written box sectors (`box_series`,
//! `select_four_series`, `box_finite_four_masses`, `finite_box_representation`).
use super::box_contour::{Z, boxc, dv, negative, nz, r, supported_kinematics};
use super::sheet::{Sheet, dd};
use super::{Constants, physical_sqrt, r_function};

/// `box_finite_four_masses` (OneLOop boxf4, Denner-Nierste-Scharf), the
/// representation used outside the contour box's kinematics.
fn four_masses_dns(p: [Z; 6], m: [Z; 4], k: &Constants<f64>) -> Z {
    let roots = m.map(|x| physical_sqrt(&x));
    let kk = |mi: usize, mj: usize, momentum: Z| dv(m[mi] + m[mj] - momentum, roots[mi] * roots[mj]);
    let k12 = kk(0, 1, p[0]);
    let k13 = kk(0, 2, p[4]);
    let k14 = kk(0, 3, p[3]);
    let k23 = kk(1, 2, p[1]);
    let k24 = kk(1, 3, p[5]);
    let k34 = kk(2, 3, p[2]);
    let (r12, _) = r_function(&k12, k);
    let (r13, d13) = r_function(&k13, k);
    let (r14, _) = r_function(&k14, k);
    let (r23, _) = r_function(&k23, k);
    let (r24, d24) = r_function(&k24, k);
    let (r34, _) = r_function(&k34, k);
    let a = dv(k34, r24) + r13 * k12 - dv(k14 * r13, r24) - k23;
    if !nz(a) {
        return Z::new(0.0, 0.0);
    }
    let b = d13 * d24 + k12 * k34 - k14 * k23;
    let c = dv(k12, r13) + r24 * k34 - dv(k14 * r24, r13) - k23;
    let discriminant = physical_sqrt(&(b * b - r(4.0) * a * c));
    let two_a = r(2.0) * a;
    let x1 = dv(b + discriminant, two_a);
    let x2 = dv(b - discriminant, two_a);
    let base = (k23 - r13 * k12 - r24 * k34 + r13 * r24 * k14).re * a.re;
    let sign = |x: f64| if negative(x) { -1 } else { 1 };
    let qx1 = Sheet::with_sign(x1, sign(base * x2.re));
    let qx2 = Sheet::with_sign(x2, sign(base * x1.re));
    let [q12, q13, q14, q23, q24, q34] = [r12, r13, r14, r23, r24, r34].map(Sheet::lower);
    let y1 = qx1.div(&q24);
    let y2 = qx2.div(&q24);
    let mut finite = dv(
        dd(&y1.mul(&q12), &y2.mul(&q12), k) * r12 + dv(dd(&y1.div(&q12), &y2.div(&q12), k), r12),
        r24,
    );
    let ratio = dv(r13, r24);
    let qr = Sheet::with_sign(ratio, -sign(r24.re));
    let z1 = qx1.mul(&qr);
    let z2 = qx2.mul(&qr);
    finite = finite
        - (dd(&z1.mul(&q23), &z2.mul(&q23), k) * r23 + dv(dd(&z1.div(&q23), &z2.div(&q23), k), r23))
            * ratio;
    let w1 = qx1.mul(&q13);
    let w2 = qx2.mul(&q13);
    finite = finite
        + (dd(&w1.mul(&q34), &w2.mul(&q34), k) * r34 + dv(dd(&w1.div(&q34), &w2.div(&q34), k), r34))
            * r13;
    finite = finite
        - dd(&qx1.mul(&q14), &qx2.mul(&q14), k) * r14
        - dv(dd(&qx1.div(&q14), &qx2.div(&q14), k), r14);
    -dv(finite, a * roots[0] * roots[1] * roots[2] * roots[3])
}

/// `box_series`: exchange the channel pair when s12 or s23 vanishes exactly.
#[inline]
fn hard_kinematics(p: [Z; 6], m: [Z; 4]) -> ([Z; 6], [Z; 4]) {
    if nz(p[4]) && nz(p[5]) {
        return (p, m);
    }
    let abs = p.map(|v| <f64 as crate::native::NativeFloat>::hypot(&v.re, &v.im));
    let min = |a: f64, b: f64| if negative(a - b) { a } else { b };
    let min13 = min(abs[0], abs[2]);
    let min24 = min(abs[1], abs[3]);
    let min56 = min(abs[4], abs[5]);
    if negative(min24 - min13) && negative(min56 - min13) {
        ([4, 1, 5, 3, 0, 2].map(|i| p[i]), [0, 2, 1, 3].map(|i| m[i]))
    } else if negative(min13 - min24) && negative(min56 - min24) {
        ([0, 5, 2, 4, 3, 1].map(|i| p[i]), [0, 1, 3, 2].map(|i| m[i]))
    } else {
        (p, m)
    }
}

/// `finite_box_representation`: OneLOop's contour representation for complex
/// masses and these timelike configurations, where it supports the kinematics.
#[inline]
fn use_contour(p: [Z; 6], m: [Z; 4]) -> bool {
    let complex_mass = m.iter().any(|x| x.im != 0.0);
    let all_external_nonnegative = p[..4].iter().all(|x| !negative(x.re));
    let channels_nonnegative = !negative(p[4].re) && !negative(p[5].re);
    (complex_mass || all_external_nonnegative || channels_nonnegative) && supported_kinematics(p)
}

/// The finite part of a covered D0 sector, or `None` for the generated code.
pub(super) fn d0(a: &[Z], k: &Constants<f64>) -> Option<Z> {
    let p: [Z; 6] = core::array::from_fn(|i| a[i]);
    let m: [Z; 4] = core::array::from_fn(|i| a[6 + i]);
    let (p, m) = hard_kinematics(p, m);
    let sector = (0..4).fold(0, |s, i| s + if nz(m[i]) { 8 >> i } else { 0 });
    match sector {
        15 => {
            // box_finite_four_masses
            let nonvacuum = p.iter().any(|x| nz(*x)) || m[1..].iter().any(|x| nz(*x - m[0]));
            if !nonvacuum {
                return Some(dv(r(1.0), r(6.0) * m[0] * m[0]));
            }
            Some(if use_contour(p, m) {
                boxc(p, m, k)
            } else {
                four_masses_dns(p, m, k)
            })
        }
        7 | 11 | 13 | 14 => {
            // box_three_masses_ir in the canonical labeling with m1 = 0.
            let shift = match sector {
                7 => 0,
                11 => 1,
                13 => 2,
                _ => 3,
            };
            const ORDER: [[usize; 6]; 4] = [
                [0, 1, 2, 3, 4, 5],
                [1, 2, 3, 0, 5, 4],
                [2, 3, 0, 1, 4, 5],
                [3, 0, 1, 2, 5, 4],
            ];
            let p = ORDER[shift].map(|i| p[i]);
            let m: [Z; 4] = core::array::from_fn(|i| m[(i + shift) % 4]);
            // The infrared sector (both massless-adjacent legs on shell) and
            // the non-contour finite formula use the generated code.
            let finite = nz(p[0] - m[1]) || nz(p[3] - m[3]);
            (finite && use_contour(p, m)).then(|| boxc(p, m, k))
        }
        _ => None,
    }
}
