//! C0 with three nonzero internal masses (`register_triangles` sector 7).
use super::super::primitives::{C, NativeFloat, mul};
use super::sheet::{Sheet, dd};
use super::{Constants, add, div, is_zero, negative, norm, physical_sqrt, r_function, real, scale, sub};

/// The finite part of sector 7, or `None` where the generated code is used.
pub(crate) fn three_masses<T: NativeFloat>(
    p: [&C<T>; 3],
    m: [&C<T>; 3],
    c: &Constants<T>,
) -> Option<C<T>> {
    if p.iter().all(|z| is_zero(z)) {
        return None;
    }
    let lambda = kallen(p, c);
    if negative(&lambda) {
        return None;
    }
    Some(select(p, m, c))
}

/// Re[(p3-p1-p2)² - 4 p1 p2] for real momenta.
#[inline]
fn kallen<T: NativeFloat>(p: [&C<T>; 3], c: &Constants<T>) -> T {
    let d = p[2].re.clone() - &p[0].re - &p[1].re;
    d.clone() * &d - c.four.clone() * &p[0].re * &p[1].re
}

#[inline]
fn rotate<'a, T: NativeFloat>(
    p: [&'a C<T>; 3],
    m: [&'a C<T>; 3],
) -> ([&'a C<T>; 3], [&'a C<T>; 3]) {
    ([p[2], p[0], p[1]], [m[2], m[0], m[1]])
}

/// The zero test of the leading coefficient's resultant (exact polynomial).
#[inline]
fn resultant<T: NativeFloat>(p: [&C<T>; 3], m: [&C<T>; 3]) -> C<T> {
    let u = sub(m[2], p[1]);
    let v = sub(m[0], p[0]);
    let w = sub(&add(m[0], m[2]), p[2]);
    let first = mul(&mul(m[0], &u), &u);
    let second = mul(&mul(&w, &u), &v);
    let third = mul(&mul(m[2], &v), &v);
    add(&sub(&first, &second), &third)
}

fn select<T: NativeFloat>(p: [&C<T>; 3], m: [&C<T>; 3], c: &Constants<T>) -> C<T> {
    // trif3's relabeling when the middle propagator has the largest width
    // (always on the real axis).
    let middle_largest =
        !negative(&(m[0].im.clone() - &m[1].im)) && !negative(&(m[2].im.clone() - &m[1].im));
    let (p, m) = if middle_largest { rotate(p, m) } else { (p, m) };
    // first(p, m)
    if !is_zero(&resultant(p, m)) {
        return ordered(p, m, c);
    }
    // second(R(p, m))
    let (p, m) = rotate(p, m);
    if !is_zero(&resultant(p, m)) {
        return ordered(p, m, c);
    }
    let (p, m) = rotate(p, m);
    ordered(p, m, c)
}

/// `triangle_finite_three_masses_ordered` (OneLOop trif3, DNS representation).
fn ordered<T: NativeFloat>(p: [&C<T>; 3], m: [&C<T>; 3], c: &Constants<T>) -> C<T> {
    let sqrt_1 = physical_sqrt(m[0]);
    let sqrt_2 = physical_sqrt(m[1]);
    let sqrt_3 = physical_sqrt(m[2]);
    let k12 = div(&sub(&add(m[0], m[1]), p[0]), &mul(&sqrt_1, &sqrt_2));
    let k13 = div(&sub(&add(m[0], m[2]), p[2]), &mul(&sqrt_1, &sqrt_3));
    let k23 = div(&sub(&add(m[1], m[2]), p[1]), &mul(&sqrt_2, &sqrt_3));
    let (r12, _) = r_function(&k12, c);
    let (r13, d13) = r_function(&k13, c);
    let (r23, _) = r_function(&k23, c);
    let a = add(
        &sub(&div(&sqrt_2, &sqrt_3), &k23),
        &mul(&r13, &sub(&k12, &div(&sqrt_2, &sqrt_1))),
    );
    if is_zero(&a) {
        return real(c.one.zero());
    }
    let b = sub(
        &add(&div(&d13, &sqrt_2), &div(&k12, &sqrt_3)),
        &div(&k23, &sqrt_1),
    );
    let lambda = kallen(p, c);
    let masses = mul(&mul(m[0], m[1]), m[2]);
    let discriminant = physical_sqrt(&div(&real(lambda), &masses));
    let two_a = scale(&a, &c.two);
    let cc = div(
        &sub(&div(&sqrt_1, &sqrt_3), &super::super::primitives::inverse(&r13)),
        &mul(&sqrt_1, &sqrt_2),
    );
    let (root_1, root_2) = if is_zero(&cc) {
        (div(&b, &a), real(c.one.zero()))
    } else {
        // (b ± d)/(2a). Where one sum cancels, its root follows from the
        // product relation x1*x2 = c/a instead (as OneLOop's solabc): identical
        // in exact arithmetic, but a root cancelled to an exact zero must not
        // reach the logarithms. Without cancellation the direct formula keeps
        // exactly coincident (Gram-zero) roots coincident.
        let plus = add(&b, &discriminant);
        let minus = sub(&b, &discriminant);
        let (np, nm) = (norm(&plus), norm(&minus));
        let sixteenth = c.half.clone() * &c.half * &c.half * &c.half;
        if nm < np.clone() * &sixteenth {
            let root_1 = div(&plus, &two_a);
            (root_1.clone(), div(&div(&cc, &a), &root_1))
        } else if np < nm.clone() * &sixteenth {
            let root_2 = div(&minus, &two_a);
            (div(&div(&cc, &a), &root_2), root_2)
        } else {
            (div(&plus, &two_a), div(&minus, &two_a))
        }
    };
    let debug = std::env::var("HW_DEBUG").is_ok();
    if debug { eprintln!("ordered p={p:?} m={m:?}\n k12={k12:?} k13={k13:?} k23={k23:?}\n r12={r12:?} r13={r13:?} r23={r23:?}\n a={a:?} b={b:?} cc={cc:?} disc={discriminant:?}\n roots {root_1:?} {root_2:?}"); }
    let qx1 = Sheet::upper(root_1);
    let qx2 = Sheet::upper(root_2.clone());
    let q12 = Sheet::lower(r12.clone());
    let q13 = Sheet::lower(r13.clone());
    let q23 = Sheet::lower(r23.clone());
    let qm1 = Sheet::lower(sqrt_1.clone());
    let qm2 = Sheet::lower(sqrt_2.clone());
    let qm3 = Sheet::lower(sqrt_3.clone());

    let z1 = qx1.mul(&qm2);
    let z2 = qx2.mul(&qm2);
    if debug { eprintln!(" dd1 {:?} {:?}", dd(&z1.mul(&q12), &z2.mul(&q12), c), dd(&z1.div(&q12), &z2.div(&q12), c)); }
    let mut finite = mul(
        &add(
            &mul(&dd(&z1.mul(&q12), &z2.mul(&q12), c), &r12),
            &div(&dd(&z1.div(&q12), &z2.div(&q12), c), &r12),
        ),
        &sqrt_2,
    );
    let temp = q13.mul(&qm2);
    let t1 = qx1.mul(&temp);
    let t2 = qx2.mul(&temp);
    finite = sub(
        &finite,
        &mul(
            &mul(
                &add(
                    &mul(&dd(&t1.mul(&q23), &t2.mul(&q23), c), &r23),
                    &div(&dd(&t1.div(&q23), &t2.div(&q23), c), &r23),
                ),
                &r13,
            ),
            &sqrt_2,
        ),
    );
    if debug { eprintln!(" after2 {finite:?} t1={t1:?} t2={t2:?} q23={q23:?}"); }
    let y1 = qx1.mul(&q13);
    let y2 = qx2.mul(&q13);
    finite = add(
        &finite,
        &sub(
            &mul(&mul(&dd(&y1.mul(&qm3), &y2.mul(&qm3), c), &r13), &sqrt_3),
            &mul(&dd(&qx1.mul(&qm1), &qx2.mul(&qm1), c), &sqrt_1),
        ),
    );
    if debug { eprintln!(" after3 {finite:?}"); }
    if !is_zero(&root_2) {
        let log_m3 = qm3.log(c);
        let log_m1 = qm1.log(c);
        let first = mul(
            &y1.div(&y2).log_over_one_minus(c),
            &sub(&add(&y1.log(c), &y2.log(c)), &scale(&log_m3, &c.two)),
        );
        let second = mul(
            &qx1.div(&qx2).log_over_one_minus(c),
            &sub(&add(&qx1.log(c), &qx2.log(c)), &scale(&log_m1, &c.two)),
        );
        finite = add(&finite, &div(&sub(&first, &second), &scale(&root_2, &c.two)));
    }
    div(
        &finite,
        &mul(&mul(&mul(&a, &sqrt_1), &sqrt_2), &sqrt_3),
    )
}
