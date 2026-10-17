//! The 't Hooft--Veltman representation for a negative Källén discriminant.
//! All continuation terms are native Symbolica expressions.
use super::*;
use sheet_exact::{SheetAtom as Q, i, im, negative, sign_nonnegative};

/// 't Hooft--Veltman eta function from the imaginary-part signs of its
/// arguments, with a vanishing imaginary part resolved by the caller.
fn eta(sa: &Atom, sb: &Atom, sc: &Atom) -> Atom {
    i() * Symbol::PI.to_atom() * sc * (1 + sa * sb) * (1 - sa * sc) / 2
}

/// Sign of the imaginary part of `z`, or `real_axis` when it vanishes.
fn im_sign(z: &Atom, real_axis: &Atom) -> Atom {
    let imaginary = im(z);
    if_nonzero_else(&imaginary, sign_nonnegative(&imaginary), real_axis.clone())
}

/// int_0^1 [log(a y² + b y + c) - log(a y0² + b y0 + c)]/(y - y0) dy with
/// b = m1 - m2 - a and c = m2, where `z = (a y0² + b y0 + c)/a`.
///
/// Real quantities take their lips from the Feynman prescription, masses
/// m - i delta with a common delta: c and z move by -i delta and -i delta/a,
/// while the roots y± = (-b ± d)/(2a) move by ±i delta/d. Distinct real
/// roots therefore separate upwards and downwards, as do coincident roots at
/// a threshold, whose shift is ±sqrt(i delta/a).
fn s3(a: &Atom, m1: &Atom, m2: &Atom, y: &Atom, z: &Atom) -> Atom {
    let b = m1 - m2 - a;
    let disc = (&b * &b - Atom::num(4) * a * m2).sqrt();
    let roots = [
        (-&b + &disc) / (Atom::num(2) * a),
        (-&b - &disc) / (Atom::num(2) * a),
    ];
    // Imaginary-part signs of -y± (and of 1 - y±, y0 - y±) on the real axis.
    let lips = [Atom::num(-1), Atom::num(1)];
    let term = |r: &Atom, lip: &Atom| {
        let qy = Q::with_sign(y - r, lip.clone());
        (&Q::with_sign(-r, lip.clone()) / &qy).dilog()
            - (&Q::with_sign(1 - r, lip.clone()) / &qy).dilog()
    };
    let below_a = -sign_nonnegative(a);
    let z_sign = im_sign(z, &below_a);
    let correction = eta(
        &im_sign(&(-&roots[0]), &lips[0]),
        &im_sign(&(-&roots[1]), &lips[1]),
        &im_sign(&(m2 / a), &below_a),
    ) - eta(
        &im_sign(&(y - &roots[0]), &lips[0]),
        &im_sign(&(y - &roots[1]), &lips[1]),
        &z_sign,
    ) - Atom::num(2) * i() * Symbol::PI.to_atom() * negative(a) * (1 - &z_sign)
        / 2;
    // The sheet correction is locally constant. Its inactive logarithm may
    // have a removable projective singularity (for example y=0 on the
    // Euclidean sheet), so preserve the same lazy branch used by the other
    // sheet primitives rather than evaluating zero times a nonfinite value.
    term(&roots[0], &lips[0])
        + term(&roots[1], &lips[1])
        + if_nonzero(&correction, &correction * Q::upper((y - 1) / y).log())
}

fn s3_call(a: &Atom, m1: &Atom, m2: &Atom, t1: Atom, t2: Atom, lambda: &Atom) -> Atom {
    // The Källén function is negative here and OneLOop takes its square
    // root on the upper lip; spell this out instead of trusting the IEEE
    // sign of a zero imaginary part.
    let lambda_root = sheet_exact::i() * (-lambda).sqrt();
    let alpha = (-&t2 + &lambda_root) / (Atom::num(2) * a);
    let b = m1 - m2 - a;
    let y = -(&t1 + &b * alpha) / &lambda_root;
    // a y² + b y + m2 at this y equals a u²/lambda + m2 - b²/(4a), with
    // u = t1 - b t2/(2a). Evaluating that real form keeps rounding noise out
    // of the sign of Im z, which selects the continuation below.
    let u = t1 - &b * t2 / (Atom::num(2) * a);
    let z = &u * &u / lambda + m2 / a - &b * &b / (Atom::num(4) * a * a);
    symbolica::symbol!("__olo_hv_s3").call((a, m1, m2, &y, &z))
}

fn body(p: [&Atom; 3], m: [&Atom; 3]) -> Atom {
    let lambda = (p[2] - p[0] - p[1]).pow(2) - Atom::num(4) * p[0] * p[1];
    let root = sheet_exact::i() * (-&lambda).sqrt();
    -(s3_call(
        p[0],
        m[0],
        m[1],
        m[1] - m[2] + p[1],
        p[2] - p[0] - p[1],
        &lambda,
    ) - s3_call(
        p[2],
        m[0],
        m[2],
        -(m[0] - m[1]) + p[2] - p[1],
        p[1] - p[0] - p[2],
        &lambda,
    ) + s3_call(
        p[1],
        m[1],
        m[2],
        -(m[0] - m[1]) + p[2] - p[1],
        p[0] + p[1] - p[2],
        &lambda,
    )) / root
}

pub(crate) fn register(map: &mut FunctionMap) {
    let symbols = [
        symbolica::symbol!("__olo_hv_a"),
        symbolica::symbol!("__olo_hv_b"),
        symbolica::symbol!("__olo_hv_c"),
        symbolica::symbol!("__olo_hv_d"),
        symbolica::symbol!("__olo_hv_e"),
        symbolica::symbol!("__olo_hv_f"),
    ];
    let x = symbols.map(Symbol::to_atom);
    map.add_function_with_options(
        symbolica::symbol!("__olo_hv_s3"),
        symbols[..5].to_vec(),
        s3(&x[0], &x[1], &x[2], &x[3], &x[4]),
        native_function_options(),
    )
    .unwrap();
    map.add_function_with_options(
        symbolica::symbol!("__olo_triangle_hv"),
        symbols.to_vec(),
        body([&x[0], &x[1], &x[2]], [&x[3], &x[4], &x[5]]),
        native_function_options(),
    )
    .unwrap();
}

pub(crate) fn expression(p: [&Atom; 3], m: [&Atom; 3]) -> Atom {
    let mut p = p.map(Clone::clone);
    let mut m = m.map(Clone::clone);
    // Sorting an edge also exchanges the opposite propagator labels.
    for (a, b, c, d) in [(0, 1, 0, 2), (1, 2, 0, 1), (0, 1, 0, 2)] {
        let swap = negative(&(&p[b] * &p[b] - &p[a] * &p[a]));
        let pa = if_nonzero_else(&swap, p[b].clone(), p[a].clone());
        p[b] = if_nonzero_else(&swap, p[a].clone(), p[b].clone());
        p[a] = pa;
        let mc = if_nonzero_else(&swap, m[d].clone(), m[c].clone());
        m[d] = if_nonzero_else(&swap, m[c].clone(), m[d].clone());
        m[c] = mc;
    }
    symbolica::symbol!("__olo_triangle_hv").call(&p.into_iter().chain(m).collect::<Vec<_>>())
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::prelude::Complex;
    #[test]
    fn hv_continuation_regression() {
        super::super::tests::large_stack(|| {
            let x = ["hv_p1", "hv_p2", "hv_p3", "hv_m1", "hv_m2", "hv_m3"]
                .map(|s| symbolica::symbol!(s).to_atom());
            let l = (&x[2] - &x[0] - &x[1]).pow(2) - Atom::num(4) * &x[0] * &x[1];
            let parts = [
                s3_call(
                    &x[0],
                    &x[3],
                    &x[4],
                    &x[4] - &x[5] + &x[1],
                    &x[2] - &x[0] - &x[1],
                    &l,
                ),
                s3_call(
                    &x[2],
                    &x[3],
                    &x[5],
                    -(&x[3] - &x[4]) + &x[2] - &x[1],
                    &x[1] - &x[0] - &x[2],
                    &l,
                ),
                s3_call(
                    &x[1],
                    &x[4],
                    &x[5],
                    -(&x[3] - &x[4]) + &x[2] - &x[1],
                    &x[0] + &x[1] - &x[2],
                    &l,
                ),
            ];
            let mut map = FunctionMap::new();
            sheet_exact::register(&mut map);
            register(&mut map);
            let compile = |parts: &[Atom]| {
                Atom::evaluator_multiple(&parts.iter().map(Atom::as_view).collect::<Vec<_>>(), &x)
                    .function_map(map.clone().into())
                    .direct_translation(true)
                    .build()
                    .unwrap()
                    .map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64()))
            };
            let mut e = compile(&parts);
            let input = [
                Complex::new(-350.9712022329882, 0.),
                Complex::new(-391.7841411757353, 0.),
                Complex::new(-738.0022648960137, 0.),
                Complex::new(342.9713488350603, -7.18128977419155),
                Complex::new(159.2238695848167, -5.483906079793897),
                Complex::new(417.3501586564256, -1.8009017204869213),
            ];
            let mut out = [Complex::new(0., 0.); 3];
            e.evaluate(&input, &mut out);
            let corrections = parts
                .iter()
                .map(|call| {
                    let symbolica::atom::AtomView::Fun(call) = call.as_view() else {
                        panic!("s3 call")
                    };
                    let args = call.iter().map(|a| a.to_owned()).collect::<Vec<_>>();
                    let exact = s3(&args[0], &args[1], &args[2], &args[3], &args[4]);
                    let symbolica::atom::AtomView::Add(sum) = exact.as_view() else {
                        panic!("s3 sum")
                    };
                    sum.iter()
                        .find_map(|term| {
                            if let symbolica::atom::AtomView::Fun(f) = term
                                && f.get_symbol() == Symbol::IF
                            {
                                return Some(f.iter().next().unwrap().to_owned());
                            }
                            None
                        })
                        .unwrap()
                })
                .collect::<Vec<_>>();
            let mut correction_evaluator = compile(&corrections);
            let mut correction_values = [Complex::new(0., 0.); 3];
            correction_evaluator.evaluate(&input, &mut correction_values);
            // Keep this fixture on the active correction branch; the expected
            // continued values below therefore also guard against omitting it.
            for correction in correction_values {
                assert_eq!(correction.re, 0.);
                assert!((correction.im + 2. * std::f64::consts::PI).abs() < 1e-14);
            }
            let expected = [
                (0.2132537392379703, 0.6879516129916667),
                (-0.07964316856308251, 0.3777241269650248),
                (-0.30445018539656665, 0.5897527762006338),
            ];
            for (actual, (re, im)) in out.into_iter().zip(expected) {
                assert!((actual.re - re).hypot(actual.im - im) < 2e-12);
            }
        });
    }
}
