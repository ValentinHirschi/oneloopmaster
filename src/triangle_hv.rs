//! The 't Hooft--Veltman representation for a negative Källén discriminant.
//! All continuation terms are native Symbolica expressions.
use super::*;
use sheet_exact::{SheetAtom as Q, i, im, negative, sign_nonnegative};

pub(crate) fn eta(a: &Atom, b: &Atom, c: &Atom) -> Atom {
    let sa = sign_nonnegative(&im(a));
    let sb = sign_nonnegative(&im(b));
    let sc = sign_nonnegative(&im(c));
    i() * Symbol::PI.to_atom() * &sc * (1 + &sa * &sb) * (1 - sa * sc) / 2
}

fn s3(a: &Atom, m1: &Atom, m2: &Atom, y: &Atom) -> Atom {
    let b = m1 - m2 - a;
    let disc = (&b * &b - Atom::num(4) * a * m2).sqrt();
    let roots = [
        (-&b + &disc) / (Atom::num(2) * a),
        (-&b - &disc) / (Atom::num(2) * a),
    ];
    let term = |r: &Atom| {
        let qy = Q::upper(y - r);
        (&Q::upper(-r) / &qy).dilog() - (&Q::upper(1 - r) / &qy).dilog()
    };
    let z = (y * (a * y + &b) + m2) / a;
    let correction = eta(&(-&roots[0]), &(-&roots[1]), &(m2 / a))
        - eta(&(y - &roots[0]), &(y - &roots[1]), &z)
        - Atom::num(2) * i() * Symbol::PI.to_atom() * negative(a) * negative(&im(&z));
    // The sheet correction is locally constant. Its inactive logarithm may
    // have a removable projective singularity (for example y=0 on the
    // Euclidean sheet), so preserve the same lazy branch used by the other
    // sheet primitives rather than evaluating zero times a nonfinite value.
    term(&roots[0])
        + term(&roots[1])
        + if_nonzero(&correction, &correction * Q::upper((y - 1) / y).log())
}

fn s3_call(a: &Atom, m1: &Atom, m2: &Atom, t1: Atom, t2: Atom, lambda_root: &Atom) -> Atom {
    let alpha = (-t2 + lambda_root) / (Atom::num(2) * a);
    let y = -(t1 + (m1 - m2 - a) * alpha) / lambda_root;
    symbolica::symbol!("__olo_hv_s3").call((a, m1, m2, y))
}

fn body(p: [&Atom; 3], m: [&Atom; 3]) -> Atom {
    let lambda = (p[2] - p[0] - p[1]).pow(2) - Atom::num(4) * p[0] * p[1];
    let root = lambda.sqrt();
    -(s3_call(
        p[0],
        m[0],
        m[1],
        m[1] - m[2] + p[1],
        p[2] - p[0] - p[1],
        &root,
    ) - s3_call(
        p[2],
        m[0],
        m[2],
        -(m[0] - m[1]) + p[2] - p[1],
        p[1] - p[0] - p[2],
        &root,
    ) + s3_call(
        p[1],
        m[1],
        m[2],
        -(m[0] - m[1]) + p[2] - p[1],
        p[0] + p[1] - p[2],
        &root,
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
        symbols[..4].to_vec(),
        s3(&x[0], &x[1], &x[2], &x[3]),
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
            let l = ((&x[2] - &x[0] - &x[1]).pow(2) - Atom::num(4) * &x[0] * &x[1]).sqrt();
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
                    let exact = s3(&args[0], &args[1], &args[2], &args[3]);
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
