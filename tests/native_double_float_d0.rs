//! A four-mass D0 with equal real masses (all `m² = 1`) whose native
//! DoubleFloat value had the wrong real part, +0.37226 instead of -0.063421,
//! in one of its eight relabellings and at the rounded arguments NinjaRS
//! passed (n12r12genericm1, lane 2). The imaginary part was right; binary64,
//! the expression interpreter and native Float agreed in every labelling.
//!
//! DoubleFloat products are not commutative in their last bits, so
//! `1 - z conj(z)` kept an imaginary residue of 1e-35. The exact sign test
//! `x - |x|` read it as negative and inverted a dilogarithm argument of
//! modulus 0.09. The native complex product now orders its real products.
#![cfg(feature = "generated-evaluators")]
use oneloop::{NativeEvaluator, NativeFloat, ScalarIntegral};
use symbolica::domains::{
    float::{Complex, DoubleFloat, Float, Real, RealLike, SingleFloat},
    rational::Rational,
};

/// `[p1², p2², p3², p4², s12, s23]`; the masses and `mu²` are 1.
const INVARIANTS: [(i64, i64); 6] = [
    (2171, 144),
    (211, 16),
    (-5, 4),
    (-587, 72),
    (-311, 72),
    (27, 16),
];

/// The finite part, from native Float and the expression interpreter at 256
/// bits (equal to all digits shown in every labelling).
const REFERENCE: [&str; 2] = [
    "-6.342090417471351697878784160261216715137e-2",
    "2.110820888749008107801505116588042554555e-2",
];

/// The rounded DoubleFloat invariants of the first labelling as NinjaRS
/// formed them from its momenta (exact decimal values of the pairs).
const ROUNDED: [&str; 6] = [
    "1.507638888888888888888888888888891080169e1",
    "1.318750000000000000000000000000001232595e1",
    "-1.249999999999999999999999999999981511073e0",
    "-8.152777777777777777777777777777821603384e0",
    "-4.319444444444444444444444444444430748943e0",
    "1.687499999999999999999999999999996918512e0",
];

/// The eight images of a labelling under the box symmetries: cyclic shifts
/// `(p1, p2, p3, p4, s12, s23) -> (p2, p3, p4, p1, s23, s12)` and the
/// reflection `-> (p3, p2, p1, p4, s23, s12)`.
fn labellings<V: Clone>(p: [V; 6]) -> Vec<[V; 6]> {
    let shift = |p: &[V; 6]| [1, 2, 3, 0, 5, 4].map(|i| p[i].clone());
    let reflect = |p: &[V; 6]| [2, 1, 0, 3, 5, 4].map(|i| p[i].clone());
    let mut images = Vec::new();
    let mut current = p;
    for _ in 0..4 {
        images.push(reflect(&current));
        images.push(current.clone());
        current = shift(&current);
    }
    images
}

fn finite<T: NativeFloat>(invariants: &[T; 6], bits: u32) -> Complex<T> {
    let mut evaluator =
        NativeEvaluator::<T>::with_binary_precision(ScalarIntegral::D0, bits).unwrap();
    let one = invariants[0].one();
    let mut input: Vec<_> = invariants
        .iter()
        .map(|p| Complex::new(p.clone(), one.zero()))
        .collect();
    input.extend((0..5).map(|_| Complex::new(one.clone(), one.zero())));
    let mut output = vec![Complex::new(one.zero(), one.zero()); 3];
    evaluator.evaluate(&input, &mut output).unwrap();
    output[0].clone()
}

fn to_float(value: &DoubleFloat) -> Float {
    let pair = value.into_inner();
    Float::with_val(256, pair.hi()) + Float::with_val(256, pair.lo())
}

fn relative_error(value: [Float; 2], reference: &[Float; 2]) -> f64 {
    let [re, im] = value;
    let (dre, dim) = (re - &reference[0], im - &reference[1]);
    let error = (dre.clone() * &dre + dim.clone() * &dim).sqrt();
    let scale =
        (reference[0].clone() * &reference[0] + reference[1].clone() * &reference[1]).sqrt();
    (error / scale).to_f64()
}

#[test]
fn equal_mass_four_mass_box_in_double_float_matches_every_labelling() {
    let reference = REFERENCE.map(|value| Float::parse(value, Some(256)).unwrap());
    let exact = INVARIANTS.map(|(n, d)| Rational::from((n, d)));
    let rounded = ROUNDED.map(|value| Float::parse(value, Some(256)).unwrap().to_double_float());
    let mut failures = Vec::new();
    let mut check = |route: String, value: [Float; 2], bound: f64| {
        let error = relative_error(value.clone(), &reference);
        if error.is_nan() || error > bound {
            failures.push(format!(
                "{route}: {:e} {:e} (error {error:e})",
                value[0], value[1]
            ));
        }
    };
    for (index, labelling) in labellings(exact).into_iter().enumerate() {
        let double = labelling.clone().map(|p| DoubleFloat::from(&p));
        let value = finite(&double, 106);
        check(
            format!("DoubleFloat labelling {index}"),
            [to_float(&value.re), to_float(&value.im)],
            1e-29,
        );
        let float = labelling.clone().map(|p| Float::new(128).from_rational(&p));
        let value = finite(&float, 128);
        check(
            format!("Float(128) labelling {index}"),
            [value.re, value.im],
            1e-34,
        );
        let binary = labelling.map(|p| p.to_f64());
        let value = finite(&binary, 53);
        check(
            format!("binary64 labelling {index}"),
            [
                Float::with_val(256, value.re),
                Float::with_val(256, value.im),
            ],
            1e-13,
        );
    }
    // The rounded invariants move the value by about 1e-31.
    for (index, labelling) in labellings(rounded).into_iter().enumerate() {
        let value = finite(&labelling, 106);
        check(
            format!("rounded DoubleFloat labelling {index}"),
            [to_float(&value.re), to_float(&value.im)],
            1e-29,
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
