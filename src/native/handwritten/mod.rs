//! Hand-written binary64 kernels for the dominant C0/D0 mass sectors.
//!
//! Each kernel is a direct transcription of the exact expressions, branch for
//! branch: the three-mass triangle of `triangle.rs` (OneLOop `trif3`), the
//! contour box of `box_complex.rs` (`boxc`, used for three and four massive
//! propagators where OneLOop selects it) and the four-mass Denner-Nierste-
//! Scharf box of `box_integral.rs` (`boxf4`), with the same exact zero tests,
//! comparisons, infinitesimal directions and sheets (`sheet_exact.rs`, mirrored
//! by [`sheet::Sheet`]). Unlike the generated code they evaluate one sector
//! dispatch per call, never evaluate inactive branches, and keep real values
//! and predicates real. Sectors not covered return `None` and use the
//! generated code, which also remains the implementation for `DoubleFloat` and
//! `Float`: those types track precision through the expression's own grouping.
//!
//! The triangle kernel is generic over [`NativeFloat`]; the box kernels are
//! binary64 only.
use super::primitives::{C, NativeFloat, inverse, mul, sqrt};

mod box_contour;
mod boxes;
pub(crate) mod sheet;
pub(crate) mod triangle;

/// Working-precision constants, prepared once per evaluator.
#[derive(Clone)]
pub(crate) struct Constants<T: NativeFloat> {
    pub(crate) one: T,
    pub(crate) two: T,
    pub(crate) four: T,
    pub(crate) half: T,
    pub(crate) pi: T,
    pub(crate) pi2_6: T,
    pub(crate) one_c: C<T>,
    pub(crate) two_fifty_six: T,
    small_ints: [T; 9],
    /// `B_n/(n+1)!` for `n = 1, 2, 4, ..., 22` (`sheet_exact::series_body`).
    pub(crate) bernoulli: [T; 12],
}

impl<T: NativeFloat> Constants<T> {
    pub(crate) fn new(prototype: &T) -> Self {
        let one = prototype.one();
        let two = prototype.from_i64(2);
        let pi = prototype.pi();
        let pi2_6 = pi.clone() * &pi / &prototype.from_i64(6);
        // B_n/(n+1)! for n = 1, 2, 4, ..., 22, as exact rationals.
        const BERNOULLI: [(i64, i64); 12] = [
            (-1, 2),
            (1, 6),
            (-1, 30),
            (1, 42),
            (-1, 30),
            (5, 66),
            (-691, 2730),
            (7, 6),
            (-3617, 510),
            (43867, 798),
            (-174611, 330),
            (854513, 138),
        ];
        let bernoulli = core::array::from_fn(|k| {
            let n = if k == 0 { 1 } else { 2 * k as u32 };
            let (numerator, denominator) = BERNOULLI[k];
            let factorial = symbolica::domains::integer::Integer::factorial(n + 1);
            let value = symbolica::domains::rational::Rational::from((numerator, denominator))
                / &symbolica::domains::rational::Rational::from(factorial);
            prototype.from_rational(&value)
        });
        Self {
            bernoulli,
            half: one.clone() / &two,
            four: prototype.from_i64(4),
            two_fifty_six: prototype.from_i64(256),
            one_c: C::new(one.clone(), prototype.zero()),
            small_ints: core::array::from_fn(|i| prototype.from_i64(i as i64 - 4)),
            one,
            two,
            pi,
            pi2_6,
        }
    }

    #[inline]
    pub(crate) fn from_int(&self, value: i32) -> T {
        if (-4..=4).contains(&value) {
            self.small_ints[(value + 4) as usize].clone()
        } else {
            self.one.from_i64(value as i64)
        }
    }
}

#[inline]
pub(crate) fn is_zero<T: NativeFloat>(z: &C<T>) -> bool {
    z.re.is_fully_zero() && z.im.is_fully_zero()
}

/// `sheet_exact::negative` of a real value: `x - |x|` is nonzero (x < 0,
/// +-inf or NaN, as the expression's own test).
#[inline]
pub(crate) fn negative<T: NativeFloat>(x: &T) -> bool {
    !(x.clone() - &x.norm()).is_fully_zero()
}

/// `sheet_exact::sign_nonnegative` of a real value: -1 iff `x - |x|` is nonzero.
#[inline]
pub(crate) fn sign_nonnegative<T: NativeFloat>(x: &T) -> i32 {
    if negative(x) { -1 } else { 1 }
}

#[inline]
pub(crate) fn add<T: NativeFloat>(a: &C<T>, b: &C<T>) -> C<T> {
    C::new(a.re.clone() + &b.re, a.im.clone() + &b.im)
}

#[inline]
pub(crate) fn sub<T: NativeFloat>(a: &C<T>, b: &C<T>) -> C<T> {
    C::new(a.re.clone() - &b.re, a.im.clone() - &b.im)
}

#[inline]
pub(crate) fn neg<T: NativeFloat>(a: &C<T>) -> C<T> {
    C::new(-a.re.clone(), -a.im.clone())
}

#[inline]
pub(crate) fn scale<T: NativeFloat>(a: &C<T>, f: &T) -> C<T> {
    C::new(a.re.clone() * f, a.im.clone() * f)
}

#[inline]
pub(crate) fn div<T: NativeFloat>(a: &C<T>, b: &C<T>) -> C<T> {
    mul(a, &inverse(b))
}

#[inline]
pub(crate) fn real<T: NativeFloat>(x: T) -> C<T> {
    let zero = x.zero();
    C::new(x, zero)
}

/// `physical_sqrt` (`sqrt_lower`): principal off the real axis, lower lip on
/// the negative real axis.
#[inline]
pub(crate) fn physical_sqrt<T: NativeFloat>(z: &C<T>) -> C<T> {
    if !z.im.is_fully_zero() {
        sqrt(z)
    } else if negative(&z.re) {
        C::new(z.re.zero(), -(-z.re.clone()).sqrt())
    } else {
        C::new(z.re.clone().sqrt(), z.im.clone())
    }
}

/// `triangle::r_function`: the root r of r + 1/r = k selected like OneLOop's
/// rfun, and r - 1/r, taken directly from the square root added to k (it is
/// exactly zero at a threshold and keeps its sign next to one).
#[inline]
pub(crate) fn r_function<T: NativeFloat>(k: &C<T>, c: &Constants<T>) -> (C<T>, C<T>) {
    if !k.im.is_fully_zero() {
        let two = real(c.two.clone());
        let difference = mul(&sqrt(&sub(k, &two)), &sqrt(&add(k, &two)));
        (scale(&add(k, &difference), &c.half), difference)
    } else {
        let x = &k.re;
        let discriminant = x.clone() * x - &c.four;
        if negative(&discriminant) {
            let root = (-discriminant).sqrt();
            (
                C::new(x.clone() * &c.half, -(root.clone() * &c.half)),
                C::new(x.zero(), -root),
            )
        } else {
            let root = discriminant.sqrt();
            let difference = if negative(x) { -root } else { root };
            (
                C::new((x.clone() + &difference) * &c.half, x.zero()),
                C::new(difference, x.zero()),
            )
        }
    }
}

/// C0 `[finite, pole, double pole]` for the covered sectors: three massive
/// propagators (`register_triangles` sector 7) with nonzero momenta and a
/// nonnegative Källén function (OneLOop trif3).
#[inline]
pub(crate) fn c0<T: NativeFloat>(a: &[C<T>], c: &Constants<T>) -> Option<[C<T>; 3]> {
    if a[3..6].iter().any(is_zero) {
        return None;
    }
    let finite = triangle::three_masses([&a[0], &a[1], &a[2]], [&a[3], &a[4], &a[5]], c)?;
    let zero = real(c.one.zero());
    Some([finite, zero.clone(), zero])
}

/// D0 `[finite, pole, double pole]` for the covered sectors (binary64 only).
#[inline]
pub(crate) fn d0<T: NativeFloat>(a: &[C<T>], c: &Constants<T>) -> Option<[C<T>; 3]> {
    use core::any::Any;
    let a: &[C<T>; 11] = a.try_into().ok()?;
    let a = (a as &dyn Any).downcast_ref::<[C<f64>; 11]>()?;
    let c = (c as &dyn Any).downcast_ref::<Constants<f64>>()?;
    let finite = boxes::d0(a, c)?;
    let zero = C::new(0.0, 0.0);
    let values: [C<f64>; 3] = [finite, zero, zero];
    (&values as &dyn Any).downcast_ref::<[C<T>; 3]>().cloned()
}
