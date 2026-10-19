//! Continued complex values ("sheets") of `sheet_exact.rs`, written directly.
//!
//! A [`Sheet`] is the qmplx-like triple of the exact expressions: a stored value
//! with a nonnegative real part, an integer phase in units of `i*pi`, and the
//! parity of the reflection that produced the stored value. Products, quotients,
//! logarithms, `log(x)/(1-x)`, the continued `Li2(1-x)` and its divided
//! difference follow `sheet_exact::{product_parts, log_body, log_ratio_body,
//! dilog_body, difference_body, near_difference}` (scalar registration)
//! branch for branch.
use super::super::primitives::{C, NativeFloat, inverse, log, mul, polylog2};
use super::{Constants, add, is_zero, neg, negative, scale, sign_nonnegative, sub};

#[derive(Clone, Debug)]
pub(crate) struct Sheet<T: NativeFloat> {
    /// Stored value; its real part is never negative (except for -0).
    pub(crate) v: C<T>,
    /// Integer multiple of `i*pi` added to the stored logarithm.
    pub(crate) phase: i32,
    /// Whether the physical value is the negated stored value.
    pub(crate) odd: bool,
}

impl<T: NativeFloat> Sheet<T> {
    /// `SheetAtom::with_sign`: real-axis values take `zero_sign` (+1 upper lip).
    #[inline]
    pub(crate) fn with_sign(value: C<T>, zero_sign: i32) -> Self {
        let n = negative(&value.re);
        if n {
            let s = if value.im.is_fully_zero() {
                zero_sign
            } else {
                sign_nonnegative(&value.im)
            };
            Self {
                v: neg(&value),
                phase: s,
                odd: true,
            }
        } else {
            Self {
                v: value,
                phase: 0,
                odd: false,
            }
        }
    }

    #[inline]
    pub(crate) fn upper(value: C<T>) -> Self {
        Self::with_sign(value, 1)
    }

    #[inline]
    pub(crate) fn lower(value: C<T>) -> Self {
        Self::with_sign(value, -1)
    }

    /// The physical (unreflected) value.
    #[inline]
    pub(crate) fn value(&self) -> C<T> {
        if self.odd { neg(&self.v) } else { self.v.clone() }
    }

    #[inline]
    pub(crate) fn log(&self, c: &Constants<T>) -> C<T> {
        if self.v.im.is_fully_zero() {
            let mut l = log(&self.v);
            if self.phase != 0 {
                l.im += c.pi.clone() * &c.from_int(self.phase);
            }
            l
        } else {
            let phase = self.phase + if self.odd { sign_nonnegative(&self.v.im) } else { 0 };
            let mut l = log(&self.value());
            if phase != 0 {
                l.im += c.pi.clone() * &c.from_int(phase);
            }
            l
        }
    }

    /// `log(x)/(1-x)` with the removable value -1 at x=1 on the principal sheet,
    /// the logarithm next to one through [`log_one_plus`] (`log_ratio_body`).
    #[inline]
    pub(crate) fn log_over_one_minus(&self, c: &Constants<T>) -> C<T> {
        let offset = sub(&self.v, &c.one_c);
        if !self.odd && negative(&(c.four.clone() * &norm(&offset) - &c.one)) {
            let mut l = log_one_plus(&self.v, &offset, c);
            if self.phase != 0 {
                l.im += c.pi.clone() * &c.from_int(self.phase);
            }
            self.over_one_minus(l, c)
        } else {
            self.over_one_minus(self.log(c), c)
        }
    }

    /// The sheet logarithm given the offset `w = x - 1` of the physical
    /// value, formed by the caller without rounding `x` (`log_given_body`):
    /// next to one an unreflected value's logarithm comes from `w` through
    /// [`log_one_plus`].
    #[inline]
    pub(crate) fn log_given(&self, w: &C<T>, c: &Constants<T>) -> C<T> {
        if !self.odd && negative(&(c.four.clone() * &norm(w) - &c.one)) {
            let mut l = log_one_plus(&self.v, w, c);
            if self.phase != 0 {
                l.im += c.pi.clone() * &c.from_int(self.phase);
            }
            l
        } else {
            self.log(c)
        }
    }

    /// `log(x)/(1-x)` with the plain logarithm (the divided difference's
    /// coincident limit, as `log_ratio_body` without `near`).
    #[inline]
    fn plain_log_over_one_minus(&self, c: &Constants<T>) -> C<T> {
        self.over_one_minus(self.log(c), c)
    }

    #[inline]
    fn over_one_minus(&self, logarithm: C<T>, c: &Constants<T>) -> C<T> {
        let difference = sub(&c.one_c, &self.value());
        if !is_zero(&difference) || !is_zero(&logarithm) {
            mul(&logarithm, &inverse(&difference))
        } else {
            C::new(-c.one.clone(), c.one.zero())
        }
    }

    /// Continued `Li2(1-x)` (`dilog_body` with Li2(1)=pi²/6 at x=0, 0 at x=1).
    pub(crate) fn dilog(&self, c: &Constants<T>) -> C<T> {
        if is_zero(&self.v) {
            return C::new(c.pi2_6.clone(), c.one.zero());
        }
        let argument = self.value();
        let nonunit = sub(&argument, &c.one_c);
        if is_zero(&nonunit) {
            return C::new(c.one.zero(), c.one.zero());
        }
        let odd = self.odd as i32;
        let even = self.phase - odd;
        let norm = self.v.re.clone() * &self.v.re - self.v.im.clone() * &(-self.v.im.clone());
        let inverted = negative(&(c.one.clone() - &norm));
        let (normalized, logarithm, even, odd) = if inverted {
            (inverse(&argument), neg(&log(&self.v)), -even, -odd)
        } else {
            (argument, log(&self.v), even, odd)
        };
        let mut sheet_logarithm = logarithm;
        if even + odd != 0 {
            sheet_logarithm.im += c.pi.clone() * &c.from_int(even + odd);
        }
        let one_minus = sub(&c.one_c, &normalized);
        let log_one_minus = log(&one_minus);
        let below_half = negative(&(c.two.clone() * &normalized.re - &c.one));
        let base = if below_half {
            let mut value = neg(&polylog2(&normalized));
            value.re += &c.pi2_6;
            sub(&value, &mul(&log_one_minus, &sheet_logarithm))
        } else {
            let value = polylog2(&one_minus);
            if even != 0 {
                // log(1-z) * i*pi*even
                let f = c.pi.clone() * &c.from_int(even);
                sub(
                    &value,
                    &C::new(-(log_one_minus.im.clone() * &f), log_one_minus.re.clone() * &f),
                )
            } else {
                value
            }
        };
        if inverted {
            let square = mul(&sheet_logarithm, &sheet_logarithm);
            sub(&neg(&base), &scale(&square, &c.half))
        } else {
            base
        }
    }

    #[inline]
    pub(crate) fn mul(&self, right: &Self) -> Self {
        let raw = mul(&self.v, &right.v);
        Self::combine(raw, self.phase + right.phase, self.odd ^ right.odd, right, 1)
    }

    #[inline]
    pub(crate) fn div(&self, right: &Self) -> Self {
        let raw = if is_zero(&sub(&self.v, &right.v)) {
            C::new(self.v.re.one(), self.v.re.zero())
        } else {
            mul(&self.v, &inverse(&right.v))
        };
        Self::combine(raw, self.phase - right.phase, self.odd ^ right.odd, right, -1)
    }

    #[inline]
    fn combine(raw: C<T>, phase: i32, odd: bool, right: &Self, direction: i32) -> Self {
        if negative(&raw.re) {
            Self {
                v: neg(&raw),
                phase: phase + direction * sign_nonnegative(&right.v.im),
                odd: !odd,
            }
        } else {
            Self { v: raw, phase, odd }
        }
    }
}

/// `divided_difference`: (Li2(1-x1) - Li2(1-x2))/(x1-x2), with the plain
/// `log(x)/(1-x)` at coincident physical values, and the forms of
/// `sheet_exact::near_difference` where the two values cancel by more than 16.
#[inline(always)]
pub(crate) fn dd<T: NativeFloat>(first: &Sheet<T>, second: &Sheet<T>, c: &Constants<T>) -> C<T> {
    let difference = sub(&first.value(), &second.value());
    if is_zero(&difference) {
        return first.plain_log_over_one_minus(c);
    }
    let upper = first.dilog(c);
    let numerator = sub(&upper, &second.dilog(c));
    // `negative(256 |numerator|² - |upper|²)` as one comparison: the same
    // decision, since a difference of distinct values is never zero.
    if c.two_fifty_six.clone() * &norm(&numerator) < norm(&upper) {
        if let Some(value) = near_difference(first, second, c) {
            return value;
        }
    }
    mul(&numerator, &inverse(&difference))
}

/// `re(z conj(z))`, which is `re² - im (-im)`: the same value.
#[inline]
fn norm<T: NativeFloat>(z: &C<T>) -> T {
    z.re.clone() * &z.re + z.im.clone() * &z.im
}

/// `sheet_exact::log_one_plus_body`: `log(1 + w)` with its relative precision
/// for small `|w|`, from a real logarithm and the angle of `1 + w`.
#[inline]
fn log_one_plus<T: NativeFloat>(q: &C<T>, w: &C<T>, c: &Constants<T>) -> C<T> {
    let norm = norm(w);
    if !negative(&(c.four.clone() * &norm - &c.one)) {
        return log(q);
    }
    let t = c.two.clone() * &w.re + &norm;
    let u = c.one.clone() + &t;
    let offset = u.clone() - &c.one;
    // The real logarithm of u > 0, and the angle of 1 + w (Re > 0), as the
    // native complex logarithm forms them.
    let log1p = if offset.is_fully_zero() { t } else { u.log() * &t / &offset };
    let z = add(&c.one_c, w);
    let angle = if z.im.is_fully_zero() { z.im.clone() } else { z.im.atan2(&z.re) };
    C::new(log1p * &c.half, angle)
}

/// `log(q)/(1-q)` of `q = 1 + w` on its principal sheet.
#[inline]
fn ratio_log<T: NativeFloat>(q: &C<T>, w: &C<T>, c: &Constants<T>) -> C<T> {
    if is_zero(w) {
        C::new(-c.one.clone(), c.one.zero())
    } else {
        mul(&log_one_plus(q, w, c), &inverse(&neg(w)))
    }
}

/// The sheet logarithm `log(v) + i pi phase`.
#[inline]
fn sheet_log<T: NativeFloat>(s: &Sheet<T>, c: &Constants<T>) -> C<T> {
    let mut l = log(&s.v);
    if s.phase != 0 {
        l.im += c.pi.clone() * &c.from_int(s.phase);
    }
    l
}

/// `sheet_exact::series_body`: `(f(s1) - f(s2))/(s1 - s2)` for
/// `f(s) = Li2(1 - e^-s)`, through `s^degree`.
fn series<T: NativeFloat>(s1: &C<T>, s2: &C<T>, degree: usize, c: &Constants<T>) -> C<T> {
    let mut power = c.one_c.clone();
    let mut complete = c.one_c.clone();
    let mut sum = c.one_c.clone();
    for n in 1..=degree {
        power = mul(&power, s2);
        complete = add(&mul(s1, &complete), &power);
        if n == 1 || n % 2 == 0 {
            sum = add(&sum, &scale(&complete, &c.bernoulli[n / 2]));
        }
    }
    sum
}

/// `sheet_exact::near_difference`; `None` selects the direct difference.
fn near_difference<T: NativeFloat>(
    first: &Sheet<T>,
    second: &Sheet<T>,
    c: &Constants<T>,
) -> Option<C<T>> {
    let one = &c.one;
    let less = |a: T, b: &T| negative(&(a - b));
    let zero_partner = |s: &Sheet<T>| {
        let x = s.value();
        less(c.four.clone() * &norm(&x), one).then(|| {
            let inverse_x = inverse(&x);
            let limit = neg(&mul(&polylog2(&x), &inverse_x));
            sub(
                &limit,
                &mul(
                    &mul(&log_one_plus(&sub(&c.one_c, &x), &neg(&x), c), &inverse_x),
                    &sheet_log(s, c),
                ),
            )
        })
    };
    if is_zero(&second.v) {
        return zero_partner(first);
    }
    if is_zero(&first.v) {
        return zero_partner(second);
    }
    if first.phase != second.phase || first.odd != second.odd {
        return None;
    }
    let [x1, x2] = [first.value(), second.value()];
    let [m1, m2] = [norm(&x1), norm(&x2)];
    // |x|² = |v|² bitwise: Sheet::dilog's inversion test.
    let (inverted_1, inverted_2) = (negative(&(one.clone() - &m1)), negative(&(one.clone() - &m2)));
    let straddle = inverted_1 != inverted_2;
    let inversion = !straddle && inverted_2;
    if straddle && second.phase * (second.phase * second.phase - 1) != 0 {
        return None;
    }
    let (n1, n2, c1, c2, gap_n, inverse_n2) = if inversion {
        let [n1, n2] = [inverse(&x1), inverse(&x2)];
        let c1 = mul(&sub(&x1, &c.one_c), &n1);
        let c2 = mul(&sub(&x2, &c.one_c), &n2);
        let gap_n = mul(&mul(&sub(&x2, &x1), &n1), &n2);
        (n1, n2, c1, c2, gap_n, x2.clone())
    } else {
        let inverse_n2 = inverse(&x2);
        let (c1, c2) = (sub(&c.one_c, &x1), sub(&c.one_c, &x2));
        let gap_n = sub(&x1, &x2);
        (x1.clone(), x2.clone(), c1, c2, gap_n, inverse_n2)
    };
    let below = negative(&(c.two.clone() * &n2.re - one));
    let sign = if inversion { -c.one.clone() } else { c.one.clone() };
    let variable = |n: &C<T>, k: &C<T>, s: &Sheet<T>| {
        if below {
            neg(&log_one_plus(k, &neg(n), c))
        } else {
            neg(&scale(&log(&s.v), &sign))
        }
    };
    let [s1, s2] = [variable(&n1, &c1, first), variable(&n2, &c2, second)];
    let [r1, r2] = [norm(&s1), norm(&s2)];
    let largest = if less(r1.clone(), &r2) { r2 } else { r1 };
    // The first omitted terms of the expansions through s^10 and s^22.
    let base = largest / &c.from_int(39);
    let square = base.clone() * &base;
    let sixth = square.clone() * &square * &square;
    let eight = c.from_int(4) * &c.two;
    let truncated = |power: &T| !(one.clone() + &(eight.clone() * power) - one).is_fully_zero();
    if truncated(&(sixth.clone() * &sixth)) {
        return None;
    }
    let degree = if truncated(&sixth) { 22 } else { 10 };
    let sum = series(&s1, &s2, degree, c);
    let rho = mul(
        &ratio_log(&mul(&n1, &inverse_n2), &mul(&gap_n, &inverse_n2), c),
        &inverse_n2,
    );
    let inverse_c2 = inverse(&c2);
    let sigma = mul(
        &ratio_log(&mul(&c1, &inverse_c2), &neg(&mul(&gap_n, &inverse_c2)), c),
        &inverse_c2,
    );
    let even = (if inversion { -1 } else { 1 }) * (second.phase - second.odd as i32);
    let normalized_difference = if below {
        sub(
            &mul(&sigma, &sub(&sum, &scale(&sheet_log(second, c), &sign))),
            &mul(&s1, &rho),
        )
    } else {
        let difference = mul(&sum, &rho);
        if even != 0 {
            let f = c.pi.clone() * &c.from_int(even);
            // i pi even sigma
            sub(
                &difference,
                &C::new(-(sigma.im.clone() * &f), sigma.re.clone() * &f),
            )
        } else {
            difference
        }
    };
    Some(if inversion {
        let inverse_v2 = inverse(&second.v);
        let ratio = mul(&first.v, &inverse_v2);
        let offset = mul(&sub(&first.v, &second.v), &inverse_v2);
        let logarithms = add(&sheet_log(first, c), &sheet_log(second, c));
        add(
            &mul(&mul(&n1, &n2), &normalized_difference),
            &mul(
                &scale(&logarithms, &c.half),
                &mul(&ratio_log(&ratio, &offset, c), &n2),
            ),
        )
    } else {
        normalized_difference
    })
}
