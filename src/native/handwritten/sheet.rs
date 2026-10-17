//! Continued complex values ("sheets") of `sheet_exact.rs`, written directly.
//!
//! A [`Sheet`] is the qmplx-like triple of the exact expressions: a stored value
//! with a nonnegative real part, an integer phase in units of `i*pi`, and the
//! parity of the reflection that produced the stored value. Products, quotients,
//! logarithms, `log(x)/(1-x)`, the continued `Li2(1-x)` and its divided
//! difference follow `sheet_exact::{product_parts, log_body, log_ratio_body,
//! dilog_body, difference_body}` (scalar registration) branch for branch.
use super::super::primitives::{C, NativeFloat, inverse, log, mul, polylog2};
use super::{Constants, is_zero, neg, negative, scale, sign_nonnegative, sub};

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

    /// `log(x)/(1-x)` with the removable value -1 at x=1 on the principal sheet.
    #[inline]
    pub(crate) fn log_over_one_minus(&self, c: &Constants<T>) -> C<T> {
        let difference = sub(&c.one_c, &self.value());
        let logarithm = self.log(c);
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

/// `divided_difference`: (Li2(1-x1) - Li2(1-x2))/(x1-x2), with the first
/// sheet's `log(x)/(1-x)` at coincident physical values.
#[inline]
pub(crate) fn dd<T: NativeFloat>(first: &Sheet<T>, second: &Sheet<T>, c: &Constants<T>) -> C<T> {
    let difference = sub(&first.value(), &second.value());
    if is_zero(&difference) {
        first.log_over_one_minus(c)
    } else {
        mul(
            &sub(&first.dilog(c), &second.dilog(c)),
            &inverse(&difference),
        )
    }
}
