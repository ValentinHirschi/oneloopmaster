//! Binary64 bubble and tadpole kernels.
//!
//! The algorithm of the [parent module](super), specialized to `f64`: the same
//! sectors and exact zero tests, the same mass ordering, lips, root pair and
//! series radius, the same `bnlog` representation for every rank. What
//! changes is how the values are formed:
//!
//! - Real kinematics stay real. A real root `x` outside `[0, 1]` has a real
//!   `y = log(1-1/x)`, formed by one `ln_1p`; on the lip (`0 < x < 1`) the
//!   imaginary part `±pi (1 - x^(n+1))` is added to a real direct formula.
//! - For real masses every root and `-1/x` is one quotient of the inputs
//!   (`p²/(m1-p²)`, `(m0-m1)/m1`, `-2p²/(Δ-b)`, ...), not a product with a
//!   reciprocal. Complex masses keep the generic code's roots exactly.
//! - Real masses below threshold have complex-conjugate roots, for which
//!   `bnlog(n, x1) + bnlog(n, x2) = 2 Re bnlog(n, x1)` (no root lies on the
//!   cut). One complex `bnlog` replaces two, and the sum is exactly real.
//! - The series of `bnlog` is evaluated with Estrin's scheme instead of
//!   Horner's rule, on the first 8 coefficients of each rank for
//!   `|y|_1 <= 0.01`, the first 16 for `|y|_1 <= 0.2` (enough there for the
//!   `2^-57` truncation bound of ranks zero to three, see [`Binary64::new`]),
//!   and otherwise on all of them (15 to 24 per rank, zero padded to 16 or
//!   24). Rank four needs 9 and 17 coefficients at those radii, so it keeps
//!   the same three lengths with its own, smaller radii (`0.0085` and `0.17`;
//!   8 and 16 coefficients reach `0.0090` and `0.1736` there). Horner's rule
//!   over the minimal count is a dependent chain of up to 23 complex
//!   multiply-adds per rank; Estrin's has 3 to 5 levels, and real `y` needs
//!   no complex arithmetic at all.
//!
//! Each rank's value depends only on the row and that rank, so the rank-two
//! coefficients of `B111` are bitwise those of `B11`, and the rank-three
//! coefficients of `B1111` bitwise those of `B111`. Rank four changes nothing
//! in the `B11` and `B111` kernels: its series length is chosen apart from the
//! others only when `RANK == 4` (a compile-time condition), their per-rank
//! arrays keep their lengths (the `W` and `P` parameters of
//! [`Binary64::bubble`]), and the rank-four kernel is not inlined into the
//! callers of the others. Their machine code is that of the crate before rank
//! four (compared instruction by instruction on x86-64).
use super::{BnlogSeries, MAX_BUBBLE_RANK, SERIES_RADIUS};
use crate::native::primitives::{C, inverse, sqrt};
use std::f64::consts::PI;

type C64 = C<f64>;

/// Coefficients evaluated for `|y|_1 <= TINY_RADIUS`.
const TINY: usize = 8;
/// Coefficients evaluated for `|y|_1 <= SHORT_RADIUS`.
const SHORT: usize = 16;
/// Coefficients evaluated otherwise (zero padded).
const LONG: usize = 24;
/// Largest `|y|_1` for which 8 coefficients meet the bound of ranks 0 to 3.
const TINY_RADIUS: f64 = 0.01;
/// Largest `|y|_1` for which 16 coefficients meet the bound of ranks 0 to 3.
const SHORT_RADIUS: f64 = 0.2;
/// [`TINY_RADIUS`] of rank four (8 coefficients reach `|y|_1 = 0.0090`).
const TINY_RADIUS_4: f64 = 0.0085;
/// [`SHORT_RADIUS`] of rank four (16 coefficients reach `|y|_1 = 0.1736`).
const SHORT_RADIUS_4: f64 = 0.17;

/// How many coefficients of the series a row needs.
#[derive(Clone, Copy)]
enum Length {
    Tiny,
    Short,
    Long,
}

impl Length {
    #[inline(always)]
    fn of(size: f64) -> Self {
        if size <= TINY_RADIUS {
            Self::Tiny
        } else if size <= SHORT_RADIUS {
            Self::Short
        } else {
            Self::Long
        }
    }

    /// [`Self::of`] for rank four.
    #[inline(always)]
    fn of_rank_four(size: f64) -> Self {
        if size <= TINY_RADIUS_4 {
            Self::Tiny
        } else if size <= SHORT_RADIUS_4 {
            Self::Short
        } else {
            Self::Long
        }
    }

    /// The length of rank `n` of a `RANK` family, given `length = of(size)`:
    /// `of(size)` except for rank four (resolved at compile time for
    /// `RANK <= 3`).
    #[inline(always)]
    fn of_rank<const RANK: usize>(n: usize, size: f64, length: Self) -> Self {
        if RANK == 4 && n == 4 {
            Self::of_rank_four(size)
        } else {
            length
        }
    }
}

const ZERO: C64 = C64 { re: 0.0, im: 0.0 };
/// Ranks `0..=MAX_BUBBLE_RANK`.
const RANKS: usize = MAX_BUBBLE_RANK + 1;
/// `1/(n+1)`, `1/(n+1)^2`, `n+1` and the harmonic numbers `H_(n+1)`.
const INVERSE: [f64; RANKS] = [1.0, 0.5, 1.0 / 3.0, 0.25, 0.2];
const INVERSE_SQUARED: [f64; RANKS] = [1.0, 0.25, 1.0 / 9.0, 0.0625, 0.04];
const COUNTS: [f64; RANKS] = [1.0, 2.0, 3.0, 4.0, 5.0];
const HARMONIC: [f64; RANKS] = [1.0, 1.5, 11.0 / 6.0, 25.0 / 12.0, 137.0 / 60.0];
/// The `1/epsilon` coefficients of `B0, B1, B11, B111, B1111`.
const POLES: [f64; RANKS] = [1.0, -0.5, 1.0 / 3.0, -0.25, 0.2];

#[inline(always)]
fn complex(re: f64, im: f64) -> C64 {
    C64 { re, im }
}

#[inline(always)]
fn is_zero(z: C64) -> bool {
    z.re == 0.0 && z.im == 0.0
}

#[inline(always)]
fn size(z: C64) -> f64 {
    z.re.abs() + z.im.abs()
}

#[inline(always)]
fn add(a: C64, b: C64) -> C64 {
    complex(a.re + b.re, a.im + b.im)
}

#[inline(always)]
fn sub(a: C64, b: C64) -> C64 {
    complex(a.re - b.re, a.im - b.im)
}

#[inline(always)]
fn neg(a: C64) -> C64 {
    complex(-a.re, -a.im)
}

#[inline(always)]
fn scale(a: C64, factor: f64) -> C64 {
    complex(a.re * factor, a.im * factor)
}

#[inline(always)]
fn mul(a: C64, b: C64) -> C64 {
    complex(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
}

/// `log(m/mu²)` on the Feynman lip (real negative ratios take `-i pi`).
#[inline(always)]
fn mass_log(mass: C64, mu_squared: f64) -> C64 {
    if mass.im == 0.0 {
        let ratio = mass.re / mu_squared;
        if ratio < 0.0 {
            complex((-ratio).ln(), -PI)
        } else {
            complex(ratio.ln(), 0.0)
        }
    } else {
        // The generic complex logarithm, as in the parent module.
        crate::native::primitives::log(&complex(mass.re / mu_squared, mass.im / mu_squared))
    }
}

/// Principal `log(1+w)` of a complex `w`, as the parent module's `log1p`.
#[inline(always)]
fn log1p(w: C64) -> C64 {
    let real = 1.0 + w.re;
    let modulus = if w.re.abs() < 0.5 && w.im.abs() < 0.5 {
        (w.re * (2.0 + w.re) + w.im * w.im).ln_1p() * 0.5
    } else {
        <f64 as crate::native::NativeFloat>::hypot(&real, &w.im).ln()
    };
    complex(modulus, w.im.atan2(real))
}

/// Powers `y, y², y⁴, y⁸, y¹⁶` of Estrin's scheme.
#[derive(Clone, Copy)]
struct Powers<T>([T; 5]);

impl Powers<f64> {
    #[inline(always)]
    fn new(y: f64) -> Self {
        let y2 = y * y;
        let y4 = y2 * y2;
        let y8 = y4 * y4;
        Self([y, y2, y4, y8, y8 * y8])
    }
}

impl Powers<C64> {
    #[inline(always)]
    fn new(y: C64) -> Self {
        let y2 = mul(y, y);
        let y4 = mul(y2, y2);
        let y8 = mul(y4, y4);
        Self([y, y2, y4, y8, mul(y8, y8)])
    }
}

/// One Estrin level over real values: `a_i = c_(2i) + c_(2i+1) q`.
#[inline(always)]
fn level<const N: usize, const M: usize>(c: &[f64; N], q: f64) -> [f64; M] {
    let mut a = [0.0; M];
    for i in 0..M {
        a[i] = c[2 * i] + c[2 * i + 1] * q;
    }
    a
}

/// `sum_k c_k y^k` over 8 real coefficients at real `y` (Estrin's scheme).
#[inline(always)]
fn estrin_tiny(c: &[f64; TINY], y: &Powers<f64>) -> f64 {
    let [y1, y2, y4, _, _] = y.0;
    let a: [f64; 4] = level(c, y1);
    (a[0] + a[1] * y2) + (a[2] + a[3] * y2) * y4
}

/// `sum_k c_k y^k` over 16 real coefficients at real `y` (Estrin's scheme).
#[inline(always)]
fn estrin_short(c: &[f64; SHORT], y: &Powers<f64>) -> f64 {
    let [y1, y2, y4, y8, _] = y.0;
    let a: [f64; 8] = level(c, y1);
    let b: [f64; 4] = level(&a, y2);
    (b[0] + b[1] * y4) + (b[2] + b[3] * y4) * y8
}

/// `sum_k c_k y^k` over 24 real coefficients at real `y` (Estrin's scheme).
#[inline(always)]
fn estrin_long(c: &[f64; LONG], y: &Powers<f64>) -> f64 {
    let [y1, y2, y4, y8, y16] = y.0;
    let a: [f64; 12] = level(c, y1);
    let b: [f64; 6] = level(&a, y2);
    let d: [f64; 3] = level(&b, y4);
    (d[0] + d[1] * y8) + d[2] * y16
}

/// First Estrin level at complex `y`: `c_(2i) + c_(2i+1) y`.
#[inline(always)]
fn pair(c0: f64, c1: f64, y: C64) -> C64 {
    complex(c0 + c1 * y.re, c1 * y.im)
}

/// `a + b q` for complex values.
#[inline(always)]
fn fold(a: C64, b: C64, q: C64) -> C64 {
    add(a, mul(b, q))
}

/// The first Estrin level at complex `y` (real coefficients).
#[inline(always)]
fn first_level<const N: usize, const M: usize>(c: &[f64; N], y: C64) -> [C64; M] {
    let mut a = [ZERO; M];
    for i in 0..M {
        a[i] = pair(c[2 * i], c[2 * i + 1], y);
    }
    a
}

/// A further Estrin level over complex values.
#[inline(always)]
fn complex_level<const N: usize, const M: usize>(c: &[C64; N], q: C64) -> [C64; M] {
    let mut a = [ZERO; M];
    for i in 0..M {
        a[i] = fold(c[2 * i], c[2 * i + 1], q);
    }
    a
}

/// [`estrin_tiny`] at complex `y`.
#[inline(always)]
fn estrin_tiny_complex(c: &[f64; TINY], y: &Powers<C64>) -> C64 {
    let [y1, y2, y4, _, _] = y.0;
    let a: [C64; 4] = first_level(c, y1);
    fold(fold(a[0], a[1], y2), fold(a[2], a[3], y2), y4)
}

/// [`estrin_short`] at complex `y`.
#[inline(always)]
fn estrin_short_complex(c: &[f64; SHORT], y: &Powers<C64>) -> C64 {
    let [y1, y2, y4, y8, _] = y.0;
    let a: [C64; 8] = first_level(c, y1);
    let b: [C64; 4] = complex_level(&a, y2);
    fold(fold(b[0], b[1], y4), fold(b[2], b[3], y4), y8)
}

/// [`estrin_long`] at complex `y`.
#[inline(always)]
fn estrin_long_complex(c: &[f64; LONG], y: &Powers<C64>) -> C64 {
    let [y1, y2, y4, y8, y16] = y.0;
    let a: [C64; 12] = first_level(c, y1);
    let b: [C64; 6] = complex_level(&a, y2);
    let d: [C64; 3] = complex_level(&b, y4);
    fold(fold(d[0], d[1], y8), d[2], y16)
}

/// Binary64 constants of the tensor kernels: the series of `bnlog` per rank.
#[derive(Clone)]
pub(super) struct Binary64 {
    /// The first [`TINY`] coefficients of each rank.
    tiny: [[f64; TINY]; RANKS],
    /// The first [`SHORT`] coefficients of each rank (zero padded).
    short: [[f64; SHORT]; RANKS],
    /// All coefficients of each rank (zero padded).
    long: [[f64; LONG]; RANKS],
}

impl Binary64 {
    /// The coefficients of the parent module's series at 53 bits. Checks the
    /// three truncations: with all coefficients every rank meets its bound up
    /// to [`SERIES_RADIUS`], with [`SHORT`] up to [`SHORT_RADIUS`] and with
    /// [`TINY`] up to [`TINY_RADIUS`] (for rank four up to [`SHORT_RADIUS_4`]
    /// and [`TINY_RADIUS_4`]).
    pub(super) fn new() -> Result<Self, String> {
        let mut tiny = [[0.0; TINY]; RANKS];
        let mut short = [[0.0; SHORT]; RANKS];
        let mut long = [[0.0; LONG]; RANKS];
        for rank in 0..=MAX_BUBBLE_RANK {
            let (short_radius, tiny_radius) = if rank == 4 {
                (SHORT_RADIUS_4, TINY_RADIUS_4)
            } else {
                (SHORT_RADIUS, TINY_RADIUS)
            };
            let series = BnlogSeries::<f64>::new(rank, 53, &0.0)?;
            let needed = |radius: f64| {
                series
                    .radii
                    .iter()
                    .position(|&bound| radius <= bound)
                    .map_or(series.coefficients.len(), |index| index + 1)
            };
            if series.coefficients.len() > LONG
                || needed(SERIES_RADIUS) > series.coefficients.len()
                || needed(short_radius) > SHORT
                || needed(tiny_radius) > TINY
            {
                return Err("binary64 bnlog series exceeds its fixed length".into());
            }
            for (k, coefficient) in series.coefficients.iter().enumerate() {
                long[rank][k] = *coefficient;
                if k < SHORT {
                    short[rank][k] = *coefficient;
                }
                if k < TINY {
                    tiny[rank][k] = *coefficient;
                }
            }
        }
        Ok(Self { tiny, short, long })
    }

    /// `S_n(y)` of rank `n` (`bnlog(n,x) = -(n+1) y (xy)^(n+1) S_n(y)`).
    #[inline(always)]
    fn series(&self, rank: usize, y: &Powers<f64>, length: Length) -> f64 {
        match length {
            Length::Tiny => estrin_tiny(&self.tiny[rank], y),
            Length::Short => estrin_short(&self.short[rank], y),
            Length::Long => estrin_long(&self.long[rank], y),
        }
    }

    #[inline(always)]
    fn series_complex(&self, rank: usize, y: &Powers<C64>, length: Length) -> C64 {
        match length {
            Length::Tiny => estrin_tiny_complex(&self.tiny[rank], y),
            Length::Short => estrin_short_complex(&self.short[rank], y),
            Length::Long => estrin_long_complex(&self.long[rank], y),
        }
    }

    /// `bnlog(n, x)` for `n = 0..=RANK` and real `x != 0`, given `w = -1/x`
    /// (formed by the caller from the inputs). For `0 < x < 1`, `upper_lip`
    /// selects `+i pi` for `log(1-1/x)`. `W` and `P` are the array lengths
    /// of [`Self::bubble`].
    #[inline(always)]
    fn bnlog_real<const RANK: usize, const W: usize, const P: usize>(
        &self,
        x: f64,
        w: f64,
        upper_lip: bool,
        out: &mut [C64; W],
    ) {
        if x == 1.0 {
            for (value, harmonic) in out.iter_mut().zip(HARMONIC).take(RANK + 1) {
                *value = complex(-harmonic, 0.0);
            }
            return;
        }
        let lip = x > 0.0 && x < 1.0;
        let (y, angle) = if lip {
            (((1.0 - x) / x).ln(), if upper_lip { PI } else { -PI })
        } else {
            (w.ln_1p(), 0.0)
        };
        if y == 0.0 && angle == 0.0 {
            // 1/x underflowed: the x -> infinity limit.
            for value in out.iter_mut().take(RANK + 1) {
                *value = ZERO;
            }
            return;
        }
        if !lip && y.abs() < SERIES_RADIUS {
            let length = Length::of(y.abs());
            let powers = Powers::<f64>::new(y);
            let xy = x * y;
            let mut power = xy;
            for (n, value) in out.iter_mut().enumerate().take(RANK + 1) {
                if n > 0 {
                    power *= xy;
                }
                let length = Length::of_rank::<RANK>(n, y.abs(), length);
                let series = self.series(n, &powers, length);
                *value = complex(-COUNTS[n] * (series * y * power), 0.0);
            }
        } else {
            // (1 - x^(n+1)) y - sum_(k=0..n) x^k/(n+1-k), y = Y + i angle.
            let mut powers = [1.0; P];
            for k in 1..=RANK + 1 {
                powers[k] = powers[k - 1] * x;
            }
            for (n, value) in out.iter_mut().enumerate().take(RANK + 1) {
                let mut polynomial = 0.0;
                for (k, power) in powers.iter().enumerate().take(n + 1) {
                    polynomial += power * INVERSE[n - k];
                }
                let one_minus = 1.0 - powers[n + 1];
                *value = complex(one_minus * y - polynomial, one_minus * angle);
            }
        }
    }

    /// `bnlog(n, x)` for `n = 0..=RANK` and complex `x` (a real `x` takes
    /// [`Self::bnlog_real`] with its lip).
    #[inline(always)]
    fn bnlog<const RANK: usize, const W: usize, const P: usize>(
        &self,
        x: C64,
        upper_lip: bool,
        out: &mut [C64; W],
    ) {
        if x.im == 0.0 {
            return self.bnlog_real::<RANK, W, P>(x.re, -1.0 / x.re, upper_lip, out);
        }
        let reciprocal = inverse(&x);
        let y = log1p(complex(-reciprocal.re, -reciprocal.im));
        if is_zero(y) {
            for value in out.iter_mut().take(RANK + 1) {
                *value = ZERO;
            }
            return;
        }
        let magnitude = size(y);
        if magnitude < SERIES_RADIUS {
            let length = Length::of(magnitude);
            let powers = Powers::<C64>::new(y);
            let xy = mul(x, y);
            let mut power = xy;
            for (n, value) in out.iter_mut().enumerate().take(RANK + 1) {
                if n > 0 {
                    power = mul(power, xy);
                }
                let length = Length::of_rank::<RANK>(n, magnitude, length);
                let series = self.series_complex(n, &powers, length);
                *value = scale(mul(mul(series, y), power), -COUNTS[n]);
            }
        } else {
            let mut powers = [complex(1.0, 0.0); P];
            for k in 1..=RANK + 1 {
                powers[k] = mul(powers[k - 1], x);
            }
            for (n, value) in out.iter_mut().enumerate().take(RANK + 1) {
                let mut polynomial = ZERO;
                for (k, power) in powers.iter().enumerate().take(n + 1) {
                    polynomial = add(polynomial, scale(*power, INVERSE[n - k]));
                }
                let one_minus = complex(1.0 - powers[n + 1].re, -powers[n + 1].im);
                *value = sub(mul(one_minus, y), polynomial);
            }
        }
    }

    /// `[B0, B1, B11, B111, B1111]`, [`Self::bubble`] at rank four. Out of
    /// line, so that the callers of the lower ranks keep their code.
    #[inline(never)]
    pub(super) fn bubble_rank_four(&self, input: &[C64; 4], output: &mut [C64]) {
        self.bubble::<4, 5, 6>(input, output);
    }

    /// `[B0, B1, B11, B111, B1111]` through `RANK`; see the parent's `bubble`.
    /// `W` (at least `RANK+1`) is the length of the per-rank arrays and `P =
    /// W+1` that of the root powers: 4 and 5 through rank three (`B11` as
    /// `bubble::<2, 4, 5>`, `B111` as `bubble::<3, 4, 5>`), as before rank
    /// four existed, so that those kernels compile to the same code.
    #[inline]
    pub(super) fn bubble<const RANK: usize, const W: usize, const P: usize>(
        &self,
        input: &[C64; 4],
        output: &mut [C64],
    ) {
        const { assert!(RANK < W && W <= RANKS && P == W + 1) };
        debug_assert!(output.len() == 3 * (RANK + 1));
        let momentum = input[0].re;
        let mu_squared = input[3].re;
        if momentum == 0.0 && is_zero(input[1]) && is_zero(input[2]) {
            for value in output.iter_mut() {
                *value = ZERO;
            }
            return;
        }
        // OneLOop keeps the order when |m1| >= |m0| and otherwise computes
        // B(p; m1, m0), reconstructing the requested order at the end.
        let switch = size(input[2]) < size(input[1]);
        let (m0, m1) = if switch {
            (input[2], input[1])
        } else {
            (input[1], input[2])
        };
        let real_masses = m0.im == 0.0 && m1.im == 0.0;
        let mut roots = [ZERO; W];
        let mut finite = [ZERO; W];
        let reconstruct;
        if momentum == 0.0 {
            let lna = neg(mass_log(m1, mu_squared));
            if m0 == m1 {
                for (n, value) in finite.iter_mut().enumerate().take(RANK + 1) {
                    *value = scale(lna, INVERSE[n]);
                }
            } else {
                // x1 = m1/(m1-m0), on the lip opposite to Re(m1-m0); a
                // massless m0 is the unit root.
                let difference = sub(m1, m0);
                let lip = difference.re < 0.0;
                if is_zero(m0) {
                    self.bnlog_real::<RANK, W, P>(1.0, -1.0, lip, &mut roots);
                } else if real_masses {
                    let (x, w) = (m1.re / difference.re, -difference.re / m1.re);
                    self.bnlog_real::<RANK, W, P>(x, w, lip, &mut roots);
                } else {
                    let x = crate::native::primitives::mul(&m1, &inverse(&difference));
                    self.bnlog::<RANK, W, P>(x, lip, &mut roots);
                }
                for (n, value) in finite.iter_mut().enumerate().take(RANK + 1) {
                    *value = scale(sub(lna, roots[n]), INVERSE[n]);
                }
            }
            // The zero-momentum formulas give B(0; m1, m0).
            reconstruct = !switch;
        } else if is_zero(m0) {
            if m1.im == 0.0 && m1.re == momentum {
                // On-shell: ((n+1) lna + 2)/(n+1)^2.
                let lna = neg(mass_log(m1, mu_squared));
                for (n, value) in finite.iter_mut().enumerate().take(RANK + 1) {
                    let shifted = complex(lna.re * COUNTS[n] + 2.0, lna.im * COUNTS[n]);
                    *value = scale(shifted, INVERSE_SQUARED[n]);
                }
            } else {
                let shifted_mass = complex(m1.re - momentum, m1.im);
                let lna = neg(mass_log(shifted_mass, mu_squared));
                // x1 = (p²-m1)/p², -1/x1 = p²/(m1-p²); the lip has the sign
                // of p². A massless m1 gives x1 = p²/p² = 1 exactly.
                let lip = momentum > 0.0;
                if shifted_mass.im == 0.0 {
                    let x = -shifted_mass.re / momentum;
                    let w = momentum / shifted_mass.re;
                    self.bnlog_real::<RANK, W, P>(x, w, lip, &mut roots);
                } else {
                    let x = complex(-shifted_mass.re / momentum, -shifted_mass.im / momentum);
                    self.bnlog::<RANK, W, P>(x, lip, &mut roots);
                }
                for (n, value) in finite.iter_mut().enumerate().take(RANK + 1) {
                    // (((lna - bnlog) (n+1)) + 1)/(n+1)^2.
                    let difference = sub(lna, roots[n]);
                    let shifted =
                        complex(difference.re * COUNTS[n] + 1.0, difference.im * COUNTS[n]);
                    *value = scale(shifted, INVERSE_SQUARED[n]);
                }
            }
            reconstruct = switch;
        } else {
            let lna = neg(mass_log(m0, mu_squared));
            // Roots of p² x² + (m1 - m0 - p²) x + m0 = 0 (solabc).
            let twice_a = momentum * 2.0;
            if real_masses {
                let linear = m1.re - m0.re - momentum;
                let discriminant = linear * linear - m0.re * (momentum * 4.0);
                if discriminant < 0.0 {
                    // Complex-conjugate roots: bnlog(x1) + bnlog(x2) is
                    // 2 Re bnlog(x1), where x1 = (-b + i sqrt(-D))/(2a).
                    let x = complex(-linear / twice_a, (-discriminant).sqrt() / twice_a);
                    self.bnlog::<RANK, W, P>(x, true, &mut roots);
                    for (n, value) in finite.iter_mut().enumerate().take(RANK + 1) {
                        let sum = complex(lna.re - 2.0 * roots[n].re, lna.im);
                        *value = scale(sum, INVERSE[n]);
                    }
                } else {
                    let root = discriminant.sqrt();
                    let plus = root - linear;
                    let minus = -(root + linear);
                    let twice_c = m0.re * 2.0;
                    // x1 is the +sqrt root (upper lip), x2 the -sqrt root.
                    let (x1, w1, x2, w2) = if plus.abs() >= minus.abs() {
                        (
                            plus / twice_a,
                            -twice_a / plus,
                            twice_c / plus,
                            -plus / twice_c,
                        )
                    } else {
                        (
                            twice_c / minus,
                            -minus / twice_c,
                            minus / twice_a,
                            -twice_a / minus,
                        )
                    };
                    let mut second = [ZERO; W];
                    self.bnlog_real::<RANK, W, P>(x1, w1, true, &mut roots);
                    self.bnlog_real::<RANK, W, P>(x2, w2, false, &mut second);
                    for (n, value) in finite.iter_mut().enumerate().take(RANK + 1) {
                        *value = scale(sub(sub(lna, roots[n]), second[n]), INVERSE[n]);
                    }
                }
            } else {
                let linear = complex(m1.re - m0.re - momentum, m1.im - m0.im);
                let four_a_c = scale(m0, momentum * 4.0);
                let discriminant = sqrt(&sub(mul(linear, linear), four_a_c));
                let plus = sub(discriminant, linear);
                let minus = neg(add(discriminant, linear));
                // As the generic code: complex masses keep its roots exactly.
                let twice_c = scale(m0, 2.0);
                let (x1, x2) = if size(plus) >= size(minus) {
                    (
                        scale(plus, 1.0 / twice_a),
                        crate::native::primitives::mul(&twice_c, &inverse(&plus)),
                    )
                } else {
                    (
                        crate::native::primitives::mul(&twice_c, &inverse(&minus)),
                        scale(minus, 1.0 / twice_a),
                    )
                };
                let mut second = [ZERO; W];
                self.bnlog::<RANK, W, P>(x1, true, &mut roots);
                self.bnlog::<RANK, W, P>(x2, false, &mut second);
                for (n, value) in finite.iter_mut().enumerate().take(RANK + 1) {
                    *value = scale(sub(sub(lna, roots[n]), second[n]), INVERSE[n]);
                }
            }
            reconstruct = switch;
        }
        // Odd ranks carry the sign (-1)^n.
        for value in finite.iter_mut().skip(1).step_by(2) {
            *value = neg(*value);
        }
        if reconstruct {
            // B(p; m0, m1) from B(p; m1, m0) (the q -> -q-p substitution):
            // B_n = (-1)^n sum_k C(n,k) B_k', with a plus sign at rank four.
            let (b0, b1, b11, b111) = (finite[0], finite[1], finite[2], finite[3]);
            finite[1] = neg(add(b0, b1));
            finite[2] = add(add(b11, scale(b1, 2.0)), b0);
            finite[3] = neg(add(add(b111, scale(add(b11, b1), 3.0)), b0));
            if RANK >= 4 {
                // finite[W - 1] is finite[4] here (W = 5).
                let b1111 = finite[W - 1];
                finite[W - 1] = add(
                    add(add(b1111, scale(add(b111, b1), 4.0)), scale(b11, 6.0)),
                    b0,
                );
            }
        }
        for (n, value) in finite.into_iter().enumerate().take(RANK + 1) {
            output[3 * n] = value;
            output[3 * n + 1] = complex(POLES[n], 0.0);
            output[3 * n + 2] = ZERO;
        }
    }

    /// `[A0, A00, A0000]`; see the parent's `tadpole`.
    #[inline]
    pub(super) fn tadpole(&self, input: &[C64; 2], output: &mut [C64; 9]) {
        *output = [ZERO; 9];
        let mass = input[0];
        if is_zero(mass) {
            return;
        }
        let a0 = sub(mass, mul(mass, mass_log(mass, input[1].re)));
        // tadpn: A00 = m/4 (A0 + m/2) and A0000 = m²/24 (A0 + 5m/6), with
        // poles m²/4 and m³/24.
        let first = scale(mass, 0.25);
        let second = scale(mul(first, mass), INVERSE[2] * 0.5);
        output[0] = a0;
        output[1] = mass;
        output[3] = mul(first, add(a0, scale(mass, 0.5)));
        output[4] = mul(first, mass);
        output[6] = mul(second, add(a0, scale(mass, 0.5 + INVERSE[2])));
        output[7] = mul(second, mass);
    }
}
