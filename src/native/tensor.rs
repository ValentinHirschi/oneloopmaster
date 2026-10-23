//! Direct Rust bubble and tadpole tensor coefficients.
//!
//! This is a port of OneLOop's `bub0`/`bub1`/`bub11`/`bub111`/`bub1111`,
//! `bnlog`, `solabc` and `tadpn` (A. van Hameren, `avh_olo.f90`, after Denner and
//! Dittmaier, Nucl. Phys. B734 (2006) 62) to the generic [`NativeFloat`]
//! domains. It keeps OneLOop's algorithm: the mass ordering by `|Re m|+|Im m|`,
//! the root-pair representation `lna - bnlog(n, x1) - bnlog(n, x2)` for every
//! rank, the symmetric reconstruction after a mass swap, and the logarithmic
//! series of `bnlog` for large roots. Its deliberate differences follow this
//! crate's native backend:
//!
//! - Degenerate sectors (zero momentum, a zero mass, equal masses, an on-shell
//!   massless bubble, a unit root) use exact zero tests, not `10*EPSN` relative
//!   tolerances. Nearby points use the regular formulas, which stay finite.
//! - Real-axis lips are explicit, exactly where OneLOop's infinitesimal
//!   `IEPS = i*EPSN^2` regulator selects them; there is no finite regulator and
//!   near-real roots are not snapped to the axis.
//! - `y = log(1-1/x)` is formed with `log1p` and shared by all ranks of a
//!   root. The series is used for `|y|_1 < 1/2` with a precision-dependent
//!   term count, which bounds the direct formula's cancellation (about
//!   `|x|^(n+1)`) by roughly 1.6 digits at rank three (2 at rank four).
//!   OneLOop switches near `|y| = 1/(8(n+1))`, where the rank-three direct
//!   formula loses about six.
//! - The reconstruction after a mass swap is `B_n = (-1)^n sum_k C(n,k) B_k'`
//!   at every rank. OneLOop's `bub1111` (3.4 and 3.6) negates the rank-four
//!   sum in each of its three swap branches, so its `B1111` has the wrong sign
//!   whenever the swap is taken; this port does not.
//!
//! Inputs are `p²` (real), `m0²`, `m1²` with `Im <= 0`, and a positive real
//! `mu²` (the scale is passed squared, unlike OneLOop's `rmu`). Coefficients are
//! `[finite, 1/epsilon, 1/epsilon²]` in OneLOop's normalization, the same as the
//! scalar masters of this crate.
//!
//! `f64` evaluators use the binary64 specialization of the same algorithm in
//! [`binary64`] (real arithmetic for real kinematics, one `bnlog` for
//! complex-conjugate roots, Estrin's scheme for the series); `DoubleFloat`
//! and `Float` use the generic code below.
use super::primitives::{C, NativeFloat, inverse, log, mul, sqrt};
use symbolica::domains::{
    float::{Float, Real, RealLike, SingleFloat},
    integer::Integer,
    rational::Rational,
};

mod binary64;

/// Highest supported longitudinal bubble rank (`B1111`).
pub(crate) const MAX_BUBBLE_RANK: usize = 4;

/// Upper bound of the L1 norm of `y = log(1-1/x)` for the series of `bnlog`.
const SERIES_RADIUS: f64 = 0.5;

/// Truncated series of `bnlog(n, x)` in `y = log(1-1/x)`:
/// `bnlog(n,x) = -(n+1) y (x y)^(n+1) sum_k c_k y^k`, where
/// `c_k = (m-1)/m! sum_j (-1)^j binom(n,j) (j+1)^(m-2)` with `m = n+2+k`
/// (OneLOop's `update_coeff`, `coeff(m, n)`).
#[derive(Clone)]
struct BnlogSeries<T> {
    coefficients: Vec<T>,
    /// `radii[j]`: with `j+1` coefficients the first omitted term is below
    /// `2^-(bits+4) |c_0|` whenever `|y|_1 <= radii[j]`.
    radii: Vec<f64>,
}

impl<T: NativeFloat> BnlogSeries<T> {
    fn new(rank: usize, bits: u32, prototype: &T) -> Result<Self, String> {
        let target = -(f64::from(bits) + 4.0) * std::f64::consts::LN_2;
        let log_magnitude = |value: &Rational| {
            let float = Float::new(64).from_rational(value);
            float.norm().log().to_f64()
        };
        let mut exact = Vec::new();
        let mut logs = Vec::new();
        let mut radii = Vec::new();
        for m in (rank as u64 + 2).. {
            if m > 20_000 {
                return Err("bnlog series exceeds its coefficient budget".into());
            }
            let mut sum = Integer::zero();
            for j in 0..=rank as i64 {
                let term = Integer::binom(rank as i64, j) * Integer::from(j + 1).pow(m - 2);
                if j % 2 == 0 {
                    sum += term;
                } else {
                    sum -= term;
                }
            }
            let coefficient = Rational::new(
                Integer::from(m as i64 - 1) * sum,
                Integer::factorial(u32::try_from(m).map_err(|_| "series order overflow")?),
            );
            logs.push(log_magnitude(&coefficient));
            exact.push(coefficient);
            let count = exact.len() - 1;
            if count > 0 {
                // With `count` coefficients, c_count is the first omitted term.
                let radius = ((logs[0] - logs[count] + target) / count as f64).exp();
                radii.push(radius);
                if radius >= SERIES_RADIUS * 1.0625 {
                    break;
                }
            }
        }
        exact.pop();
        Ok(Self {
            coefficients: exact.iter().map(|c| prototype.from_rational(c)).collect(),
            radii,
        })
    }

    #[inline]
    fn evaluate(&self, y: &C<T>, size: f64) -> C<T> {
        let count = self
            .radii
            .iter()
            .position(|&radius| size <= radius)
            .map_or(self.coefficients.len(), |index| index + 1);
        let mut value = C::new(self.coefficients[count - 1].clone(), y.re.zero());
        for coefficient in self.coefficients[..count - 1].iter().rev() {
            value = mul(&value, y);
            value.re += coefficient;
        }
        value
    }
}

/// Precision-specific constants of the tensor kernels.
#[derive(Clone)]
pub(crate) struct TensorKernel<T> {
    zero: T,
    one: T,
    two: T,
    half: T,
    pi: T,
    /// `k+1` and `1/(k+1)` for `k = 0..=4`.
    counts: [T; MAX_BUBBLE_RANK + 1],
    inverse: [T; MAX_BUBBLE_RANK + 1],
    /// `1/(k+1)^2` for `k = 0..=4`.
    inverse_squared: [T; MAX_BUBBLE_RANK + 1],
    /// Harmonic numbers `H_(n+1)`: `bnlog(n, 1) = -H_(n+1)`.
    harmonic: [T; MAX_BUBBLE_RANK + 1],
    series: [BnlogSeries<T>; MAX_BUBBLE_RANK + 1],
}

impl<T: NativeFloat> TensorKernel<T> {
    pub(crate) fn new(bits: u32) -> Result<Self, String> {
        let prototype = T::zero_at(bits);
        let rational = |num: i64, den: i64| prototype.from_rational(&Rational::new(num, den));
        Ok(Self {
            zero: prototype.zero(),
            one: prototype.one(),
            two: prototype.from_i64(2),
            half: rational(1, 2),
            pi: prototype.pi(),
            counts: std::array::from_fn(|k| prototype.from_i64(k as i64 + 1)),
            inverse: std::array::from_fn(|k| rational(1, k as i64 + 1)),
            inverse_squared: std::array::from_fn(|k| rational(1, (k as i64 + 1).pow(2))),
            harmonic: [
                rational(1, 1),
                rational(3, 2),
                rational(11, 6),
                rational(25, 12),
                rational(137, 60),
            ],
            series: [
                BnlogSeries::new(0, bits, &prototype)?,
                BnlogSeries::new(1, bits, &prototype)?,
                BnlogSeries::new(2, bits, &prototype)?,
                BnlogSeries::new(3, bits, &prototype)?,
                BnlogSeries::new(4, bits, &prototype)?,
            ],
        })
    }

    #[inline]
    fn complex(&self, re: T) -> C<T> {
        C::new(re, self.zero.clone())
    }

    #[inline]
    fn scale(value: &C<T>, factor: &T) -> C<T> {
        C::new(value.re.clone() * factor, value.im.clone() * factor)
    }

    /// `log(m/mu²)` on the Feynman lip: real negative ratios take `-i pi`
    /// (OneLOop's `logc(qonv(m/rmu2,-1))`).
    #[inline]
    fn mass_log(&self, mass: &C<T>, mu_squared: &T) -> C<T> {
        if mass.im.is_fully_zero() {
            let ratio = mass.re.clone() / mu_squared;
            if ratio < self.zero {
                C::new((-ratio).log(), -self.pi.clone())
            } else {
                C::new(ratio.log(), self.zero.clone())
            }
        } else {
            log(&C::new(
                mass.re.clone() / mu_squared,
                mass.im.clone() / mu_squared,
            ))
        }
    }

    /// Principal `log(1+w)`, retaining the relative precision of small `w`.
    #[inline]
    fn log1p(&self, w: &C<T>) -> C<T> {
        let real = self.one.clone() + &w.re;
        let modulus = if w.re.norm() < self.half && w.im.norm() < self.half {
            // |1+w|^2 - 1 = w.re (2 + w.re) + w.im^2.
            let shifted = w.re.clone() * (self.two.clone() + &w.re) + w.im.clone() * &w.im;
            shifted.log1p() * &self.half
        } else {
            NativeFloat::hypot(&real, &w.im).log()
        };
        C::new(modulus, w.im.atan2(&real))
    }

    /// `bnlog(n, x) = (n+1) int_0^1 t^n log(1 - t/x) dt` for `n = 0..=rank`.
    ///
    /// For real `0 < x < 1`, `1 - 1/x` lies on the negative axis and
    /// `upper_lip` selects `+i pi` (OneLOop's `sgn` of `Im x`).
    fn bnlog(&self, x: &C<T>, upper_lip: bool, rank: usize, out: &mut [C<T>; MAX_BUBBLE_RANK + 1]) {
        let real = x.im.is_fully_zero();
        if real && x.re == self.one {
            for (value, harmonic) in out.iter_mut().zip(&self.harmonic).take(rank + 1) {
                *value = self.complex(-harmonic.clone());
            }
            return;
        }
        let y = if real && x.re > self.zero && x.re < self.one {
            let magnitude = ((self.one.clone() - &x.re) / &x.re).log();
            C::new(
                magnitude,
                if upper_lip {
                    self.pi.clone()
                } else {
                    -self.pi.clone()
                },
            )
        } else {
            let reciprocal = inverse(x);
            self.log1p(&C::new(-reciprocal.re, -reciprocal.im))
        };
        if y.re.is_fully_zero() && y.im.is_fully_zero() {
            // 1/x underflowed: the x -> infinity limit.
            for value in out.iter_mut().take(rank + 1) {
                *value = self.complex(self.zero.clone());
            }
            return;
        }
        let size = y.re.norm().to_f64() + y.im.norm().to_f64();
        if size < SERIES_RADIUS {
            let xy = mul(x, &y);
            let mut power = xy.clone();
            for (n, value) in out.iter_mut().enumerate().take(rank + 1) {
                if n > 0 {
                    power = mul(&power, &xy);
                }
                let series = self.series[n].evaluate(&y, size);
                let product = mul(&mul(&series, &y), &power);
                *value = Self::scale(&product, &-self.counts[n].clone());
            }
        } else {
            // (1 - x^(n+1)) y - sum_(k=0..n) x^k/(n+1-k).
            let mut powers: [C<T>; MAX_BUBBLE_RANK + 2] =
                std::array::from_fn(|_| self.complex(self.one.clone()));
            for k in 1..=rank + 1 {
                powers[k] = mul(&powers[k - 1], x);
            }
            for (n, value) in out.iter_mut().enumerate().take(rank + 1) {
                let mut polynomial = self.complex(self.zero.clone());
                for (k, power) in powers.iter().enumerate().take(n + 1) {
                    polynomial += &Self::scale(power, &self.inverse[n - k]);
                }
                let one_minus = C::new(
                    self.one.clone() - &powers[n + 1].re,
                    -powers[n + 1].im.clone(),
                );
                *value = mul(&one_minus, &y) - &polynomial;
            }
        }
    }

    /// `[B0, B1, B11, B111, B1111]` through `rank`, as `[finite, 1/eps,
    /// 1/eps²]` triples. `input` is `[p², m0², m1², mu²]`, already validated.
    pub(crate) fn bubble(&self, input: &[C<T>], rank: usize, output: &mut [C<T>]) {
        debug_assert!(rank <= MAX_BUBBLE_RANK && output.len() == 3 * (rank + 1));
        let momentum = &input[0].re;
        let mu_squared = &input[3].re;
        let is_zero = |z: &C<T>| z.re.is_fully_zero() && z.im.is_fully_zero();
        for value in output.iter_mut() {
            *value = self.complex(self.zero.clone());
        }
        if momentum.is_fully_zero() && is_zero(&input[1]) && is_zero(&input[2]) {
            return;
        }
        let size = |z: &C<T>| z.re.norm() + z.im.norm();
        // OneLOop keeps the order when |m1| >= |m0| and otherwise computes
        // B(p; m1, m0), reconstructing the requested order at the end.
        let switch = size(&input[2]) < size(&input[1]);
        let (m0, m1) = if switch {
            (&input[2], &input[1])
        } else {
            (&input[1], &input[2])
        };

        let mut roots: [[C<T>; MAX_BUBBLE_RANK + 1]; 2] =
            std::array::from_fn(|_| std::array::from_fn(|_| self.complex(self.zero.clone())));
        let mut finite: [C<T>; MAX_BUBBLE_RANK + 1] =
            std::array::from_fn(|_| self.complex(self.zero.clone()));
        let reconstruct;
        if momentum.is_fully_zero() {
            let lna = -self.mass_log(m1, mu_squared);
            if m0 == m1 {
                for (n, value) in finite.iter_mut().enumerate().take(rank + 1) {
                    *value = Self::scale(&lna, &self.inverse[n]);
                }
            } else {
                // x1 = m1/(m1-m0); OneLOop's -i*EPSN^2*|m1| in the numerator
                // puts real roots on the lip opposite to Re(m1-m0). A massless
                // m0 gives the unit root exactly: m1 * (1/m1) need not be 1.
                let difference = m1.clone() - m0;
                let x1 = if is_zero(m0) {
                    self.complex(self.one.clone())
                } else {
                    mul(m1, &inverse(&difference))
                };
                self.bnlog(&x1, difference.re < self.zero, rank, &mut roots[0]);
                for (n, value) in finite.iter_mut().enumerate().take(rank + 1) {
                    *value = Self::scale(&(lna.clone() - &roots[0][n]), &self.inverse[n]);
                }
            }
            // The zero-momentum formulas give B(0; m1, m0).
            reconstruct = !switch;
        } else if is_zero(m0) {
            if m1.im.is_fully_zero() && m1.re == *momentum {
                // On-shell: ((n+1) lna + 2)/(n+1)^2.
                let lna = -self.mass_log(m1, mu_squared);
                for (n, value) in finite.iter_mut().enumerate().take(rank + 1) {
                    let shifted = C::new(
                        lna.re.clone() * &self.counts[n] + &self.two,
                        lna.im.clone() * &self.counts[n],
                    );
                    *value = Self::scale(&shifted, &self.inverse_squared[n]);
                }
            } else {
                let shifted_mass = C::new(m1.re.clone() - momentum, m1.im.clone());
                let lna = -self.mass_log(&shifted_mass, mu_squared);
                // x1 = (p²-m1)/p²; the regulator +i*EPSN^2*|m1|/p² selects
                // the lip with the sign of p². Divide rather than scale by
                // 1/p², so that a massless m1 gives the unit root exactly.
                let x1 = C::new(
                    -shifted_mass.re.clone() / momentum,
                    -shifted_mass.im.clone() / momentum,
                );
                self.bnlog(&x1, *momentum > self.zero, rank, &mut roots[0]);
                for (n, value) in finite.iter_mut().enumerate().take(rank + 1) {
                    // (((lna - bnlog) (n+1)) + 1)/(n+1)^2.
                    let difference = lna.clone() - &roots[0][n];
                    let shifted = C::new(
                        difference.re * &self.counts[n] + &self.one,
                        difference.im * &self.counts[n],
                    );
                    *value = Self::scale(&shifted, &self.inverse_squared[n]);
                }
            }
            reconstruct = switch;
        } else {
            let lna = -self.mass_log(m0, mu_squared);
            // Roots of p² x² + (m1 - m0 - p²) x + m0 = 0 (solabc).
            let linear = C::new(m1.re.clone() - &m0.re - momentum, m1.im.clone() - &m0.im);
            let four_a_c = Self::scale(m0, &(momentum.clone() * self.two.clone() * &self.two));
            let discriminant = sqrt(&(mul(&linear, &linear) - &four_a_c));
            let plus = discriminant.clone() - &linear;
            let minus = -(discriminant + &linear);
            let twice_cc = Self::scale(m0, &self.two);
            let twice_a = momentum.clone() * &self.two;
            let (x1, x2) = if size(&plus) >= size(&minus) {
                let x1 = Self::scale(&plus, &(self.one.clone() / &twice_a));
                (x1, mul(&twice_cc, &inverse(&plus)))
            } else {
                let x2 = Self::scale(&minus, &(self.one.clone() / &twice_a));
                (mul(&twice_cc, &inverse(&minus)), x2)
            };
            // The +sqrt root takes the upper lip, the -sqrt root the lower one.
            let [first, second] = &mut roots;
            self.bnlog(&x1, true, rank, first);
            self.bnlog(&x2, false, rank, second);
            for (n, value) in finite.iter_mut().enumerate().take(rank + 1) {
                let sum = lna.clone() - &roots[0][n] - &roots[1][n];
                *value = Self::scale(&sum, &self.inverse[n]);
            }
            reconstruct = switch;
        }
        // Odd ranks carry the sign (-1)^n.
        for value in finite.iter_mut().skip(1).step_by(2) {
            *value = -value.clone();
        }
        if reconstruct {
            // B(p; m0, m1) from B(p; m1, m0) (the q -> -q-p substitution):
            // B_n = (-1)^n sum_k C(n,k) B_k', with a plus sign at rank four.
            let [b0, b1, b11, b111, b1111] = finite.clone();
            finite[1] = -(b0.clone() + &b1);
            if rank >= 2 {
                finite[2] = b11.clone() + &Self::scale(&b1, &self.two) + &b0;
            }
            if rank >= 3 {
                finite[3] =
                    -(b111.clone() + &Self::scale(&(b11.clone() + &b1), &self.counts[2]) + &b0);
            }
            if rank >= 4 {
                finite[4] = b1111
                    + &Self::scale(&(b111 + &b1), &self.counts[3])
                    + &Self::scale(&b11, &(self.counts[2].clone() * &self.two))
                    + &b0;
            }
        }
        let poles = [
            self.one.clone(),
            -self.half.clone(),
            self.inverse[2].clone(),
            -self.inverse[3].clone(),
            self.inverse[4].clone(),
        ];
        for (n, (value, pole)) in finite.into_iter().zip(poles).enumerate().take(rank + 1) {
            output[3 * n] = value;
            output[3 * n + 1] = self.complex(pole);
        }
    }

    /// `[A0, A00, A0000]` (OneLOop's `tadpn`), from `[m², mu²]`.
    pub(crate) fn tadpole(&self, input: &[C<T>], output: &mut [C<T>]) {
        debug_assert!(output.len() == 9);
        for value in output.iter_mut() {
            *value = self.complex(self.zero.clone());
        }
        let mass = &input[0];
        if mass.re.is_fully_zero() && mass.im.is_fully_zero() {
            return;
        }
        let a0 = mass.clone() - &mul(mass, &self.mass_log(mass, &input[1].re));
        // tadpn: A00 = m/4 (A0 + m/2) and A0000 = m²/24 (A0 + 5m/6), with
        // poles m²/4 and m³/24.
        let first = Self::scale(mass, &self.inverse[3]);
        let second = Self::scale(&mul(&first, mass), &(self.inverse[2].clone() * &self.half));
        let first_shift = self.half.clone();
        let second_shift = self.half.clone() + &self.inverse[2];
        let a00 = mul(&first, &(a0.clone() + &Self::scale(mass, &first_shift)));
        let a0000 = mul(&second, &(a0.clone() + &Self::scale(mass, &second_shift)));
        output[0] = a0;
        output[1] = mass.clone();
        output[3] = a00;
        output[4] = mul(&first, mass);
        output[6] = a0000;
        output[7] = mul(&second, mass);
    }
}

/// Tensor-coefficient families evaluated together, like OneLOop's `olo_an`
/// and `olo_b11`/`olo_bn`.
///
/// Each row of input is the family's [`arity`](Self::arity) arguments with the
/// squared scale last: `[m², mu²]` for [`An`](Self::An) and `[p², m0², m1², mu²]`
/// for the bubbles, where the bubble denominators are `(q²-m0²)((q+p)²-m1²)`.
/// Each row of output is [`coefficients`](Self::coefficients) consecutive
/// `[finite, 1/epsilon, 1/epsilon²]` triples, in the order listed per variant.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum TensorIntegral {
    /// `[A0, A00, A0000]` (OneLOop's `olo_an` at rank four).
    An,
    /// `[B0, B1, B11]`: the coefficients of `olo_b11` that a rank-`n`
    /// reduction such as Ninja's consumes (`B00` is not computed).
    B11,
    /// `[B0, B1, B11, B111]`, for reductions with numerator rank `n+1`.
    B111,
    /// `[B0, B1, B11, B111, B1111]`, for reductions with numerator rank
    /// `n+2`. `B1111 = 1/(5 eps) - int_0^1 x^4 log(D(x)/mu²) dx`; its
    /// first four coefficients are those of [`B111`](Self::B111), bitwise.
    B1111,
}

impl TensorIntegral {
    /// Family name.
    pub fn name(self) -> &'static str {
        match self {
            Self::An => "An",
            Self::B11 => "B11",
            Self::B111 => "B111",
            Self::B1111 => "B1111",
        }
    }

    /// Number of numeric arguments per row, including `mu²`.
    pub fn arity(self) -> usize {
        match self {
            Self::An => 2,
            Self::B11 | Self::B111 | Self::B1111 => 4,
        }
    }

    /// Number of tensor coefficients per row; each is a Laurent triple.
    pub fn coefficients(self) -> usize {
        match self {
            Self::An | Self::B11 => 3,
            Self::B111 => 4,
            Self::B1111 => 5,
        }
    }

    /// Complex outputs per row, `3 * coefficients()`.
    pub fn outputs(self) -> usize {
        3 * self.coefficients()
    }

    /// Names of the coefficients, in output order.
    pub fn coefficient_names(self) -> &'static [&'static str] {
        match self {
            Self::An => &["A0", "A00", "A0000"],
            Self::B11 => &["B0", "B1", "B11"],
            Self::B111 => &["B0", "B1", "B11", "B111"],
            Self::B1111 => &["B0", "B1", "B11", "B111", "B1111"],
        }
    }

    fn scalar_domain(self) -> crate::ScalarIntegral {
        match self {
            Self::An => crate::ScalarIntegral::A0,
            Self::B11 | Self::B111 | Self::B1111 => crate::ScalarIntegral::B0,
        }
    }
}

/// A reusable direct Rust evaluator of tensor coefficients.
///
/// `T` may be `f64`, Symbolica's `DoubleFloat`, or arbitrary-precision `Float`.
/// The arithmetic is a hand-written port of OneLOop's tensor routines (see the
/// [module documentation](self)); it needs neither generated code nor a
/// runtime evaluator, so it is available without `generated-evaluators`.
/// Inputs and the domain contract are those of [`crate::NativeEvaluator`] for
/// A0 (`An`) and B0 (bubbles): finite components, real `p²`, `Im m² <= 0`, and
/// a positive real `mu²` passed **squared**. Shape or domain errors leave the
/// output unchanged. Batches are flat row-major arrays; every row is validated
/// before any output is written, and no worker threads are spawned.
#[derive(Clone)]
pub struct NativeTensorEvaluator<T: NativeFloat> {
    family: TensorIntegral,
    poles: crate::PoleSelection,
    bits: u32,
    kernel: TensorKernel<T>,
    /// The binary64 kernels, present exactly for `T = f64`.
    binary64: Option<binary64::Binary64>,
    arguments: Vec<C<T>>,
}

impl<T: NativeFloat> NativeTensorEvaluator<T> {
    /// Construct at the numeric type's default working precision.
    pub fn new(family: TensorIntegral) -> Result<Self, String> {
        Self::with_binary_precision(family, T::DEFAULT_BITS)
    }

    /// Construct with an explicit bit count. Fixed-precision types require
    /// their native bit count; use `Float` for arbitrary working precision.
    pub fn with_binary_precision(family: TensorIntegral, bits: u32) -> Result<Self, String> {
        if bits == 0 {
            return Err("binary precision must be positive".into());
        }
        if let Some(fixed) = T::FIXED_BITS
            && bits != fixed
        {
            return Err(format!(
                "this numeric type has fixed precision {fixed} bits"
            ));
        }
        Ok(Self {
            family,
            poles: crate::PoleSelection::BOTH,
            bits,
            kernel: TensorKernel::new(bits)?,
            binary64: Self::binary64_kernel()?,
            arguments: (0..family.arity())
                .map(|_| C::new(T::zero_at(bits), T::zero_at(bits)))
                .collect(),
        })
    }

    fn binary64_kernel() -> Result<Option<binary64::Binary64>, String> {
        if std::any::TypeId::of::<T>() == std::any::TypeId::of::<f64>() {
            binary64::Binary64::new().map(Some)
        } else {
            Ok(None)
        }
    }

    /// Select the binary64 kernels (the `f64` default) or the generic code
    /// for `f64`; other types always use the generic code. For validation only.
    #[doc(hidden)]
    pub fn set_binary64(&mut self, enabled: bool) -> Result<(), String> {
        self.binary64 = if enabled {
            Self::binary64_kernel()?
        } else {
            None
        };
        Ok(())
    }

    /// Family represented by this evaluator.
    pub fn family(&self) -> TensorIntegral {
        self.family
    }

    /// Working precision in bits (not decimal digits).
    pub fn binary_precision(&self) -> u32 {
        self.bits
    }

    /// Retain the chosen UV/IR poles after the complete tensor calculation.
    /// All finite coefficients, including rational terms, are unchanged.
    pub fn with_pole_selection(mut self, selection: crate::PoleSelection) -> Self {
        self.set_pole_selection(selection);
        self
    }

    /// Change the selection for future rows and batches on this evaluator.
    pub fn set_pole_selection(&mut self, selection: crate::PoleSelection) {
        self.poles = selection;
    }

    pub fn pole_selection(&self) -> crate::PoleSelection {
        self.poles
    }

    /// Evaluate one row into `family().outputs()` values.
    pub fn evaluate(&mut self, input: &[C<T>], output: &mut [C<T>]) -> Result<(), String> {
        self.validate_shape(input.len(), output.len(), 1)?;
        super::validate_point(self.family.scalar_domain(), input)?;
        self.evaluate_row(input, output);
        Ok(())
    }

    /// Evaluate `rows` flat row-major points. Every row is validated before
    /// any output is written; empty batches are valid.
    pub fn evaluate_batch(
        &mut self,
        input: &[C<T>],
        output: &mut [C<T>],
        rows: usize,
    ) -> Result<(), String> {
        self.validate_shape(input.len(), output.len(), rows)?;
        let arity = self.family.arity();
        for point in input.chunks_exact(arity) {
            super::validate_point(self.family.scalar_domain(), point)?;
        }
        for (point, result) in input
            .chunks_exact(arity)
            .zip(output.chunks_exact_mut(self.family.outputs()))
        {
            self.evaluate_row(point, result);
        }
        Ok(())
    }

    #[inline]
    fn evaluate_row(&mut self, input: &[C<T>], output: &mut [C<T>]) {
        if let Some(kernel) = &self.binary64
            && binary64_row(kernel, self.family, input, output)
        {
            self.select_poles(output);
            return;
        }
        let input = if T::FIXED_BITS.is_some() {
            input
        } else {
            for (working, value) in self.arguments.iter_mut().zip(input) {
                working.re = value.re.at_precision(self.bits);
                working.im = value.im.at_precision(self.bits);
            }
            &self.arguments
        };
        match self.family {
            TensorIntegral::An => self.kernel.tadpole(input, output),
            TensorIntegral::B11 => self.kernel.bubble(input, 2, output),
            TensorIntegral::B111 => self.kernel.bubble(input, 3, output),
            TensorIntegral::B1111 => self.kernel.bubble(input, 4, output),
        }
        self.select_poles(output);
    }

    fn select_poles(&self, output: &mut [C<T>]) {
        if self.poles == crate::PoleSelection::BOTH
            || output
                .iter()
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return;
        }
        for (rank, triple) in output.chunks_exact_mut(3).enumerate() {
            let coefficient = if self.family == TensorIntegral::An {
                crate::PoleCoefficient::TadpoleTensor
            } else {
                crate::PoleCoefficient::BubbleMoment { rank }
            };
            self.poles
                .apply(coefficient, triple.try_into().expect("Laurent triple"))
                .expect("native tensor families stop at rank four");
        }
    }

    #[inline]
    fn validate_shape(&self, inputs: usize, outputs: usize, rows: usize) -> Result<(), String> {
        let expected_inputs = rows
            .checked_mul(self.family.arity())
            .ok_or("input size overflow")?;
        let expected_outputs = rows
            .checked_mul(self.family.outputs())
            .ok_or("output size overflow")?;
        if inputs != expected_inputs || outputs != expected_outputs {
            return Err(format!(
                "expected {expected_inputs} inputs and {expected_outputs} outputs for {rows} rows; received {inputs} and {outputs}"
            ));
        }
        Ok(())
    }
}

/// One validated row through the binary64 kernels; `false` (nothing written)
/// unless `T` is `f64` and the shapes are the family's.
#[inline]
fn binary64_row<T: NativeFloat>(
    kernel: &binary64::Binary64,
    family: TensorIntegral,
    input: &[C<T>],
    output: &mut [C<T>],
) -> bool {
    use std::any::Any;
    fn cast<'a, T: 'static, const N: usize, const M: usize>(
        input: &'a [C<T>],
        output: &'a mut [C<T>],
    ) -> Option<(&'a [C<f64>; N], &'a mut [C<f64>; M])> {
        let input: &[C<T>; N] = input.try_into().ok()?;
        let output: &mut [C<T>; M] = output.try_into().ok()?;
        Some((
            (input as &dyn Any).downcast_ref()?,
            (output as &mut dyn Any).downcast_mut()?,
        ))
    }
    match family {
        TensorIntegral::An => cast::<T, 2, 9>(input, output)
            .map(|(input, output)| kernel.tadpole(input, output))
            .is_some(),
        TensorIntegral::B11 => cast::<T, 4, 9>(input, output)
            .map(|(input, output)| kernel.bubble::<2, 4, 5>(input, output))
            .is_some(),
        TensorIntegral::B111 => cast::<T, 4, 12>(input, output)
            .map(|(input, output)| kernel.bubble::<3, 4, 5>(input, output))
            .is_some(),
        TensorIntegral::B1111 => cast::<T, 4, 15>(input, output)
            .map(|(input, output)| kernel.bubble_rank_four(input, output))
            .is_some(),
    }
}
