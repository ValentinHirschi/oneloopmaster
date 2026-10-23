//! Ultraviolet and infrared components of one-loop Laurent coefficients.
//!
//! A massless bubble at zero invariant is zero with a single dimensional
//! regulator, but has residues +1 (UV) and -1 (IR) if their origins are kept
//! separate. The longitudinal tensor moments have UV residues
//! `integral_0^1 (-x)^rank dx = (-1)^rank/(rank+1)`.
//!
//! Pole selection is an operation on a completed integral. It must not be
//! inserted into a recurrence that also uses pole coefficients to calculate
//! finite terms. The ordinary finite coefficient is always retained.

use crate::ScalarIntegral;
use symbolica::domains::float::FloatLike;

/// Which pole origins to retain. Both are enabled by default.
///
/// This is a per-evaluator choice; no process-global state is changed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PoleSelection {
    pub uv: bool,
    pub ir: bool,
}

impl Default for PoleSelection {
    fn default() -> Self {
        Self::BOTH
    }
}

/// The scalar or tensor coefficient whose poles are being selected.
///
/// The bubble convention is `(q^2-m0^2)((q+p)^2-m1^2)`. `BubbleMoment`
/// labels B0, B1, B11, B111 or B1111 by the number of longitudinal indices.
/// `TadpoleTensor` includes A0, A00 and A0000; each has only UV poles.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PoleCoefficient {
    Scalar(ScalarIntegral),
    TadpoleTensor,
    BubbleMoment { rank: usize },
}

impl PoleSelection {
    pub const BOTH: Self = Self { uv: true, ir: true };
    pub const UV: Self = Self {
        uv: true,
        ir: false,
    };
    pub const IR: Self = Self {
        uv: false,
        ir: true,
    };
    pub const NONE: Self = Self {
        uv: false,
        ir: false,
    };

    /// Select poles in `[finite, simple pole, double pole]` order.
    ///
    /// Call this after successfully evaluating the combined result, including
    /// any finite rational terms. It also serves expression-based or external
    /// evaluators, so the pole conventions have a single public implementation.
    /// The numerical type and its precision are preserved; fractions are formed
    /// from integers in that type, never from rounded binary64 constants.
    ///
    /// A rank above four is unsupported and leaves the complete array unchanged.
    /// The caller remains responsible for checking that the underlying integral
    /// succeeded: selection is not an input validator or a numerical rescue.
    pub fn apply<T: FloatLike>(
        self,
        coefficient: PoleCoefficient,
        laurent: &mut [T; 3],
    ) -> Result<(), String> {
        if matches!(coefficient, PoleCoefficient::BubbleMoment { rank } if rank > 4) {
            return Err("bubble pole selection supports ranks zero through four".into());
        }
        if self == Self::BOTH {
            return Ok(());
        }
        let uv = match coefficient {
            PoleCoefficient::Scalar(ScalarIntegral::A0) | PoleCoefficient::TadpoleTensor => {
                laurent[1].clone()
            }
            PoleCoefficient::Scalar(ScalarIntegral::B0) => laurent[1].one(),
            PoleCoefficient::BubbleMoment { rank } => {
                let sign = if rank % 2 == 0 { 1 } else { -1 };
                laurent[1].from_i64(sign) / laurent[1].from_i64((rank + 1) as i64)
            }
            PoleCoefficient::Scalar(
                ScalarIntegral::DB0 | ScalarIntegral::C0 | ScalarIntegral::D0,
            ) => laurent[1].zero(),
        };
        if !self.uv {
            laurent[1] = laurent[1].clone() - &uv;
        }
        if !self.ir {
            laurent[1] = if self.uv { uv } else { laurent[1].zero() };
            laurent[2] = laurent[2].zero();
        }
        Ok(())
    }
}
