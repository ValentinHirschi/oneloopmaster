//! Exact relabelling of soft boxes before numerical evaluation.
//!
//! The scalar box depends only on its set of propagators. The box expressions
//! (OneLOop's d0) find the soft limits on the four external invariants only: a
//! massless line i with two lines k on their mass shells, (q_i - q_k)^2 = m_k
//! exactly, across its legs. A relabelling of the same box can put such a pair
//! with a massive k on a channel (s12 or s23), where the finite formulas apply
//! and miss the pole. Every numerical D0 entry point therefore evaluates such a
//! point in the other labelling, of the two that exchange the channel pair,
//! whose channels are nonzero and hold no such pair, using exact tests on the
//! input values. Every other point keeps its labelling. Symbolic uses of the
//! D0 expressions (transparent definitions, inspection, Dual evaluators) do
//! not pass these entry points; callers evaluating them numerically can apply
//! [`regular_d0_labelling`] to each point first.
use crate::NativeFloat;
use symbolica::domains::float::Complex;

/// Slot of the invariant `(q_i - q_k)^2` in `[p1², p2², p3², p4², s12, s23]`.
const SLOT: [[usize; 4]; 4] = [[6, 0, 4, 3], [0, 6, 1, 5], [4, 1, 6, 2], [3, 5, 2, 6]];

/// Argument orders of the labellings (1, 3, 2, 4) and (1, 2, 4, 3): the
/// channel exchanges of `box_series`, with the masses and the scale.
const EXCHANGES: [[usize; 11]; 2] = [
    [4, 1, 5, 3, 0, 2, 6, 8, 7, 9, 10],
    [0, 5, 2, 4, 3, 1, 6, 7, 9, 8, 10],
];

fn is_zero<T: NativeFloat>(z: &Complex<T>) -> bool {
    z.re.is_fully_zero() && z.im.is_fully_zero()
}

/// `value - mass` vanishes exactly, as the expressions' zero tests decide.
fn on_shell<T: NativeFloat>(value: &Complex<T>, mass: &Complex<T>) -> bool {
    (value.re.clone() - &mass.re).is_fully_zero() && (value.im.clone() - &mass.im).is_fully_zero()
}

/// The argument order in which to evaluate a D0 point, if its labelling puts
/// a soft line's on-shell pair on a channel. `point` holds the eleven D0
/// arguments `[p1², p2², p3², p4², s12, s23, m1², m2², m3², m4², mu²]`.
pub(crate) fn soft_order<T: NativeFloat>(point: &[Complex<T>]) -> Option<&'static [usize; 11]> {
    if point.len() != 11 {
        return None;
    }
    let (p, rest) = point.split_at(6);
    let m = &rest[..4];
    // Only a massless line can be soft.
    if !m.iter().any(is_zero) {
        return None;
    }
    let pair = |i: usize, k: usize| &p[SLOT[i][k]];
    // False iff massless line i has its massive channel partner k and one of
    // its leg partners j, l on shell.
    let soft_free = |i: usize, k: usize, j: usize, l: usize| {
        !is_zero(&m[i])
            || is_zero(&m[k])
            || !on_shell(pair(i, k), &m[k])
            || (!on_shell(pair(i, j), &m[j]) && !on_shell(pair(i, l), &m[l]))
    };
    // No channel of the labelling (a, b, c, d) holds such a pair.
    let labelling_free = |[a, b, c, d]: [usize; 4]| {
        soft_free(a, c, b, d)
            && soft_free(c, a, b, d)
            && soft_free(b, d, a, c)
            && soft_free(d, b, a, c)
    };
    if labelling_free([0, 1, 2, 3]) {
        None
    } else if !is_zero(&p[0]) && !is_zero(&p[2]) && labelling_free([0, 2, 1, 3]) {
        Some(&EXCHANGES[0])
    } else if !is_zero(&p[3]) && !is_zero(&p[1]) && labelling_free([0, 1, 3, 2]) {
        Some(&EXCHANGES[1])
    } else {
        None
    }
}

/// Relabel one D0 point `[p1², p2², p3², p4², s12, s23, m1², m2², m3², m4²,
/// mu²]` in place, as every numerical D0 entry point of this crate does: when
/// a channel holds an exactly on-shell pair of a soft line, permute the point
/// into the labelling of the same box, of the two that exchange the channel
/// pair, whose channels are nonzero and hold no such pair. Returns whether it
/// permuted. Apply it to points evaluated with the symbolic D0 expressions
/// ([`crate::OneLoopExpressions::d0`], the transparent definitions), which do
/// not relabel. Points of any other length are left unchanged.
pub fn regular_d0_labelling<T: NativeFloat>(point: &mut [Complex<T>]) -> bool {
    match soft_order(point) {
        Some(order) => {
            let permuted: [Complex<T>; 11] = core::array::from_fn(|k| point[order[k]].clone());
            for (target, value) in point.iter_mut().zip(permuted) {
                *target = value;
            }
            true
        }
        None => false,
    }
}

/// `point` itself, or its relabelled copy in `storage` (see [`soft_order`]).
pub(crate) fn regular_point<'a, T: NativeFloat>(
    point: &'a [Complex<T>],
    storage: &'a mut Vec<Complex<T>>,
) -> &'a [Complex<T>] {
    match soft_order(point) {
        Some(order) => {
            storage.clear();
            storage.extend(order.iter().map(|&k| point[k].clone()));
            storage
        }
        None => point,
    }
}

/// Rows of a flat row-major D0 batch, relabelled into `storage` when any row
/// needs it (see [`soft_order`]).
pub(crate) fn regular_rows<'a, T: NativeFloat>(
    rows: &'a [Complex<T>],
    storage: &'a mut Vec<Complex<T>>,
) -> &'a [Complex<T>] {
    if rows.len() % 11 != 0 || !rows.chunks_exact(11).any(|row| soft_order(row).is_some()) {
        return rows;
    }
    storage.clear();
    for row in rows.chunks_exact(11) {
        match soft_order(row) {
            Some(order) => storage.extend(order.iter().map(|&k| row[k].clone())),
            None => storage.extend_from_slice(row),
        }
    }
    storage
}
