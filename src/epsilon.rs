//! Positive Laurent orders in the common 1/r_Gamma normalization.
//!
//! B0[eps] = zeta(2) + 1/2 int_0^1 log^2((m0*(1-x)+m1*x-p*x*(1-x)-i0)/mu2) dx.
//! Integrating the two root logarithms gives logarithms and dilogarithms.
//! Endpoint, linear and coincident-root limits are evaluated separately.
use super::*;

fn log_with_sign(z: &Atom, sign: &Atom) -> Atom {
    let real_sign: Atom = (sign + sign.conj()) / 2;
    if_nonzero_else(
        &(&real_sign - Symbol::ABS.call((&real_sign,))),
        physical_log(z),
        upper_log(z),
    )
}

fn dilog_with_sign(z: &Atom, sign: &Atom) -> Atom {
    if z.is_zero() || z == &Atom::num(1) {
        return dilog_atom(z);
    }
    let one_minus: Atom = 1 - z;
    let below_cut: Atom = &one_minus - Symbol::ABS.call((&one_minus,));
    let continued: Atom = Symbol::PI.to_atom().pow(2) / 6
        - dilog_atom(&one_minus)
        - z.log() * log_with_sign(&one_minus, &(-sign));
    if_nonzero_else(
        &off_real_axis(z),
        dilog_atom(z),
        if_nonzero_else(&below_cut, continued, dilog_atom(z)),
    )
}

fn log_moments(root: &Atom, log: &Atom) -> (Atom, Atom) {
    ((1 - root) * log - 1, (1 - root) * (log * log - 2 * log) + 2)
}

// Integral log(1-x/s)/(x-r) dx from zero to one. The log correction
// continues Li2 along the straight parameter path, including a cut crossing.
fn cross_log(p: &Atom, r: &Atom, s: &Atom, lr: &Atom, ls: &Atom) -> Atom {
    let difference: Atom = r - s;
    let z0: Atom = r / &difference;
    let z1: Atom = (r - 1) / &difference;
    // Give both masses a common -i*eta. These are Im(dz/deta) on the
    // real axis and therefore determine the otherwise ambiguous cut lips.
    let denominator: Atom = p * difference.pow(3);
    let sign0: Atom = (-r - s) / &denominator;
    let sign1: Atom = (2 - r - s) / &denominator;
    let l0: Atom = log_with_sign(&(1 - &z0), &(-&sign0));
    let l1: Atom = log_with_sign(&(1 - &z1), &(-&sign1));
    -&l0 * lr + (ls - &l1 + &l0) * log_with_sign(&z1, &sign1) - dilog_with_sign(&z1, &sign1)
        + dilog_with_sign(&z0, &sign0)
}

fn weighted_log_square(mass: &Atom, scale: &Atom) -> Atom {
    if mass.is_zero() {
        return Atom::Zero;
    }
    let l: Atom = physical_log(&(mass / scale));
    if_nonzero(mass, mass * (&l * &l - 2 * &l + 2))
}

fn one_massless(p: &Atom, mass: &Atom, scale: &Atom) -> Atom {
    let shifted: Atom = mass - p;
    let lmass: Atom = physical_log(&(mass / scale));
    let on_shell: Atom = &lmass * &lmass - 4 * &lmass + 8;
    if shifted.is_zero() {
        return on_shell;
    }
    let root: Atom = -&shifted / p;
    let l0: Atom = physical_log(&(&shifted / scale));
    let l: Atom = log_with_sign(&(1 - Atom::num(1) / &root), p);
    let (j, k) = log_moments(&root, &l);
    let mixed: Atom =
        2 + (&root - 1) * &l - &root * dilog_with_sign(&(Atom::num(1) / &root), &(-p));
    let ordinary: Atom = &l0 * &l0 + 2 * &l0 * (j - 1) + k + 2 + 2 * mixed;
    if_nonzero_else(&shifted, ordinary, on_shell)
}

fn massive(p: &Atom, a: &Atom, b: &Atom, scale: &Atom) -> Atom {
    let linear: Atom = b - a - p;
    let discriminant: Atom = &linear * &linear - 4 * p * a;
    let delta: Atom = discriminant.sqrt();
    let r: Atom = (-&linear + &delta) / (2 * p);
    let s: Atom = (-&linear - &delta) / (2 * p);
    let l0: Atom = physical_log(&(a / scale));
    let l1: Atom = physical_log(&(b / scale));
    // At a double root, integrate log^2(p*(x-r)^2) directly. This also
    // covers an integrable real threshold inside the parameter interval.
    let coincident: Atom = (1 - &r) * (&l1 * &l1 - 4 * &l1 + 8) + &r * (&l0 * &l0 - 4 * &l0 + 8);
    if discriminant.is_zero() {
        return coincident;
    }
    let lr: Atom = upper_log(&(1 - Atom::num(1) / &r));
    let ls: Atom = physical_log(&(1 - Atom::num(1) / &s));
    let (jr, kr) = log_moments(&r, &lr);
    let (js, ks) = log_moments(&s, &ls);
    let mixed: Atom = &lr * &ls
        - &jr
        - &js
        - &r * cross_log(p, &r, &s, &lr, &ls)
        - &s * cross_log(p, &s, &r, &ls, &lr);
    let ordinary: Atom = &l0 * &l0 + 2 * &l0 * (jr + js) + kr + ks + 2 * mixed;
    if_nonzero_else(&discriminant, ordinary, coincident)
}

/// B0's O(epsilon) coefficient, with squared arguments and scale.
///
/// The normalization divides by r_Gamma, as for C0/D0. External p is
/// real, masses have non-positive imaginary parts, and mu_squared is
/// positive real. The completely scaleless bubble is zero.
pub fn b0_epsilon(p: &Atom, a: &Atom, b: &Atom, mu_squared: &Atom) -> Atom {
    fn zero_momentum(a: &Atom, b: &Atom, scale: &Atom) -> Atom {
        if a.is_zero() && b.is_zero() {
            return Atom::Zero;
        }
        let zeta: Atom = Symbol::PI.to_atom().pow(2) / 6;
        if a.is_zero() {
            return if_nonzero(b, &zeta + weighted_log_square(b, scale) / (2 * b));
        }
        if b.is_zero() {
            return if_nonzero(a, &zeta + weighted_log_square(a, scale) / (2 * a));
        }
        let difference: Atom = b - a;
        let equal: Atom = &zeta + physical_log(&(a / scale)).pow(2) / 2;
        if difference.is_zero() {
            return equal;
        }
        let unequal: Atom = &zeta
            + (weighted_log_square(b, scale) - weighted_log_square(a, scale)) / (2 * &difference);
        if_nonzero_else(
            a,
            if_nonzero_else(&difference, unequal.clone(), equal),
            if_nonzero_else(b, unequal, Atom::Zero),
        )
    }
    if p.is_zero() {
        return zero_momentum(a, b, mu_squared);
    }
    let zeta: Atom = Symbol::PI.to_atom().pow(2) / 6;
    let l: Atom = physical_log(&(-p / mu_squared));
    let massless: Atom = 4 - 2 * &l + &l * &l / 2;
    if a.is_zero() && b.is_zero() {
        return if_nonzero(p, massless);
    }
    let nonzero_p: Atom = if a.is_zero() {
        if_nonzero_else(b, &zeta + one_massless(p, b, mu_squared) / 2, massless)
    } else if b.is_zero() {
        if_nonzero_else(a, &zeta + one_massless(p, a, mu_squared) / 2, massless)
    } else {
        if_nonzero_else(
            a,
            if_nonzero_else(
                b,
                &zeta + massive(p, a, b, mu_squared) / 2,
                &zeta + one_massless(p, a, mu_squared) / 2,
            ),
            if_nonzero_else(b, &zeta + one_massless(p, b, mu_squared) / 2, massless),
        )
    };
    if_nonzero_else(p, nonzero_p, zero_momentum(a, b, mu_squared))
}
