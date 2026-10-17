//! Tadpole and bubble coefficients.
use super::*;

/// OneLOop's `bnlog(rank, x)` at `x = mass/difference`, the zero-momentum
/// root of the tensor-bubble anchors, with `difference = mass - other_mass`.
///
/// OneLOop forms `x = (mass - i eps |mass|)/difference`, so a real `x` in
/// `(0, 1)`, where `log(1 - 1/x)` sits on its cut, approaches from below
/// for a positive `difference` and from above for a negative one. The latter
/// occurs for real masses of opposite sign.
pub(crate) fn bnlog(rank: usize, mass: &Atom, difference: &Atom) -> Atom {
    debug_assert!(rank <= 4);
    let root = mass / difference;
    if mass.is_zero() {
        return Atom::num(0);
    }

    let harmonic = || {
        -(1..=rank + 1)
            .map(|index| Atom::num(1) / index)
            .fold(Atom::num(0), |sum, term| sum + term)
    };
    let one_minus_root = 1 - &root;
    let mut root_power = Atom::num(1);
    let mut polynomial = Atom::num(0);
    for power in 0..=rank {
        polynomial += &root_power / (rank + 1 - power);
        root_power *= &root;
    }
    let argument: Atom = 1 - Atom::num(1) / &root;
    let real_difference = grouped(&((difference + difference.conj()) / 2));
    let below = &real_difference - Symbol::ABS.call((&real_difference,));
    let logarithm = if_nonzero_else(&below, upper_log(&argument), physical_log(&argument));
    let regular = (1 - &root_power) * logarithm - polynomial;
    if_nonzero_else(
        &root,
        if_nonzero_else(&one_minus_root, regular, harmonic()),
        Atom::num(0),
    )
}

pub(crate) fn b0_at_zero_momentum(mass_0: &Atom, mass_1: &Atom, mu_squared: &Atom) -> Atom {
    let difference = mass_0 - mass_1;
    let unequal =
        1 - (x_log_ratio(mass_0, mu_squared) - x_log_ratio(mass_1, mu_squared)) / &difference;
    let equal = if mass_0.is_zero() {
        Atom::num(0)
    } else {
        if_nonzero(mass_0, -physical_log(&(mass_0 / mu_squared)))
    };
    if_nonzero_else(&difference, unequal, equal)
}

pub(crate) fn b0_with_first_mass_zero(momentum: &Atom, mass_1: &Atom, mu_squared: &Atom) -> Atom {
    let shifted_mass = mass_1 - momentum;
    2 - (x_log_ratio(mass_1, mu_squared) - x_log_ratio(&shifted_mass, mu_squared)) / momentum
}

/// The labelled roots `x` of `p x² + (m1 - m0 - p) x + m0` and `x - 1`,
/// each from its own quadratic so that neither loses precision when a root
/// approaches zero or one. Both share the discriminant `λ(p, m0, m1)`.
fn bubble_roots(
    momentum: &Atom,
    mass_0: &Atom,
    mass_1: &Atom,
    discriminant: &Atom,
) -> [(Atom, Atom); 2] {
    let linear = mass_1 - mass_0 - momentum;
    let shifted_linear = mass_1 - mass_0 + momentum;
    [
        labelled_quadratic_roots(momentum, &linear, mass_0, discriminant),
        labelled_quadratic_roots(momentum, &shifted_linear, mass_1, discriminant),
    ]
}

pub(crate) fn b0_with_nonzero_momentum_and_mass(
    momentum: &Atom,
    mass_0: &Atom,
    mass_1: &Atom,
    mu_squared: &Atom,
) -> Atom {
    let linear = mass_1 - mass_0 - momentum;
    let discriminant = (&linear * &linear - Atom::num(4) * momentum * mass_0).sqrt();
    let [(root_plus, root_minus), (shifted_plus, shifted_minus)] =
        bubble_roots(momentum, mass_0, mass_1, &discriminant);

    // bnlog(0, x) = (1 - x) log(1 - 1/x) - 1 = -y log(y/x) - 1 with y = x - 1.
    // Forming 1 - 1/x from a root within an ulp of one cancels completely,
    // e.g. log(0) for m1 << m0; the product y log(y/x) vanishes with y.
    // With the Feynman -i0 in the quadratic, real roots approach from
    // opposite sides: root_plus from above and root_minus from below.
    let term = |shifted: &Atom, root: &Atom, upper: bool| {
        let ratio = shifted / root;
        let logarithm = if upper {
            upper_log(&ratio)
        } else {
            physical_log(&ratio)
        };
        if_nonzero(shifted, shifted * logarithm)
    };
    Atom::num(2) - physical_log(&(mass_0 / mu_squared))
        + term(&shifted_plus, &root_plus, true)
        + term(&shifted_minus, &root_minus, false)
}

pub(crate) fn b0_with_nonzero_momentum(
    momentum: &Atom,
    mass_0: &Atom,
    mass_1: &Atom,
    mu_squared: &Atom,
) -> Atom {
    if mass_0.is_zero() {
        return b0_with_first_mass_zero(momentum, mass_1, mu_squared);
    }
    // B0 is symmetric in its masses; OneLOop orders them so that a massless
    // line is always handled by the one-mass formula. The two-mass roots at
    // m1 = 0 are exactly m0/p and 1, but the rounded quadratic formula can
    // return 1 - ulp there, turning bnlog's log(1 - 1/x) into log(0).
    let second_massless = b0_with_first_mass_zero(momentum, mass_0, mu_squared);
    if mass_1.is_zero() {
        return second_massless;
    }
    let massive = b0_with_nonzero_momentum_and_mass(momentum, mass_0, mass_1, mu_squared);
    let first_massless = b0_with_first_mass_zero(momentum, mass_1, mu_squared);
    if_nonzero_else(
        mass_0,
        if_nonzero_else(mass_1, massive, second_massless),
        first_massless,
    )
}

pub(crate) fn ordinary_log3(value: &Atom) -> Atom {
    let delta = value - 1;
    let regular = ((physical_log(value) / &delta - 1) / &delta + Atom::num(1) / 2) / &delta;
    if_nonzero_else(&delta, regular, Atom::num(1) / 3)
}

pub(crate) fn db0_at_zero_momentum(mass_0: &Atom, mass_1: &Atom) -> Atom {
    if mass_1.is_zero() {
        return if mass_0.is_zero() {
            Atom::num(Coefficient::Indeterminate)
        } else {
            Atom::num(1) / (Atom::num(2) * mass_0)
        };
    }
    if mass_0.is_zero() {
        return Atom::num(1) / (Atom::num(2) * mass_1);
    }

    let difference = mass_0 - mass_1;
    let ratio = mass_0 / mass_1;
    // Original dbub0 evaluates olog3(m0/m1,0) after ordering the masses.
    // For opposite-sign real masses its logarithm is on the upper lip when
    // the denominator is positive and on the lower lip when it is negative.
    // Preserve that ratio convention explicitly and symmetrically under mass
    // exchange; both forms are the same principal function off the axis.
    let denominator_real = mass_1 + mass_1.conj();
    let negative_denominator = &denominator_real - Symbol::ABS.call((&denominator_real,));
    let continued_log3 = if_nonzero_else(
        &negative_denominator,
        ordinary_log3(&ratio),
        ordinary_log3(&ratio.conj()).conj(),
    );
    let unequal = (Atom::num(1) / 2 - &ratio * continued_log3) / mass_1;
    let equal = Atom::num(1) / (Atom::num(6) * mass_1);
    let both_nonzero = if_nonzero_else(&difference, unequal, equal);
    let first_massless = Atom::num(1) / (Atom::num(2) * mass_1);
    // Retain a symbolic denominator in the undefined fully scaleless branch:
    // runtime substitution must not silently turn the derivative into zero.
    let second_massless = if_nonzero_else(
        mass_0,
        Atom::num(1) / (Atom::num(2) * mass_0),
        Atom::num(1) / (mass_0 + mass_1),
    );

    if_nonzero_else(
        mass_1,
        if_nonzero_else(mass_0, both_nonzero, first_massless),
        second_massless,
    )
}

pub(crate) fn db0_with_one_massless(momentum: &Atom, mass: &Atom, mu_squared: &Atom) -> Atom {
    let shifted_mass = mass - momentum;
    let logarithmic_difference =
        x_log_ratio(mass, mu_squared) - x_log_ratio(&shifted_mass, mu_squared);
    &logarithmic_difference / (momentum * momentum)
        - (physical_log(&(&shifted_mass / mu_squared)) + 1) / momentum
}

pub(crate) fn db0_with_two_masses(
    momentum: &Atom,
    mass_0: &Atom,
    mass_1: &Atom,
    mu_squared: &Atom,
) -> Atom {
    let mass_difference = mass_1 - mass_0;
    let linear = &mass_difference - momentum;
    let discriminant_squared = &linear * &linear - Atom::num(4) * momentum * mass_0;
    let discriminant = discriminant_squared.sqrt();
    let [(root_plus, root_minus), (shifted_plus, shifted_minus)] =
        bubble_roots(momentum, mass_0, mass_1, &discriminant);
    // The ratios 1 - 1/x = (x - 1)/x are positive for real roots outside
    // [0,1], avoiding a subtraction of negative-axis logs with different IEEE
    // signed zeros. Inside the interval, the plus/minus roots require
    // upper/lower lips.
    //
    // They are the roots (k ± d)/(2 m0) of m0 r² - k r + m1 = 0, with
    // k = m0 + m1 - p. Both are evaluated directly unless one of them cancels
    // (|k ± d| below a quarter of |k ∓ d|, e.g. for m1 << m0); that one is
    // then the quotient of the cancellation-free roots, y/x. Near the
    // pseudo-threshold, where d is small, the direct pair keeps the rounding
    // errors of both ratios correlated in the divided difference
    // (log r+ - log r-)/d below, which a Vieta quotient would decorrelate.
    // Both forms are exact; the choice only affects rounding.
    let k = mass_0 + mass_1 - momentum;
    let ratio_numerators = [&k + &discriminant, &k - &discriminant];
    let sizes = ratio_numerators
        .each_ref()
        .map(|numerator| Symbol::ABS.call((numerator,)));
    // Nonzero exactly when 4 |own| < |other|.
    let cancels = |own: &Atom, other: &Atom| {
        let margin = grouped(&(Atom::num(4) * own - other));
        &margin - Symbol::ABS.call((&margin,))
    };
    let [plus_numerator, minus_numerator] = ratio_numerators;
    let plus_ratio = if_nonzero_else(
        &cancels(&sizes[0], &sizes[1]),
        shifted_plus / root_plus,
        plus_numerator / (Atom::num(2) * mass_0),
    );
    let minus_ratio = if_nonzero_else(
        &cancels(&sizes[1], &sizes[0]),
        shifted_minus / root_minus,
        minus_numerator / (Atom::num(2) * mass_0),
    );
    let inverse_quadratic_integral =
        (upper_log(&plus_ratio) - physical_log(&minus_ratio)) / &discriminant;
    let momentum_squared = momentum * momentum;
    let logarithmic_term = &mass_difference
        * (physical_log(&(mass_1 / mu_squared)) - physical_log(&(mass_0 / mu_squared)))
        / (Atom::num(2) * &momentum_squared);
    let inverse_quadratic_coefficient = (momentum * (mass_0 + mass_1)
        - &mass_difference * &mass_difference)
        / (Atom::num(2) * momentum_squared);

    let distinct_roots = -Atom::num(1) / momentum
        + logarithmic_term
        + inverse_quadratic_coefficient * inverse_quadratic_integral;
    // Original dbub0's coincident-root convention. The regular pseudo-
    // threshold is a genuine finite limit; at the normal threshold this
    // instead reproduces OneLOop's prescribed value, not the divergent
    // one-sided derivative limit. Root location cannot be classified solely
    // by s-m0-m1 when real mass squares can be negative.
    // The public Fortran wrapper orders the two masses by |Re(m)|+|Im(m)|
    // before applying this prescription. Unlike the regular divided
    // difference, its normal-threshold value is sensitive to the root's
    // label, so preserve that ordering when choosing the coincident root.
    let size = |mass: &Atom| {
        Symbol::ABS.call((mass + mass.conj(),))
            + Symbol::ABS.call(((mass - mass.conj()) / sheet_exact::i(),))
    };
    let size_difference = size(mass_1) - size(mass_0);
    let swap = &size_difference - Symbol::ABS.call((&size_difference,));
    let ordered_difference = if_nonzero_else(&swap, -&mass_difference, mass_difference.clone());
    let root = (momentum + ordered_difference) / (Atom::num(2) * momentum);
    let ratio = &root / (&root - 1);
    let negative_momentum = momentum - Symbol::ABS.call((momentum,));
    let logarithm = if_nonzero_else(&negative_momentum, upper_log(&ratio), physical_log(&ratio));
    let coincident_roots = ((Atom::num(2) * root - 1) * logarithm - 2) / momentum;
    if_nonzero_else(&discriminant_squared, distinct_roots, coincident_roots)
}

pub(crate) fn db0_regular_nonzero_momentum(
    momentum: &Atom,
    mass_0: &Atom,
    mass_1: &Atom,
    mu_squared: &Atom,
) -> Atom {
    if mass_0.is_zero() {
        return db0_with_one_massless(momentum, mass_1, mu_squared);
    }
    if mass_1.is_zero() {
        return db0_with_one_massless(momentum, mass_0, mu_squared);
    }

    if_nonzero_else(
        mass_0,
        if_nonzero_else(
            mass_1,
            db0_with_two_masses(momentum, mass_0, mass_1, mu_squared),
            db0_with_one_massless(momentum, mass_0, mu_squared),
        ),
        db0_with_one_massless(momentum, mass_1, mu_squared),
    )
}

pub(crate) fn db0_at_nonzero_momentum(
    momentum: &Atom,
    mass_0: &Atom,
    mass_1: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let regular = db0_regular_nonzero_momentum(momentum, mass_0, mass_1, mu_squared);
    // The original on-shell branch takes the logarithm of the reciprocal
    // scale ratio on the upper lip, including negative real shell masses.
    let shell_log = upper_log(&(mu_squared / momentum));
    let shell_finite = -(Atom::num(1) + shell_log / 2) / momentum;
    let shell_pole = -Atom::num(1) / (Atom::num(2) * momentum);

    let finite = if_nonzero_else(
        mass_0,
        if_nonzero_else(
            mass_1,
            regular.clone(),
            if_nonzero_else(&(momentum - mass_0), regular.clone(), shell_finite.clone()),
        ),
        if_nonzero_else(&(momentum - mass_1), regular, shell_finite),
    );
    let pole = if_nonzero_else(
        mass_0,
        if_nonzero_else(
            mass_1,
            Atom::num(0),
            if_nonzero_else(&(momentum - mass_0), Atom::num(0), shell_pole.clone()),
        ),
        if_nonzero_else(&(momentum - mass_1), Atom::num(0), shell_pole),
    );

    LaurentSeries::new(finite, pole, Atom::num(0))
}

/// Constructs the scalar tadpole A0 as exact Symbolica expressions.
///
/// `mass_squared` is the squared propagator mass and `mu_squared` is the
/// renormalization scale squared. A zero mass gives the scaleless result zero.
pub fn a0(mass_squared: &Atom, mu_squared: &Atom) -> LaurentSeries {
    if mass_squared.is_zero() {
        return LaurentSeries::new(Atom::num(0), Atom::num(0), Atom::num(0));
    }

    let logarithm = physical_log(&(mass_squared / mu_squared));
    let finite = if_nonzero(mass_squared, mass_squared * (1 - logarithm));

    LaurentSeries::new(finite, mass_squared.clone(), Atom::num(0))
}

/// Constructs the O(epsilon) coefficient of the scalar tadpole A0,
/// `m²(1 - L + L²/2 + π²/6)` with `L = log(m²/mu²)`.
///
/// It uses the triangle and box normalization, dividing by
/// `Γ(1+ε)Γ(1-ε)²/Γ(1-2ε)`, so a reduction may combine it with their double
/// poles. The alternative `exp(γ_E ε)` convention agrees through the
/// finite term but gives `π²m²/12` less here. A zero mass gives zero.
pub fn a0_epsilon(mass_squared: &Atom, mu_squared: &Atom) -> Atom {
    if mass_squared.is_zero() {
        return Atom::num(0);
    }

    let logarithm = physical_log(&(mass_squared / mu_squared));
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    let series = Atom::num(1) - &logarithm + &logarithm * &logarithm / 2 + pi_squared / 6;
    if_nonzero(mass_squared, mass_squared * series)
}

/// Constructs all tadpole coefficients through tensor rank four.
///
/// Coefficients are returned in the order A0, A00, A0000, truncated to the
/// requested rank. Odd ranks have the same coefficient set as the preceding
/// even rank, matching the OneLOop interface.
pub fn an(
    rank: usize,
    mass_squared: &Atom,
    mu_squared: &Atom,
) -> Result<Vec<LaurentSeries>, UnsupportedRank> {
    if rank > 4 {
        return Err(UnsupportedRank(rank));
    }

    if mass_squared.is_zero() {
        return Ok((0..rank / 2 + 1)
            .map(|_| LaurentSeries::new(Atom::num(0), Atom::num(0), Atom::num(0)))
            .collect());
    }

    let logarithm = physical_log(&(mass_squared / mu_squared));
    let mut coefficients = Vec::with_capacity(rank / 2 + 1);
    coefficients.push(a0(mass_squared, mu_squared));

    if rank >= 2 {
        let finite = mass_squared * mass_squared * (Atom::num(3) / 2 - &logarithm) / 4;
        let pole = mass_squared * mass_squared / 4;
        coefficients.push(LaurentSeries::new(
            if_nonzero(mass_squared, finite),
            pole,
            Atom::num(0),
        ));
    }

    if rank >= 4 {
        let mass_cubed = mass_squared * mass_squared * mass_squared;
        let finite = &mass_cubed * (Atom::num(11) / 6 - logarithm) / 24;
        coefficients.push(LaurentSeries::new(
            if_nonzero(mass_squared, finite),
            mass_cubed / 24,
            Atom::num(0),
        ));
    }

    Ok(coefficients)
}

/// Constructs the scalar bubble B0 as exact Symbolica expressions.
///
/// `momentum_squared` is real, the squared masses must have non-positive
/// imaginary parts, and `mu_squared` is the renormalization scale squared.
/// Numerical proximity branches from the standalone evaluator are deliberately
/// omitted. Evaluation uses the selected Symbolica numeric domain; fixed-f64
/// evaluation can suffer severe cancellation. This crate supplies no automatic
/// precision escalation.
pub fn b0(
    momentum_squared: &Atom,
    mass_0_squared: &Atom,
    mass_1_squared: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    if momentum_squared.is_zero() && mass_0_squared.is_zero() && mass_1_squared.is_zero() {
        return LaurentSeries::new(Atom::num(0), Atom::num(0), Atom::num(0));
    }

    let zero_momentum = b0_at_zero_momentum(mass_0_squared, mass_1_squared, mu_squared);
    let nonzero_momentum =
        b0_with_nonzero_momentum(momentum_squared, mass_0_squared, mass_1_squared, mu_squared);
    let finite = if_nonzero_else(momentum_squared, nonzero_momentum, zero_momentum);
    let simple_pole = if_nonzero_else(
        momentum_squared,
        Atom::num(1),
        if_nonzero_else(
            mass_0_squared,
            Atom::num(1),
            if_nonzero(mass_1_squared, Atom::num(1)),
        ),
    );

    LaurentSeries::new(finite, simple_pole, Atom::num(0))
}

/// Constructs dB0/dp² as exact Symbolica expressions.
///
/// The regular nonzero-momentum branch uses an analytic derivative of B0.
/// Explicit zero-momentum and massless on-shell limits preserve the full
/// Laurent structure of OneLOop.
/// At exact normal thresholds, the coincident-root branch reproduces OneLOop's
/// prescribed value, not the divergent one-sided derivative limit.
///
/// The fully scaleless point `p² = m0² = m1² = 0` is undefined and unsupported.
/// Direct construction marks it indeterminate; substituting that point into a
/// generic expression produces a nonfinite value.
pub fn db0(
    momentum_squared: &Atom,
    mass_0_squared: &Atom,
    mass_1_squared: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    if momentum_squared.is_zero() && mass_0_squared.is_zero() && mass_1_squared.is_zero() {
        let indeterminate = Atom::num(Coefficient::Indeterminate);
        return LaurentSeries::new(indeterminate.clone(), indeterminate, Atom::num(0));
    }

    if momentum_squared.is_zero() {
        return LaurentSeries::new(
            db0_at_zero_momentum(mass_0_squared, mass_1_squared),
            Atom::num(0),
            Atom::num(0),
        );
    }

    let nonzero =
        db0_at_nonzero_momentum(momentum_squared, mass_0_squared, mass_1_squared, mu_squared);
    let [nonzero_finite, nonzero_pole, _] = nonzero.into_coefficients();
    let finite = if_nonzero_else(
        momentum_squared,
        nonzero_finite,
        db0_at_zero_momentum(mass_0_squared, mass_1_squared),
    );
    let pole = if_nonzero_else(momentum_squared, nonzero_pole, Atom::num(0));

    LaurentSeries::new(finite, pole, Atom::num(0))
}

pub(crate) fn b1_zero_momentum_anchor(mass_0: &Atom, mass_1: &Atom, mu_squared: &Atom) -> Atom {
    let logarithm = physical_log(&(mass_1 / mu_squared));
    let difference = mass_1 - mass_0;
    let unequal = (&logarithm + bnlog(1, mass_1, &difference)) / 2;
    if_nonzero_else(&difference, unequal, logarithm / 2)
}

pub(crate) fn b1_at_zero_momentum(mass_0: &Atom, mass_1: &Atom, mu_squared: &Atom) -> Atom {
    if mass_0.is_zero() && mass_1.is_zero() {
        return Atom::num(0);
    }
    if mass_1.is_zero() {
        return b1_zero_momentum_anchor(mass_1, mass_0, mu_squared);
    }

    let anchored_at_mass_1 = -b0_at_zero_momentum(mass_0, mass_1, mu_squared)
        - b1_zero_momentum_anchor(mass_0, mass_1, mu_squared);
    if mass_0.is_zero() {
        return anchored_at_mass_1;
    }

    let anchored_at_mass_0 = b1_zero_momentum_anchor(mass_1, mass_0, mu_squared);
    if_nonzero_else(
        mass_1,
        anchored_at_mass_1,
        if_nonzero(mass_0, anchored_at_mass_0),
    )
}

pub(crate) fn b1_at_nonzero_momentum(
    momentum: &Atom,
    mass_0: &Atom,
    mass_1: &Atom,
    mu_squared: &Atom,
) -> Atom {
    let a0_mass_0 = a0(mass_0, mu_squared).into_coefficients()[0].clone();
    let a0_mass_1 = a0(mass_1, mu_squared).into_coefficients()[0].clone();
    let b0_finite = b0(momentum, mass_0, mass_1, mu_squared).into_coefficients()[0].clone();

    (a0_mass_0 - a0_mass_1 + (mass_1 - mass_0 - momentum) * b0_finite) / (Atom::num(2) * momentum)
}

/// Constructs the rank-one longitudinal bubble coefficient B1.
///
/// The generic branch is the exact Passarino--Veltman reduction to A0 and B0.
/// Its removable zero-momentum limit is expressed directly with native logs
/// and conditionals, so the returned expression remains fully transparent.
pub fn b1(
    momentum_squared: &Atom,
    mass_0_squared: &Atom,
    mass_1_squared: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    if momentum_squared.is_zero() && mass_0_squared.is_zero() && mass_1_squared.is_zero() {
        return LaurentSeries::new(Atom::num(0), Atom::num(0), Atom::num(0));
    }

    let zero_momentum = b1_at_zero_momentum(mass_0_squared, mass_1_squared, mu_squared);
    let nonzero_momentum =
        b1_at_nonzero_momentum(momentum_squared, mass_0_squared, mass_1_squared, mu_squared);
    let finite = if_nonzero_else(momentum_squared, nonzero_momentum, zero_momentum);
    let simple_pole = if_nonzero_else(
        momentum_squared,
        -Atom::num(1) / 2,
        if_nonzero_else(
            mass_0_squared,
            -Atom::num(1) / 2,
            if_nonzero(mass_1_squared, -Atom::num(1) / 2),
        ),
    );

    LaurentSeries::new(finite, simple_pole, Atom::num(0))
}

/// Constructs the transverse rank-two bubble coefficient B00.
///
/// The factor `1/(d-1)` is expanded through finite order for `d=4-2 epsilon`,
/// including the rational finite term generated by the simple pole.
pub fn b00(
    momentum_squared: &Atom,
    mass_0_squared: &Atom,
    mass_1_squared: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let a0_mass_1 = a0(mass_1_squared, mu_squared).into_coefficients();
    let b0_coefficients =
        b0(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let b1_coefficients =
        b1(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let kinematic = mass_1_squared - mass_0_squared - momentum_squared;
    let numerator_finite = &a0_mass_1[0] + Atom::num(2) * mass_0_squared * &b0_coefficients[0]
        - &kinematic * &b1_coefficients[0];
    let numerator_pole = &a0_mass_1[1] + Atom::num(2) * mass_0_squared * &b0_coefficients[1]
        - kinematic * &b1_coefficients[1];

    LaurentSeries::new(
        numerator_finite / 6 + &numerator_pole / 9,
        numerator_pole / 6,
        Atom::num(0),
    )
}

pub(crate) fn b11_zero_momentum_anchor(mass_0: &Atom, mass_1: &Atom, mu_squared: &Atom) -> Atom {
    let logarithm = physical_log(&(mass_1 / mu_squared));
    let difference = mass_1 - mass_0;
    let unequal = (-&logarithm - bnlog(2, mass_1, &difference)) / 3;
    if_nonzero_else(&difference, unequal, -logarithm / 3)
}

pub(crate) fn b11_at_zero_momentum(mass_0: &Atom, mass_1: &Atom, mu_squared: &Atom) -> Atom {
    if mass_0.is_zero() && mass_1.is_zero() {
        return Atom::num(0);
    }
    if mass_1.is_zero() {
        return b11_zero_momentum_anchor(mass_1, mass_0, mu_squared);
    }

    let anchored_at_mass_1 = b11_zero_momentum_anchor(mass_0, mass_1, mu_squared)
        + Atom::num(2) * b1_zero_momentum_anchor(mass_0, mass_1, mu_squared)
        + b0_at_zero_momentum(mass_0, mass_1, mu_squared);
    if mass_0.is_zero() {
        return anchored_at_mass_1;
    }

    let anchored_at_mass_0 = b11_zero_momentum_anchor(mass_1, mass_0, mu_squared);
    if_nonzero_else(
        mass_1,
        anchored_at_mass_1,
        if_nonzero(mass_0, anchored_at_mass_0),
    )
}

/// Constructs the longitudinal rank-two bubble coefficient B11.
///
/// Away from zero momentum this is the exact reduction to A0, B1, and B00.
/// The removable zero-momentum limit is a native logarithmic expression.
pub fn b11(
    momentum_squared: &Atom,
    mass_0_squared: &Atom,
    mass_1_squared: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    if momentum_squared.is_zero() && mass_0_squared.is_zero() && mass_1_squared.is_zero() {
        return LaurentSeries::new(Atom::num(0), Atom::num(0), Atom::num(0));
    }

    let zero_finite = b11_at_zero_momentum(mass_0_squared, mass_1_squared, mu_squared);
    let a0_mass_1 = a0(mass_1_squared, mu_squared).into_coefficients();
    let b1_coefficients =
        b1(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let b00_coefficients =
        b00(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let kinematic = mass_1_squared - mass_0_squared - momentum_squared;
    let denominator = Atom::num(2) * momentum_squared;
    let nonzero_finite = (&a0_mass_1[0] + &kinematic * &b1_coefficients[0]
        - Atom::num(2) * &b00_coefficients[0])
        / &denominator;
    let nonzero_pole = (&a0_mass_1[1] + kinematic * &b1_coefficients[1]
        - Atom::num(2) * &b00_coefficients[1])
        / denominator;
    let finite = if_nonzero_else(momentum_squared, nonzero_finite, zero_finite);
    let pole = if_nonzero_else(
        momentum_squared,
        nonzero_pole,
        if_nonzero_else(
            mass_0_squared,
            Atom::num(1) / 3,
            if_nonzero(mass_1_squared, Atom::num(1) / 3),
        ),
    );

    LaurentSeries::new(finite, pole, Atom::num(0))
}

/// Constructs the mixed rank-three bubble coefficient B001.
///
/// The factor `1/d` is expanded through finite order for `d=4-2 epsilon`,
/// so the pole contributes its required finite rational term.
pub fn b001(
    momentum_squared: &Atom,
    mass_0_squared: &Atom,
    mass_1_squared: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let a0_mass_1 = a0(mass_1_squared, mu_squared).into_coefficients();
    let b1_coefficients =
        b1(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let b11_coefficients =
        b11(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let kinematic = mass_1_squared - mass_0_squared - momentum_squared;
    let numerator_finite = Atom::num(2) * mass_0_squared * &b1_coefficients[0]
        - &a0_mass_1[0]
        - &kinematic * &b11_coefficients[0];
    let numerator_pole = Atom::num(2) * mass_0_squared * &b1_coefficients[1]
        - &a0_mass_1[1]
        - kinematic * &b11_coefficients[1];

    LaurentSeries::new(
        numerator_finite / 8 + &numerator_pole / 16,
        numerator_pole / 8,
        Atom::num(0),
    )
}

pub(crate) fn b111_zero_momentum_anchor(mass_0: &Atom, mass_1: &Atom, mu_squared: &Atom) -> Atom {
    let logarithm = physical_log(&(mass_1 / mu_squared));
    let difference = mass_1 - mass_0;
    let unequal = (&logarithm + bnlog(3, mass_1, &difference)) / 4;
    if_nonzero_else(&difference, unequal, logarithm / 4)
}

pub(crate) fn b111_at_zero_momentum(mass_0: &Atom, mass_1: &Atom, mu_squared: &Atom) -> Atom {
    if mass_0.is_zero() && mass_1.is_zero() {
        return Atom::num(0);
    }
    if mass_1.is_zero() {
        return b111_zero_momentum_anchor(mass_1, mass_0, mu_squared);
    }

    let anchored_at_mass_1 = -b111_zero_momentum_anchor(mass_0, mass_1, mu_squared)
        - Atom::num(3) * b11_zero_momentum_anchor(mass_0, mass_1, mu_squared)
        - Atom::num(3) * b1_zero_momentum_anchor(mass_0, mass_1, mu_squared)
        - b0_at_zero_momentum(mass_0, mass_1, mu_squared);
    if mass_0.is_zero() {
        return anchored_at_mass_1;
    }

    let anchored_at_mass_0 = b111_zero_momentum_anchor(mass_1, mass_0, mu_squared);
    if_nonzero_else(
        mass_1,
        anchored_at_mass_1,
        if_nonzero(mass_0, anchored_at_mass_0),
    )
}

/// Constructs the longitudinal rank-three bubble coefficient B111.
pub fn b111(
    momentum_squared: &Atom,
    mass_0_squared: &Atom,
    mass_1_squared: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    if momentum_squared.is_zero() && mass_0_squared.is_zero() && mass_1_squared.is_zero() {
        return LaurentSeries::new(Atom::num(0), Atom::num(0), Atom::num(0));
    }

    let zero_finite = b111_at_zero_momentum(mass_0_squared, mass_1_squared, mu_squared);
    let a0_mass_1 = a0(mass_1_squared, mu_squared).into_coefficients();
    let b11_coefficients =
        b11(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let b001_coefficients =
        b001(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let kinematic = mass_1_squared - mass_0_squared - momentum_squared;
    let denominator = Atom::num(2) * momentum_squared;
    let nonzero_finite = (-&a0_mass_1[0] + &kinematic * &b11_coefficients[0]
        - Atom::num(4) * &b001_coefficients[0])
        / &denominator;
    let nonzero_pole = (-&a0_mass_1[1] + kinematic * &b11_coefficients[1]
        - Atom::num(4) * &b001_coefficients[1])
        / denominator;
    let finite = if_nonzero_else(momentum_squared, nonzero_finite, zero_finite);
    let pole = if_nonzero_else(
        momentum_squared,
        nonzero_pole,
        if_nonzero_else(
            mass_0_squared,
            -Atom::num(1) / 4,
            if_nonzero(mass_1_squared, -Atom::num(1) / 4),
        ),
    );

    LaurentSeries::new(finite, pole, Atom::num(0))
}

pub(crate) fn a00(mass_squared: &Atom, mu_squared: &Atom) -> LaurentSeries {
    an(2, mass_squared, mu_squared)
        .expect("rank two is supported")
        .pop()
        .expect("A00 is present at rank two")
}

/// Constructs the double-metric rank-four bubble coefficient B0000.
pub fn b0000(
    momentum_squared: &Atom,
    mass_0_squared: &Atom,
    mass_1_squared: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let a00_mass_1 = a00(mass_1_squared, mu_squared).into_coefficients();
    let b00_coefficients =
        b00(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let b001_coefficients =
        b001(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let kinematic = mass_1_squared - mass_0_squared - momentum_squared;
    let numerator_finite = Atom::num(2) * mass_0_squared * &b00_coefficients[0] + &a00_mass_1[0]
        - &kinematic * &b001_coefficients[0];
    let numerator_pole = Atom::num(2) * mass_0_squared * &b00_coefficients[1] + &a00_mass_1[1]
        - kinematic * &b001_coefficients[1];

    LaurentSeries::new(
        numerator_finite / 10 + &numerator_pole / 25,
        numerator_pole / 10,
        Atom::num(0),
    )
}

/// Constructs the mixed rank-four bubble coefficient B0011.
pub fn b0011(
    momentum_squared: &Atom,
    mass_0_squared: &Atom,
    mass_1_squared: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let a0_mass_1 = a0(mass_1_squared, mu_squared).into_coefficients();
    let b11_coefficients =
        b11(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let b111_coefficients =
        b111(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let kinematic = mass_1_squared - mass_0_squared - momentum_squared;
    let numerator_finite = Atom::num(2) * mass_0_squared * &b11_coefficients[0] + &a0_mass_1[0]
        - &kinematic * &b111_coefficients[0];
    let numerator_pole = Atom::num(2) * mass_0_squared * &b11_coefficients[1] + &a0_mass_1[1]
        - kinematic * &b111_coefficients[1];

    LaurentSeries::new(
        numerator_finite / 10 + &numerator_pole / 25,
        numerator_pole / 10,
        Atom::num(0),
    )
}

pub(crate) fn b1111_zero_momentum_anchor(mass_0: &Atom, mass_1: &Atom, mu_squared: &Atom) -> Atom {
    let logarithm = physical_log(&(mass_1 / mu_squared));
    let difference = mass_1 - mass_0;
    let unequal = (-&logarithm - bnlog(4, mass_1, &difference)) / 5;
    if_nonzero_else(&difference, unequal, -logarithm / 5)
}

pub(crate) fn b1111_at_zero_momentum(mass_0: &Atom, mass_1: &Atom, mu_squared: &Atom) -> Atom {
    if mass_0.is_zero() && mass_1.is_zero() {
        return Atom::num(0);
    }
    if mass_1.is_zero() {
        return b1111_zero_momentum_anchor(mass_1, mass_0, mu_squared);
    }

    // Reflection x -> 1-x of the Feynman parameter: with B_n = (-1)^n ∫x^n f,
    // ∫(1-x)^4 f = B0 + 4 B1 + 6 B11 + 4 B111 + B1111. OneLOop's bub1111
    // negates this sum in its mass-swap reconstruction (taken when m1 is the
    // larger mass), so its zero-momentum B1111 has the wrong sign there and
    // disagrees with the p² -> 0 limit and with ∫x^4 [1/eps - log D].
    let anchored_at_mass_1 = b1111_zero_momentum_anchor(mass_0, mass_1, mu_squared)
        + Atom::num(4) * b111_zero_momentum_anchor(mass_0, mass_1, mu_squared)
        + Atom::num(6) * b11_zero_momentum_anchor(mass_0, mass_1, mu_squared)
        + Atom::num(4) * b1_zero_momentum_anchor(mass_0, mass_1, mu_squared)
        + b0_at_zero_momentum(mass_0, mass_1, mu_squared);
    if mass_0.is_zero() {
        return anchored_at_mass_1;
    }

    let anchored_at_mass_0 = b1111_zero_momentum_anchor(mass_1, mass_0, mu_squared);
    if_nonzero_else(
        mass_1,
        anchored_at_mass_1,
        if_nonzero(mass_0, anchored_at_mass_0),
    )
}

/// Constructs the longitudinal rank-four bubble coefficient B1111.
pub fn b1111(
    momentum_squared: &Atom,
    mass_0_squared: &Atom,
    mass_1_squared: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    if momentum_squared.is_zero() && mass_0_squared.is_zero() && mass_1_squared.is_zero() {
        return LaurentSeries::new(Atom::num(0), Atom::num(0), Atom::num(0));
    }

    let zero_finite = b1111_at_zero_momentum(mass_0_squared, mass_1_squared, mu_squared);
    let a0_mass_1 = a0(mass_1_squared, mu_squared).into_coefficients();
    let b111_coefficients =
        b111(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let b0011_coefficients =
        b0011(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).into_coefficients();
    let kinematic = mass_1_squared - mass_0_squared - momentum_squared;
    let denominator = Atom::num(2) * momentum_squared;
    let nonzero_finite = (&a0_mass_1[0] + &kinematic * &b111_coefficients[0]
        - Atom::num(6) * &b0011_coefficients[0])
        / &denominator;
    let nonzero_pole = (&a0_mass_1[1] + kinematic * &b111_coefficients[1]
        - Atom::num(6) * &b0011_coefficients[1])
        / denominator;
    let finite = if_nonzero_else(momentum_squared, nonzero_finite, zero_finite);
    let pole = if_nonzero_else(
        momentum_squared,
        nonzero_pole,
        if_nonzero_else(
            mass_0_squared,
            Atom::num(1) / 5,
            if_nonzero(mass_1_squared, Atom::num(1) / 5),
        ),
    );

    LaurentSeries::new(finite, pole, Atom::num(0))
}

/// Constructs all bubble coefficients through tensor rank four.
///
/// The order is `B0, B1, B00, B11, B001, B111, B0000, B0011, B1111`,
/// truncated at the requested rank.
pub fn bn(
    rank: usize,
    momentum_squared: &Atom,
    mass_0_squared: &Atom,
    mass_1_squared: &Atom,
    mu_squared: &Atom,
) -> Result<Vec<LaurentSeries>, UnsupportedRank> {
    if rank > 4 {
        return Err(UnsupportedRank(rank));
    }
    let mut result = Vec::with_capacity([1, 2, 4, 6, 9][rank]);
    result.push(b0(
        momentum_squared,
        mass_0_squared,
        mass_1_squared,
        mu_squared,
    ));
    if rank >= 1 {
        result.push(b1(
            momentum_squared,
            mass_0_squared,
            mass_1_squared,
            mu_squared,
        ));
    }
    if rank >= 2 {
        result.push(b00(
            momentum_squared,
            mass_0_squared,
            mass_1_squared,
            mu_squared,
        ));
        result.push(b11(
            momentum_squared,
            mass_0_squared,
            mass_1_squared,
            mu_squared,
        ));
    }
    if rank >= 3 {
        result.push(b001(
            momentum_squared,
            mass_0_squared,
            mass_1_squared,
            mu_squared,
        ));
        result.push(b111(
            momentum_squared,
            mass_0_squared,
            mass_1_squared,
            mu_squared,
        ));
    }
    if rank >= 4 {
        result.push(b0000(
            momentum_squared,
            mass_0_squared,
            mass_1_squared,
            mu_squared,
        ));
        result.push(b0011(
            momentum_squared,
            mass_0_squared,
            mass_1_squared,
            mu_squared,
        ));
        result.push(b1111(
            momentum_squared,
            mass_0_squared,
            mass_1_squared,
            mu_squared,
        ));
    }
    Ok(result)
}
