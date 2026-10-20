//! A sheet as (physical argument, continued logarithm), with no parity tree.
pub(crate) use crate::sheet_exact::{i, im, negative, re, sign_nonnegative, sqrt_lower};
use crate::{Atom, AtomCore, FunctionMap, Symbol, if_nonzero_else as choose};
use std::ops::{Div, Mul};

#[derive(Clone)]
pub(crate) struct SheetAtom {
    argument: Atom,
    logarithm: Atom,
}

impl SheetAtom {
    pub(crate) fn lower(value: Atom) -> Self {
        Self::with_sign(value, Atom::num(-1))
    }
    pub(crate) fn upper(value: Atom) -> Self {
        Self::with_sign(value, Atom::num(1))
    }
    pub(crate) fn with_sign(value: Atom, zero_sign: Atom) -> Self {
        let logarithm =
            symbolica::symbol!("__olo_inspect_sheet_initial_log").call((&value, zero_sign));
        Self {
            argument: value,
            logarithm,
        }
    }
    pub(crate) fn scale_positive(&self, scale: &Atom) -> Self {
        Self {
            argument: &self.argument * scale,
            logarithm: &self.logarithm + scale.log(),
        }
    }
    pub(crate) fn value(&self) -> Atom {
        self.argument.clone()
    }
    pub(crate) fn log(&self) -> Atom {
        self.logarithm.clone()
    }
    pub(crate) fn log_over_one_minus(&self) -> Atom {
        let difference = 1 - &self.argument;
        choose(&difference, &self.logarithm / &difference, Atom::num(-1))
    }
    /// The offset `x - 1` only restores digits of the numerical maps; the
    /// inspection form is the continued logarithm itself.
    pub(crate) fn log_given(&self, _: &Atom) -> Atom {
        self.log()
    }
    pub(crate) fn dilog(&self) -> Atom {
        symbolica::symbol!("__olo_inspect_sheet_dilog").call((&self.argument, &self.logarithm))
    }
    /// The complement only restores digits of the numerical maps.
    pub(crate) fn dilog_complement(&self, _: &Atom) -> Atom {
        self.dilog()
    }
    /// The plain form of a value that the numerical maps choose between
    /// algebraically identical variants by their rounding (`numerical`).
    pub(crate) fn rounding_choice(_numerical: Atom, plain: Atom) -> Atom {
        plain
    }
    pub(crate) fn wound(&self, turns: &Atom) -> Self {
        Self {
            argument: self.argument.clone(),
            logarithm: &self.logarithm + Atom::num(2) * i() * Symbol::PI.to_atom() * turns,
        }
    }
    pub(crate) fn split_winding(&self) -> (SheetAtom, Atom) {
        let principal = self.argument.log();
        let winding = (&self.logarithm - &principal) / (i() * Symbol::PI.to_atom());
        (
            SheetAtom {
                argument: self.argument.clone(),
                logarithm: principal,
            },
            winding,
        )
    }
}

impl Mul for &SheetAtom {
    type Output = SheetAtom;
    fn mul(self, rhs: Self) -> SheetAtom {
        SheetAtom {
            argument: &self.argument * &rhs.argument,
            logarithm: &self.logarithm + &rhs.logarithm,
        }
    }
}
impl Div for &SheetAtom {
    type Output = SheetAtom;
    fn div(self, rhs: Self) -> SheetAtom {
        SheetAtom {
            argument: choose(
                &(&self.argument - &rhs.argument),
                &self.argument / &rhs.argument,
                Atom::num(1),
            ),
            logarithm: &self.logarithm - &rhs.logarithm,
        }
    }
}

pub(crate) fn divided_difference(first: &SheetAtom, second: &SheetAtom) -> Atom {
    let difference = &first.argument - &second.argument;
    choose(
        &difference,
        (first.dilog() - second.dilog()) / &difference,
        first.log_over_one_minus(),
    )
}

/// The complements `1 - x` only restore digits of the numerical maps; the
/// inspection form is the divided difference itself.
pub(crate) fn complement_difference(
    first: &SheetAtom,
    _: &Atom,
    second: &SheetAtom,
    _: &Atom,
) -> Atom {
    divided_difference(first, second)
}

pub(crate) fn register(map: &mut FunctionMap) {
    crate::sheet_exact::register_for_inspection(map);
    let parameters = [
        symbolica::symbol!("__olo_inspect_z"),
        symbolica::symbol!("__olo_inspect_log"),
    ];
    let [z, logarithm] = parameters.map(Symbol::to_atom);
    let on_axis =
        Symbol::ABS.call(&z).log() + i() * Symbol::PI.to_atom() * negative(&re(&z)) * &logarithm;
    map.add_function(
        symbolica::symbol!("__olo_inspect_sheet_initial_log"),
        parameters.to_vec(),
        choose(&im(&z), z.log(), on_axis),
    )
    .unwrap();
    // This is Euler's identity on the same two half-planes as dilog_body.
    // Outside the unit circle use the identical inversion identity, retaining
    // the full continued logarithm (including every winding).
    let base = |z: &Atom, l: &Atom| {
        let log_one_minus = (Atom::num(1) - z).log();
        let below = Symbol::PI.to_atom().pow(2) / 6 - crate::dilog_atom(z) - &log_one_minus * l;
        let above = crate::dilog_atom(1 - z) - log_one_minus * (l - z.log());
        choose(&negative(&(2 * re(z) - 1)), below, above)
    };
    let ordinary = choose(
        &negative(&(1 - &z * z.conj())),
        -base(&(Atom::num(1) / &z), &(-&logarithm)) - logarithm.pow(2) / 2,
        base(&z, &logarithm),
    );
    let body = choose(
        &z,
        choose(&(&z - 1), ordinary, Atom::num(0)),
        Symbol::PI.to_atom().pow(2) / 6,
    );
    map.add_function(
        symbolica::symbol!("__olo_inspect_sheet_dilog"),
        parameters.to_vec(),
        body,
    )
    .unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sheet_exact::SheetAtom as Original;
    use symbolica::prelude::Complex;

    #[test]
    fn continued_log_representation_preserves_sheet_operations() {
        crate::tests::large_stack(|| {
            let x = symbolica::symbol!("inspect_sheet_test_x").to_atom();
            let y = symbolica::symbol!("inspect_sheet_test_y").to_atom();
            let s = symbolica::symbol!("inspect_sheet_test_s").to_atom();
            let old_x = Original::with_sign(x.clone(), s.clone());
            let old_y = Original::with_sign(y.clone(), -&s);
            let new_x = SheetAtom::with_sign(x.clone(), s.clone());
            let new_y = SheetAtom::with_sign(y.clone(), -&s);
            let old = [
                old_x.clone(),
                &old_x * &old_y,
                &old_x / &old_y,
                &(&old_x * &old_y) * &old_x,
                &(&old_x / &old_y) / &old_y,
            ];
            let new = [
                new_x.clone(),
                &new_x * &new_y,
                &new_x / &new_y,
                &(&new_x * &new_y) * &new_x,
                &(&new_x / &new_y) / &new_y,
            ];
            let mut map = FunctionMap::new();
            register(&mut map);
            let mut evaluators = old
                .iter()
                .zip(&new)
                .map(|(o, n)| {
                    Atom::evaluator_multiple(
                        &[o.value(), n.value(), o.log(), n.log(), o.dilog(), n.dilog()],
                        &[x.clone(), y.clone(), s.clone()],
                    )
                    .function_map(map.as_symbolica().clone())
                    .direct_translation(true)
                    .horner_iterations(0)
                    .build()
                    .unwrap()
                    .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()))
                })
                .collect::<Vec<_>>();
            let points = [
                Complex::new(0.2, 0.3),
                Complex::new(2., 3.),
                Complex::new(-0.2, 0.3),
                Complex::new(-2., -3.),
                Complex::new(0.2, -0.3),
                Complex::new(2., -3.),
                Complex::new(-2., 3.),
                Complex::new(-0.2, -0.3),
                Complex::new(0.2, 0.),
                Complex::new(2., 0.),
                Complex::new(-0.2, 0.),
                Complex::new(-2., 0.),
                Complex::new(0., 1.),
                Complex::new(0., -1.),
            ];
            for x in points {
                for y in points {
                    for s in [-1., 1.] {
                        for (operation, evaluator) in evaluators.iter_mut().enumerate() {
                            let mut out = [Complex::new(0., 0.); 6];
                            evaluator.evaluate(&[x, y, Complex::new(s, 0.)], &mut out);
                            for (component, pair) in out.chunks_exact(2).enumerate() {
                                let error =
                                    (pair[0].re - pair[1].re).hypot(pair[0].im - pair[1].im);
                                assert!(
                                    error < 1e-10,
                                    "operation={operation}, component={component}, x={x}, y={y}, s={s}, old={}, new={}",
                                    pair[0],
                                    pair[1]
                                );
                            }
                        }
                    }
                }
            }
        });
    }
}
