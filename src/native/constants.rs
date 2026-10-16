//! Exact constants shared by the master runtime and reusable Rust emission.

use super::primitives::{C, NativeFloat, polylog2};
use symbolica::domains::{integer::Integer, rational::Rational};

/// Exact initialization data emitted alongside generated arithmetic.
/// No strings or rational conversions are used in a warmed scalar call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[cfg_attr(not(feature = "rust-codegen"), allow(dead_code))]
pub enum ConstantSpec<'a> {
    Rational { re: &'a str, im: &'a str },
    Pi,
    Polylog2 { re: &'a str, im: &'a str },
}

/// A reusable constant table for a set of generated Rust modules.
#[derive(Clone)]
pub struct ConstantContext<T: NativeFloat> {
    constants: Vec<C<T>>,
}

impl<T: NativeFloat> ConstantContext<T> {
    /// Resolve exact constants directly at the requested binary precision.
    /// Fixed-precision types require their native bit count. This constructor
    /// does not inspect, copy, or change the precision of evaluation inputs.
    pub fn new(specifications: &[ConstantSpec<'_>], bits: u32) -> Result<Self, String> {
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
        let prototype = T::zero_at(bits);
        let rational = |text: &str| -> Result<T, String> {
            let (num, den) = text.split_once('/').unwrap_or((text, "1"));
            let num = num.parse::<Integer>().map_err(|error| error.to_string())?;
            let den = den.parse::<Integer>().map_err(|error| error.to_string())?;
            if den == 0 {
                return Err("invalid generated rational denominator".into());
            }
            Ok(prototype.from_rational(&Rational::from((num, den))))
        };
        let constants = specifications
            .iter()
            .map(|spec| {
                Ok(match spec {
                    ConstantSpec::Rational { re, im } => C::new(rational(re)?, rational(im)?),
                    ConstantSpec::Pi => C::new(prototype.pi(), prototype.zero()),
                    ConstantSpec::Polylog2 { re, im } => {
                        polylog2(&C::new(rational(re)?, rational(im)?))
                    }
                })
            })
            .collect::<Result<Vec<_>, String>>()?;
        Ok(Self { constants })
    }

    /// Read the constant at an index assigned by the emitter.
    /// Panics if the index is outside this table.
    #[inline(always)]
    pub fn constant(&self, index: usize) -> &C<T> {
        &self.constants[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::domains::float::{DoubleFloat, Float, FloatLike, SingleFloat};

    #[test]
    fn generated_constants_preserve_working_precision() {
        let specs = [ConstantSpec::Rational {
            re: "100000000000000000001/100000000000000000000",
            im: "-2/3",
        }];
        let double =
            ConstantContext::<DoubleFloat>::new(&specs, DoubleFloat::DEFAULT_BITS).unwrap();
        let wide = ConstantContext::<Float>::new(&specs, 256).unwrap();
        let machine = ConstantContext::<f64>::new(&specs, 53).unwrap();
        assert_eq!(machine.constant(0).re, 1.0);
        assert_ne!(double.constant(0).re, DoubleFloat::from(1.0));
        assert_ne!(wide.constant(0).re, Float::with_val(256, 1));
        assert_eq!(wide.constant(0).re.get_precision(), 256);
        let exact = Rational::from((
            "100000000000000000001".parse::<Integer>().unwrap(),
            "100000000000000000000".parse::<Integer>().unwrap(),
        ));
        assert_eq!(wide.constant(0).re, Float::new(256).from_rational(&exact));
    }

    #[test]
    fn invalid_constant_tables_and_precisions_are_rejected() {
        for text in ["1/0", "not-a-rational", "2/3/4"] {
            assert!(
                ConstantContext::<f64>::new(&[ConstantSpec::Rational { re: text, im: "0" }], 53)
                    .is_err()
            );
        }
        assert!(ConstantContext::<Float>::new(&[], 0).is_err());
        assert!(ConstantContext::<f64>::new(&[], 256).is_err());
        assert!(ConstantContext::<DoubleFloat>::new(&[], 53).is_err());
        assert!(ConstantContext::<Float>::new(&[], 256).is_ok());
    }
}
