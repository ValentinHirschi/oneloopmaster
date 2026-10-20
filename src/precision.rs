//! Reusable native scalar evaluation at a fixed arbitrary working precision.
use crate::{EvaluationBackend, NativeEvaluator, ScalarIntegral, masters};
use symbolica::{
    domains::float::{Complex, Float, SingleFloat},
    evaluate::ExpressionEvaluator,
};

type C = Complex<Float>;

#[derive(Clone)]
enum PreparedPrecision {
    Expression(ExpressionEvaluator<C>),
    Native(NativeEvaluator<Float>),
}

fn default_backend() -> EvaluationBackend {
    if crate::DEFAULT_BACKEND == EvaluationBackend::Native {
        EvaluationBackend::Native
    } else {
        EvaluationBackend::Expression
    }
}

/// A reusable three-coefficient native scalar evaluator.
///
/// Decimal requests use at least `ceil(digits * log2(10)) + 32` working bits.
/// This fixed guard is not adaptive and does not certify that all requested
/// digits survive cancellation or proximity to a singularity. Inputs are
/// rounded directly to the working precision, never through machine floats.
/// Increasing the precision of an already rounded input cannot recover its
/// missing digits. Outputs retain Numerica's computed precision; it is not
/// padded upward to the requested or working precision.
#[derive(Clone)]
pub struct PrecisionEvaluator {
    family: ScalarIntegral,
    decimal_digits: u32,
    binary_precision: u32,
    evaluator: PreparedPrecision,
    arguments: Vec<C>,
}

impl PrecisionEvaluator {
    /// Compile the native expressions for a positive decimal-digit request.
    /// Ordinary Python-number calls may use the separate f64 JIT API at its
    /// default 16 digits; this constructor always selects arbitrary precision.
    pub fn new(family: ScalarIntegral, decimal_digits: u32) -> Result<Self, String> {
        Self::with_backend(family, decimal_digits, default_backend())
    }

    /// Choose direct Rust or native-expression interpretation explicitly.
    /// SymJIT is binary64-only and cannot provide arbitrary precision.
    pub fn with_backend(
        family: ScalarIntegral,
        decimal_digits: u32,
        backend: EvaluationBackend,
    ) -> Result<Self, String> {
        if decimal_digits == 0 {
            return Err("decimal precision must be positive".into());
        }
        // A conservative rational upper bound for log2(10), avoiding an
        // overflow-prone float-to-integer conversion for very large requests.
        const LOG2_10_NUMERATOR: u64 = 3_321_928_095;
        const DENOMINATOR: u64 = 1_000_000_000;
        let bits = (u64::from(decimal_digits) * LOG2_10_NUMERATOR).div_ceil(DENOMINATOR) + 32;
        let bits = u32::try_from(bits)
            .map_err(|_| "decimal precision exceeds the binary precision range")?;
        Self::build(family, decimal_digits, bits, backend)
    }

    /// Compile at exactly `bits` working bits, without adding decimal guards.
    /// `precision()` reports the corresponding number of whole decimal digits;
    /// `binary_precision()` reports the supplied working precision exactly.
    pub fn with_binary_precision(family: ScalarIntegral, bits: u32) -> Result<Self, String> {
        Self::with_binary_precision_and_backend(family, bits, default_backend())
    }

    /// Specify working bits and implementation independently; no guard bits are added.
    pub fn with_binary_precision_and_backend(
        family: ScalarIntegral,
        bits: u32,
        backend: EvaluationBackend,
    ) -> Result<Self, String> {
        if bits == 0 {
            return Err("binary precision must be positive".into());
        }
        // This conversion is metadata only; no numerical input or coefficient
        // is ever converted to f64 by this evaluator.
        let decimal_digits = (f64::from(bits) * std::f64::consts::LOG10_2).floor() as u32;
        Self::build(family, decimal_digits, bits, backend)
    }

    fn build(
        family: ScalarIntegral,
        decimal_digits: u32,
        bits: u32,
        backend: EvaluationBackend,
    ) -> Result<Self, String> {
        if backend == EvaluationBackend::SymJit {
            return Err(
                "SymJIT supports binary64 only; use native or expression for arbitrary precision"
                    .into(),
            );
        }
        crate::initialization::ensure_symbolica_state();
        let converter = Float::new(bits);
        let evaluator = if backend == EvaluationBackend::Native {
            PreparedPrecision::Native(NativeEvaluator::with_binary_precision(family, bits)?)
        } else {
            PreparedPrecision::Expression(
                masters::exact_evaluator(family)
                    .clone()
                    .map_coeff_with_prec(
                        &|value| {
                            Complex::new(
                                converter.from_rational(&value.re),
                                converter.from_rational(&value.im),
                            )
                        },
                        bits,
                    ),
            )
        };
        Ok(Self {
            family,
            decimal_digits,
            binary_precision: bits,
            evaluator,
            arguments: (0..family.arity())
                .map(|_| Complex::new(Float::new(bits), Float::new(bits)))
                .collect(),
        })
    }

    /// Requested decimal digits, or whole decimal digits for the binary constructor.
    pub fn precision(&self) -> u32 {
        self.decimal_digits
    }

    /// Fixed working precision, including decimal-constructor guard bits.
    pub fn binary_precision(&self) -> u32 {
        self.binary_precision
    }

    /// Integral family represented by this evaluator.
    pub fn family(&self) -> ScalarIntegral {
        self.family
    }

    /// Implementation used by this fixed-precision workspace.
    pub fn backend(&self) -> EvaluationBackend {
        match self.evaluator {
            PreparedPrecision::Expression(_) => EvaluationBackend::Expression,
            PreparedPrecision::Native(_) => EvaluationBackend::Native,
        }
    }

    /// Evaluate one point into `[finite, simple pole, double pole]`.
    /// Inputs use the scalar family's usual ordering, with `mu_squared` last.
    /// Inputs must be finite, external squared invariants real, squared masses
    /// in the closed lower half-plane, and `mu_squared` positive and real.
    /// Validation failures leave the caller's output unchanged. Singular
    /// kinematics may still produce nonfinite coefficients.
    pub fn evaluate(&mut self, arguments: &[C], output: &mut [C]) -> Result<(), String> {
        crate::record_usage();
        self.validate(arguments.len(), output.len(), 1)?;
        self.validate_point(arguments)?;
        self.evaluate_row(arguments, output);
        Ok(())
    }

    /// Evaluate flat row-major points, writing three consecutive coefficients
    /// for each row. Empty batches are accepted without entering the evaluator.
    /// This native arbitrary-precision path is scalar, not a SIMD/f64 fallback.
    pub fn evaluate_batch(
        &mut self,
        arguments: &[C],
        output: &mut [C],
        rows: usize,
    ) -> Result<(), String> {
        crate::record_usage();
        self.validate(arguments.len(), output.len(), rows)?;
        // Validate the entire batch before evaluating any row, and before
        // rounding inputs to the chosen working precision.
        for input in arguments.chunks_exact(self.family.arity()) {
            self.validate_point(input)?;
        }
        for (input, output) in arguments
            .chunks_exact(self.family.arity())
            .zip(output.chunks_exact_mut(3))
        {
            self.evaluate_row(input, output);
        }
        Ok(())
    }

    fn validate(&self, arguments: usize, output: usize, rows: usize) -> Result<(), String> {
        let expected_arguments = rows
            .checked_mul(self.family.arity())
            .ok_or("input size overflow")?;
        let expected_output = rows.checked_mul(3).ok_or("output size overflow")?;
        if arguments != expected_arguments || output != expected_output {
            return Err(format!(
                "expected {expected_arguments} inputs and {expected_output} outputs for {rows} rows; received {arguments} and {output}"
            ));
        }
        Ok(())
    }

    fn evaluate_row(&mut self, arguments: &[C], output: &mut [C]) {
        match &mut self.evaluator {
            PreparedPrecision::Expression(evaluator) => {
                for (working, supplied) in self.arguments.iter_mut().zip(arguments) {
                    working.clone_from(supplied);
                    if working.re.prec() != self.binary_precision {
                        working.re.set_prec(self.binary_precision);
                    }
                    if working.im.prec() != self.binary_precision {
                        working.im.set_prec(self.binary_precision);
                    }
                }
                if self.family == ScalarIntegral::D0 {
                    crate::labelling::regular_d0_labelling(&mut self.arguments);
                }
                evaluator.evaluate(&self.arguments, output);
            }
            // NativeEvaluator already maintains its own rounded-input workspace.
            // Pass the original values directly instead of cloning them twice.
            PreparedPrecision::Native(evaluator) => evaluator
                .evaluate(arguments, output)
                .expect("previously validated precise inputs"),
        }
        // In particular, do not set_prec upward on computed output values:
        // their precision records may reflect cancellation in the expression.
    }

    fn validate_point(&self, arguments: &[C]) -> Result<(), String> {
        if arguments
            .iter()
            .any(|value| !value.re.is_finite() || !value.im.is_finite())
        {
            return Err("all input components must be finite".into());
        }
        let momenta = match self.family {
            ScalarIntegral::A0 => 0,
            ScalarIntegral::B0 | ScalarIntegral::DB0 => 1,
            ScalarIntegral::C0 => 3,
            ScalarIntegral::D0 => 6,
        };
        if arguments[..momenta].iter().any(|value| !value.im.is_zero()) {
            return Err("external squared invariants must be real".into());
        }
        if arguments[momenta..arguments.len() - 1]
            .iter()
            .any(|value| !value.im.is_zero() && !value.im.is_negative())
        {
            return Err("squared masses require nonpositive imaginary parts".into());
        }
        let mu = &arguments[arguments.len() - 1];
        if !mu.im.is_zero() || mu.re.is_zero() || mu.re.is_negative() {
            return Err("mu_squared must be positive and real".into());
        }
        Ok(())
    }
}
