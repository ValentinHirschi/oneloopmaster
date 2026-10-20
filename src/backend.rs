//! Explicit selection between direct Rust and inspectable evaluator backends.
use crate::{NativeEvaluator, ScalarIntegral};
use std::sync::{Mutex, OnceLock};
use symbolica::{domains::float::Complex, evaluate::ExpressionEvaluator};

/// Numerical implementation selected for manual evaluations.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EvaluationBackend {
    /// Ahead-of-time generic Rust with no runtime expression evaluator.
    Native,
    /// The prepared binary64 SymJIT O2 evaluator (no arbitrary precision).
    SymJit,
    /// Symbolica's native expression interpreter, also supporting Float.
    Expression,
}

impl EvaluationBackend {
    /// Whether this implementation was included in the build.
    pub const fn is_available(self) -> bool {
        match self {
            Self::Native => cfg!(feature = "generated-evaluators"),
            Self::SymJit | Self::Expression => true,
        }
    }

    /// Stable name used by the Python adapter and diagnostics.
    pub fn name(self) -> &'static str {
        match self {
            Self::Native => "native",
            Self::SymJit => "symjit",
            Self::Expression => "expression",
        }
    }
}

// Direct Rust is faster on each measured mixed-family workload. SymJIT remains
// selectable explicitly; uniform SIMD-friendly inputs can still favor it.
#[cfg(feature = "generated-evaluators")]
pub const DEFAULT_BACKEND: EvaluationBackend = EvaluationBackend::Native;
#[cfg(not(feature = "generated-evaluators"))]
pub const DEFAULT_BACKEND: EvaluationBackend = EvaluationBackend::Expression;

type C = Complex<f64>;
struct NativeBackend {
    evaluator: NativeEvaluator<f64>,
    group: crate::evaluators::CoefficientGroups,
}
static NATIVE: [OnceLock<Mutex<NativeBackend>>; 5] = [const { OnceLock::new() }; 5];
struct ExpressionBackend {
    evaluator: ExpressionEvaluator<C>,
    group: crate::evaluators::CoefficientGroups,
}
static EXPRESSION: [OnceLock<Mutex<ExpressionBackend>>; 5] = [const { OnceLock::new() }; 5];

fn expression(family: ScalarIntegral) -> &'static Mutex<ExpressionBackend> {
    EXPRESSION[family as usize].get_or_init(|| {
        Mutex::new(ExpressionBackend {
            evaluator: crate::masters::exact_evaluator(family)
                .clone()
                .map_coeff(&|v| C::new(v.re.to_f64(), v.im.to_f64())),
            group: crate::evaluators::CoefficientGroups::new(),
        })
    })
}

#[cfg(not(feature = "generated-evaluators"))]
pub(crate) fn expression_coefficient(family: ScalarIntegral, tag: usize, input: &[C]) -> C {
    crate::native::validate_point(family, input).expect("valid master arguments");
    let mut value = expression(family)
        .lock()
        .expect("expression evaluator cache poisoned");
    let ExpressionBackend { evaluator, group } = &mut *value;
    group.evaluate(input, tag, |output| {
        let mut relabelled = Vec::new();
        evaluator.evaluate(crate::labelling::regular_point(input, &mut relabelled), output);
    })
}

fn native(family: ScalarIntegral) -> Result<&'static Mutex<NativeBackend>, String> {
    if let Some(value) = NATIVE[family as usize].get() {
        return Ok(value);
    }
    let value = Mutex::new(NativeBackend {
        evaluator: NativeEvaluator::new(family)?,
        group: crate::evaluators::CoefficientGroups::new(),
    });
    // Concurrent creation is harmless; constants and workspaces are independent.
    let _ = NATIVE[family as usize].set(value);
    Ok(NATIVE[family as usize].get().expect("native cache was set"))
}

#[cfg(feature = "generated-evaluators")]
pub(crate) fn initialize_native_all() -> Result<(), String> {
    for family in [
        ScalarIntegral::A0,
        ScalarIntegral::B0,
        ScalarIntegral::DB0,
        ScalarIntegral::C0,
        ScalarIntegral::D0,
    ] {
        let _ = native(family)?;
    }
    Ok(())
}

/// Evaluate a point through the selected default implementation.
pub fn evaluate(family: ScalarIntegral, input: &[C], output: &mut [C]) -> Result<(), String> {
    evaluate_with_backend(family, input, output, DEFAULT_BACKEND)
}

/// Evaluate flat row-major points through the selected default implementation.
pub fn evaluate_batch(
    family: ScalarIntegral,
    input: &[C],
    output: &mut [C],
    rows: usize,
) -> Result<(), String> {
    evaluate_batch_with_backend(family, input, output, rows, DEFAULT_BACKEND)
}

/// Explicit manual selection, retained for performance comparisons and tests.
pub fn evaluate_with_backend(
    family: ScalarIntegral,
    input: &[C],
    output: &mut [C],
    backend: EvaluationBackend,
) -> Result<(), String> {
    crate::record_usage();
    evaluate_batch_with_backend(family, input, output, 1, backend)
}

/// Explicit backend selection for a flat row-major batch.
pub fn evaluate_batch_with_backend(
    family: ScalarIntegral,
    input: &[C],
    output: &mut [C],
    rows: usize,
    backend: EvaluationBackend,
) -> Result<(), String> {
    crate::record_usage();
    if backend != EvaluationBackend::Native {
        let expected_input = rows
            .checked_mul(family.arity())
            .ok_or("input size overflow")?;
        let expected_output = rows.checked_mul(3).ok_or("output size overflow")?;
        if input.len() != expected_input || output.len() != expected_output {
            return Err(format!(
                "expected {expected_input} inputs and {expected_output} outputs"
            ));
        }
        for point in input.chunks_exact(family.arity()) {
            crate::native::validate_point(family, point)?;
        }
    }
    match backend {
        EvaluationBackend::Native => {
            crate::initialization::ensure_symbolica_state();
            let mut value = native(family)?
                .lock()
                .map_err(|_| "native evaluator cache poisoned")?;
            value.group.invalidate();
            value.evaluator.evaluate_batch(input, output, rows)
        }
        EvaluationBackend::SymJit => crate::evaluators::evaluate_batch(family, input, output, rows),
        EvaluationBackend::Expression => {
            let mut value = expression(family)
                .lock()
                .map_err(|_| "expression evaluator cache poisoned")?;
            value.group.invalidate();
            let mut relabelled = Vec::new();
            for (point, values) in input
                .chunks_exact(family.arity())
                .zip(output.chunks_exact_mut(3))
            {
                let point = crate::labelling::regular_point(point, &mut relabelled);
                value.evaluator.evaluate(point, values);
            }
            Ok(())
        }
    }
}

#[cfg(feature = "generated-evaluators")]
pub(crate) fn native_coefficient(family: ScalarIntegral, tag: usize, input: &[C]) -> C {
    let (coefficient, result) = {
        let mut value = native(family)
            .expect("native constants initialized")
            .lock()
            .expect("native evaluator cache poisoned");
        let NativeBackend { evaluator, group } = &mut *value;
        let mut result = Ok(());
        let coefficient = group.evaluate(input, tag, |output| {
            result = evaluator.evaluate(input, output);
        });
        if result.is_err() {
            // The group callback has no Result channel. Do not retain the
            // placeholder output it installed after a rejected point.
            group.invalidate();
        }
        (coefficient, result)
    };
    // Invalid user input is an expected error, not a corrupt evaluator. Keep
    // the numeric callback's panic API without poisoning the shared mutex.
    result.expect("valid master arguments");
    coefficient
}
