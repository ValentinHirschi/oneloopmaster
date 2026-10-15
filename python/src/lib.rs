//! Thin Python adapters for the core native scalar evaluator.
//!
//! Standalone wheels expose numbers only. Expression interop is enabled only
//! through the `community` feature, linked into one Symbolica community kernel.
//! Identical revisions in separate extension binaries do not imply shared state.
//! Import registers the Python API only. Symbol hooks are registered with
//! Symbolica state; formulas and numerical backends are prepared on demand.
//! Community uppercase exports are lazily resolved primitive Expressions;
//! lowercase functions evaluate numerical coefficient triples in both builds.
//! Construction/evaluation remains on the calling thread; no hidden worker or
//! automatic precision escalation is introduced here.

mod precision;

use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarEvaluator, ScalarIntegral,
};
use pyo3::{
    exceptions::{PyRuntimeError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyAnyMethods, PyComplex, PyComplexMethods, PyModuleMethods},
};
use symbolica::{
    atom::Atom,
    domains::float::{Complex, Float, SingleFloat},
    evaluate::ExpressionEvaluator,
};

type Coefficients = (Py<PyAny>, Py<PyAny>, Py<PyAny>);

#[derive(Clone, Copy)]
enum BackendChoice {
    Auto,
    Native,
    SymJit,
    Expression,
}

impl BackendChoice {
    fn parse(name: &str) -> PyResult<Self> {
        match name.trim().to_ascii_lowercase().as_str() {
            "auto" => Ok(Self::Auto),
            "native" => Ok(Self::Native),
            "symjit" => Ok(Self::SymJit),
            "expression" | "symbolica" => Ok(Self::Expression),
            _ => Err(PyValueError::new_err(
                "backend must be auto, native, symjit or expression",
            )),
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Native => "native",
            Self::SymJit => "symjit",
            Self::Expression => "expression",
        }
    }

    fn resolve(self, arbitrary: bool) -> PyResult<EvaluationBackend> {
        let backend = match self {
            Self::Auto if arbitrary => match oneloop::DEFAULT_BACKEND {
                EvaluationBackend::Native => EvaluationBackend::Native,
                _ => EvaluationBackend::Expression,
            },
            Self::Auto => oneloop::DEFAULT_BACKEND,
            Self::Native => EvaluationBackend::Native,
            Self::SymJit => EvaluationBackend::SymJit,
            Self::Expression => EvaluationBackend::Expression,
        };
        if !backend.is_available() {
            return Err(PyValueError::new_err(
                "backend='native' requires the generated-evaluators Cargo feature; use auto or expression, or rebuild with that feature",
            ));
        }
        if arbitrary && backend == EvaluationBackend::SymJit {
            return Err(PyValueError::new_err(
                "backend='symjit' supports binary64 only; Decimal/large-integer inputs or prec != 16 require native, expression or auto",
            ));
        }
        Ok(backend)
    }
}

fn backend_name(backend: EvaluationBackend) -> &'static str {
    match backend {
        EvaluationBackend::Native => "native",
        EvaluationBackend::SymJit => "symjit",
        EvaluationBackend::Expression => "expression",
    }
}

enum MachineEvaluator {
    Native(NativeEvaluator<f64>),
    SymJit(ScalarEvaluator),
    Expression(ExpressionEvaluator<Complex<f64>>),
}

impl MachineEvaluator {
    fn new(family: ScalarIntegral, backend: EvaluationBackend, rebuild: bool) -> PyResult<Self> {
        oneloop::record_usage();
        match backend {
            EvaluationBackend::Native => NativeEvaluator::new(family)
                .map(Self::Native)
                .map_err(PyRuntimeError::new_err),
            EvaluationBackend::SymJit => build(family, rebuild).map(Self::SymJit),
            EvaluationBackend::Expression => {
                let parameters = (0..family.arity())
                    .map(|index| {
                        symbolica::symbol!(format!(
                            "__oneloop_python_expression_{}_{}",
                            family.name(),
                            index
                        ))
                        .to_atom()
                    })
                    .collect::<Vec<_>>();
                let master = match family {
                    ScalarIntegral::A0 => oneloop::A0(),
                    ScalarIntegral::B0 => oneloop::B0(),
                    ScalarIntegral::DB0 => oneloop::dB0(),
                    ScalarIntegral::C0 => oneloop::C0(),
                    ScalarIntegral::D0 => oneloop::D0(),
                };
                let calls = [0, -1, -2].map(|tag| {
                    let mut arguments = vec![Atom::num(tag)];
                    arguments.extend(parameters.iter().cloned());
                    master.call(&arguments)
                });
                let exact = oneloop::OneLoopExpressions::new()
                    .evaluator(&calls, &parameters)
                    .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
                Ok(Self::Expression(exact.map_coeff(&|value| {
                    Complex::new(value.re.to_f64(), value.im.to_f64())
                })))
            }
        }
    }

    fn evaluate(&mut self, input: &[Complex<f64>], output: &mut [Complex<f64>]) -> PyResult<()> {
        match self {
            Self::Native(evaluator) => evaluator
                .evaluate(input, output)
                .map_err(PyRuntimeError::new_err),
            Self::SymJit(evaluator) => evaluator
                .evaluate(input, output)
                .map_err(PyRuntimeError::new_err),
            Self::Expression(evaluator) => {
                evaluator.evaluate(input, output);
                Ok(())
            }
        }
    }

    fn evaluate_batch(
        &mut self,
        input: &[Complex<f64>],
        output: &mut [Complex<f64>],
        rows: usize,
        arity: usize,
    ) -> PyResult<()> {
        match self {
            Self::Native(evaluator) => evaluator
                .evaluate_batch(input, output, rows)
                .map_err(PyRuntimeError::new_err),
            Self::SymJit(evaluator) => evaluator
                .evaluate_batch(input, output, rows)
                .map_err(PyRuntimeError::new_err),
            Self::Expression(evaluator) => {
                for (point, result) in input.chunks_exact(arity).zip(output.chunks_exact_mut(3)) {
                    evaluator.evaluate(point, result);
                }
                Ok(())
            }
        }
    }
}

fn family(name: &str) -> PyResult<ScalarIntegral> {
    match name.trim().to_ascii_uppercase().as_str() {
        "A0" => Ok(ScalarIntegral::A0),
        "B0" => Ok(ScalarIntegral::B0),
        "DB0" => Ok(ScalarIntegral::DB0),
        "C0" => Ok(ScalarIntegral::C0),
        "D0" => Ok(ScalarIntegral::D0),
        _ => Err(PyValueError::new_err(
            "family must be A0, B0, dB0, C0 or D0",
        )),
    }
}

/// A symbolic family selector belongs to the same kernel as this module.
#[cfg(feature = "community")]
fn symbolic_family(value: &Bound<'_, PyAny>) -> PyResult<ScalarIntegral> {
    // Symbolica's extractor accepts both Symbol and a variable Expression,
    // without accepting strings or casting a composite into a made-up name.
    let symbol = value.extract::<symbolica::atom::Symbol>()?;
    let name = symbol.get_name();
    let short = name
        .strip_prefix("oneloopmaster::")
        .or_else(|| name.strip_prefix("oneloop::"))
        .ok_or_else(|| {
            PyValueError::new_err("expected an oneloopmaster/oneloop scalar master symbol")
        })?;
    match short {
        "A0" | "B0" | "dB0" | "C0" | "D0" => family(short),
        _ => Err(PyValueError::new_err("expected A0, B0, dB0, C0 or D0")),
    }
}

fn number(value: &Bound<'_, PyAny>) -> PyResult<Complex<f64>> {
    if let Ok(value) = value.cast::<PyComplex>() {
        Ok(Complex::new(value.real(), value.imag()))
    } else if let Ok(value) = value.extract::<f64>() {
        Ok(Complex::new(value, 0.0))
    } else {
        Err(PyTypeError::new_err(
            "numeric inputs must be Python real or complex numbers",
        ))
    }
}

fn validate(family: ScalarIntegral, input: &[Complex<f64>]) -> PyResult<()> {
    if input.len() != family.arity() {
        return Err(PyValueError::new_err(format!(
            "{} expects {} ordered arguments, including mu_squared last; received {}",
            family.name(),
            family.arity(),
            input.len()
        )));
    }
    if input.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
        return Err(PyValueError::new_err("all input components must be finite"));
    }
    let momenta = match family {
        ScalarIntegral::A0 => 0,
        ScalarIntegral::B0 | ScalarIntegral::DB0 => 1,
        ScalarIntegral::C0 => 3,
        ScalarIntegral::D0 => 6,
    };
    if input[..momenta].iter().any(|z| z.im != 0.0) {
        return Err(PyValueError::new_err(
            "external squared invariants must be real",
        ));
    }
    if input[momenta..input.len() - 1].iter().any(|z| z.im > 0.0) {
        return Err(PyValueError::new_err(
            "squared masses require nonpositive imaginary parts",
        ));
    }
    let mu = input[input.len() - 1];
    if mu.im != 0.0 || mu.re <= 0.0 {
        return Err(PyValueError::new_err(
            "mu_squared must be positive and real",
        ));
    }
    Ok(())
}

fn coefficients(py: Python<'_>, output: &[Complex<f64>]) -> Coefficients {
    let convert = |z: &Complex<f64>| PyComplex::from_doubles(py, z.re, z.im).unbind().into_any();
    (
        convert(&output[0]),
        convert(&output[1]),
        convert(&output[2]),
    )
}

fn arbitrary_required(arguments: &[Bound<'_, PyAny>], digits: u32) -> PyResult<bool> {
    let mut required = digits != 16;
    for value in arguments {
        // Validate every input type, even after a Decimal has selected this path.
        required |= precision::needs_arbitrary(value)?;
    }
    Ok(required)
}

fn validate_arbitrary(family: ScalarIntegral, input: &[Complex<Float>]) -> PyResult<()> {
    if input.len() != family.arity() {
        return Err(PyValueError::new_err(format!(
            "{} expects {} ordered arguments, including mu_squared last; received {}",
            family.name(),
            family.arity(),
            input.len()
        )));
    }
    if input.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
        return Err(PyValueError::new_err("all input components must be finite"));
    }
    let momenta = match family {
        ScalarIntegral::A0 => 0,
        ScalarIntegral::B0 | ScalarIntegral::DB0 => 1,
        ScalarIntegral::C0 => 3,
        ScalarIntegral::D0 => 6,
    };
    if input[..momenta].iter().any(|z| !z.im.is_zero()) {
        return Err(PyValueError::new_err(
            "external squared invariants must be real",
        ));
    }
    if input[momenta..input.len() - 1]
        .iter()
        .any(|z| !z.im.is_zero() && !z.im.is_negative())
    {
        return Err(PyValueError::new_err(
            "squared masses require nonpositive imaginary parts",
        ));
    }
    let mu = &input[input.len() - 1];
    if !mu.im.is_zero() || mu.re.is_zero() || mu.re.is_negative() {
        return Err(PyValueError::new_err(
            "mu_squared must be positive and real",
        ));
    }
    Ok(())
}

fn arbitrary_coefficients(
    py: Python<'_>,
    output: &[Complex<Float>],
    digits: u32,
) -> PyResult<Coefficients> {
    Ok((
        precision::decimal_output(py, &output[0], digits)?,
        precision::decimal_output(py, &output[1], digits)?,
        precision::decimal_output(py, &output[2], digits)?,
    ))
}

fn arbitrary_single(
    py: Python<'_>,
    evaluator: &mut PrecisionEvaluator,
    arguments: &[Bound<'_, PyAny>],
    digits: u32,
) -> PyResult<Coefficients> {
    let bits = evaluator.binary_precision();
    let input = arguments
        .iter()
        .map(|value| precision::arbitrary_number(value, bits))
        .collect::<PyResult<Vec<_>>>()?;
    validate_arbitrary(evaluator.family(), &input)?;
    let mut output =
        core::array::from_fn::<_, 3, _>(|_| Complex::new(Float::new(bits), Float::new(bits)));
    evaluator
        .evaluate(&input, &mut output)
        .map_err(PyRuntimeError::new_err)?;
    arbitrary_coefficients(py, &output, digits)
}

fn single_inputs(
    py: Python<'_>,
    family: ScalarIntegral,
    arguments: &[Bound<'_, PyAny>],
    digits: u32,
    rebuild: bool,
    backend: BackendChoice,
) -> PyResult<Coefficients> {
    if arbitrary_required(arguments, digits)? {
        let mut evaluator =
            PrecisionEvaluator::with_backend(family, digits, backend.resolve(true)?)
                .map_err(PyValueError::new_err)?;
        arbitrary_single(py, &mut evaluator, arguments, digits)
    } else {
        let input = arguments.iter().map(number).collect::<PyResult<Vec<_>>>()?;
        single(py, family, &input, rebuild, backend.resolve(false)?)
    }
}

fn build(family: ScalarIntegral, rebuild: bool) -> PyResult<ScalarEvaluator> {
    if rebuild {
        ScalarEvaluator::rebuild(family)
    } else {
        ScalarEvaluator::cached(family)
    }
    .map_err(PyRuntimeError::new_err)
}

fn single(
    py: Python<'_>,
    family: ScalarIntegral,
    input: &[Complex<f64>],
    rebuild: bool,
    backend: EvaluationBackend,
) -> PyResult<Coefficients> {
    validate(family, input)?;
    let mut output = [Complex::new(0.0, 0.0); 3];
    if rebuild {
        MachineEvaluator::new(family, backend, true)?.evaluate(input, &mut output)?;
    } else {
        oneloop::evaluate_with_backend(family, input, &mut output, backend)
            .map_err(PyRuntimeError::new_err)?;
    }
    Ok(coefficients(py, &output))
}

/// Report whether the native one-loop module has initialized its evaluation support.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hepkit as hep
/// >>> from symbolica.community.hepkit import oneloop
/// >>> initialized = oneloop.is_initialized()
#[pyfunction]
fn is_initialized() -> bool {
    oneloop::is_initialized()
}

/// Reuse a numerical evaluator for one scalar one-loop master family.
///
/// Pass a bare primitive symbol: ``oneloop.A0``, ``B0``, ``dB0``, ``C0`` or
/// ``D0``. Every evaluation takes physical arguments in that primitive's order,
/// with the squared renormalization scale last. Unlike the lowercase convenience
/// functions, ``evaluate`` requires the scale explicitly.
///
/// Results are ordered ``(finite, 1/eps, 1/eps**2)``. ``prec`` counts decimal
/// significant digits. Binary64 evaluation returns Python complex values;
/// arbitrary-precision evaluation returns DecimalComplex components. Use decimal
/// strings through Decimal to avoid rounding your inputs before evaluation.
/// Per-call precision/backend overrides do not change the stored defaults.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hepkit as hep
/// >>> from symbolica.community.hepkit import oneloop
/// >>> evaluator = oneloop.Evaluator(oneloop.A0)
/// >>> finite, pole, double_pole = evaluator.evaluate([1.0, 1.0])
/// >>> assert pole == 1+0j and double_pole == 0j
/// >>> rows = evaluator.evaluate_batch([[1.0, 1.0], [2.0, 1.0]])
/// >>> assert len(rows) == 2
#[cfg_attr(
    feature = "community",
    pyclass(
        name = "Evaluator",
        unsendable,
        module = "symbolica.community.hepkit.oneloop"
    )
)]
#[cfg_attr(
    not(feature = "community"),
    pyclass(name = "Evaluator", unsendable, module = "oneloop_native")
)]
struct Evaluator {
    family: ScalarIntegral,
    evaluator: Option<(EvaluationBackend, MachineEvaluator)>,
    prec: u32,
    backend: BackendChoice,
    arbitrary: Option<(u32, PrecisionEvaluator)>,
}

impl Evaluator {
    fn machine(&mut self, backend: EvaluationBackend) -> PyResult<&mut MachineEvaluator> {
        if self
            .evaluator
            .as_ref()
            .is_none_or(|(current, _)| *current != backend)
        {
            self.evaluator = Some((backend, MachineEvaluator::new(self.family, backend, false)?));
        }
        Ok(&mut self.evaluator.as_mut().unwrap().1)
    }

    fn arbitrary(
        &mut self,
        digits: u32,
        backend: EvaluationBackend,
    ) -> PyResult<&mut PrecisionEvaluator> {
        if self
            .arbitrary
            .as_ref()
            .is_none_or(|(current, value)| *current != digits || value.backend() != backend)
        {
            let evaluator = PrecisionEvaluator::with_backend(self.family, digits, backend)
                .map_err(PyValueError::new_err)?;
            self.arbitrary = Some((digits, evaluator));
        }
        Ok(&mut self.arbitrary.as_mut().unwrap().1)
    }
}

#[pymethods]
impl Evaluator {
    /// Prepare an evaluator and retain default precision and backend choices.
    ///
    /// ``auto`` selects a supported backend for the requested precision. ``native``
    /// uses native numerical evaluation; ``symjit`` supports binary64 only.
    /// ``expression`` (also named ``symbolica``) evaluates Symbolica formulas.
    /// Unsupported backend/precision combinations raise an error.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hepkit as hep
    /// >>> from symbolica.community.hepkit import oneloop
    /// >>> evaluator = oneloop.Evaluator(oneloop.A0)
    /// >>> assert evaluator.family == "A0" and evaluator.arity == 2
    ///
    /// Parameters
    /// ----------
    /// family : Expression
    ///     Bare primitive symbol such as oneloop.B0, not a B0(...) call.
    /// rebuild : bool, optional
    ///     Recreate the evaluator's workspace; default False.
    /// prec : int, optional
    ///     Positive number of decimal significant digits; default 16.
    /// backend : str, optional
    ///     "auto", "native", "symjit", "expression" or "symbolica"; default "auto".
    #[new]
    #[pyo3(signature = (family, rebuild=false, *, prec=16, backend="auto"))]
    fn new(
        family: &Bound<'_, PyAny>,
        rebuild: bool,
        #[pyo3(from_py_with = precision::parse_precision)] prec: u32,
        backend: &str,
    ) -> PyResult<Self> {
        oneloop::record_usage();
        #[cfg(feature = "community")]
        let family = symbolic_family(family)?;
        #[cfg(not(feature = "community"))]
        let family = self::family(&family.extract::<String>()?)?;
        let backend = BackendChoice::parse(backend)?;
        let resolved = backend.resolve(prec != 16)?;
        Ok(Self {
            family,
            evaluator: if prec == 16 {
                Some((resolved, MachineEvaluator::new(family, resolved, rebuild)?))
            } else {
                None
            },
            prec,
            backend,
            arbitrary: if prec != 16 {
                Some((
                    prec,
                    PrecisionEvaluator::with_backend(family, prec, resolved)
                        .map_err(PyValueError::new_err)?,
                ))
            } else {
                None
            },
        })
    }

    /// Primitive family name, such as "A0" or "B0".
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hepkit as hep
    /// >>> from symbolica.community.hepkit import oneloop
    /// >>> evaluator = oneloop.Evaluator(oneloop.A0)
    /// >>> assert evaluator.family == "A0"
    #[getter]
    fn family(&self) -> &'static str {
        self.family.name()
    }

    /// Number of physical arguments including the final squared scale.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hepkit as hep
    /// >>> from symbolica.community.hepkit import oneloop
    /// >>> evaluator = oneloop.Evaluator(oneloop.A0)
    /// >>> assert evaluator.arity == 2
    /// >>> assert oneloop.Evaluator(oneloop.B0).arity == 4
    #[getter]
    fn arity(&self) -> usize {
        self.family.arity()
    }

    /// Default decimal significant-digit count; per-call overrides leave it unchanged.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hepkit as hep
    /// >>> from symbolica.community.hepkit import oneloop
    /// >>> evaluator = oneloop.Evaluator(oneloop.A0)
    /// >>> assert evaluator.prec == 16
    #[getter]
    fn prec(&self) -> u32 {
        self.prec
    }

    /// Requested default backend name; "auto" remains "auto" after backend selection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hepkit as hep
    /// >>> from symbolica.community.hepkit import oneloop
    /// >>> evaluator = oneloop.Evaluator(oneloop.A0)
    /// >>> assert evaluator.backend == "auto"
    #[getter]
    fn backend(&self) -> &'static str {
        self.backend.name()
    }

    /// Evaluate one kinematic point and return (finite, simple pole, double pole).
    ///
    /// Supply exactly ``arity`` numeric arguments including the squared scale.
    /// Inputs must satisfy the selected primitive's kinematic domain. A per-call
    /// ``prec`` or ``backend`` override changes only this evaluation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hepkit as hep
    /// >>> from symbolica.community.hepkit import oneloop
    /// >>> evaluator = oneloop.Evaluator(oneloop.A0)
    /// >>> finite, pole, double_pole = evaluator.evaluate([1.0, 1.0])
    /// >>> assert pole == 1+0j
    /// >>> from decimal import Decimal
    /// >>> high_precision = evaluator.evaluate([Decimal("1"), Decimal("1")], prec=40)
    /// >>> assert high_precision[1].real == Decimal("1")
    ///
    /// Parameters
    /// ----------
    /// arguments : sequence[Number]
    ///     Ordered invariants, squared masses and squared scale for the primitive.
    /// prec : int or None, optional
    ///     Decimal significant digits; None uses the constructor default.
    /// backend : str or None, optional
    ///     Backend override; None uses the constructor default.
    #[pyo3(signature = (arguments, *, prec=None, backend=None))]
    fn evaluate(
        &mut self,
        py: Python<'_>,
        arguments: Vec<Bound<'_, PyAny>>,
        prec: Option<&Bound<'_, PyAny>>,
        backend: Option<&str>,
    ) -> PyResult<Coefficients> {
        let digits = precision::requested_precision(prec, self.prec)?;
        let backend = backend
            .map(BackendChoice::parse)
            .transpose()?
            .unwrap_or(self.backend);
        if arbitrary_required(&arguments, digits)? {
            return arbitrary_single(
                py,
                self.arbitrary(digits, backend.resolve(true)?)?,
                &arguments,
                digits,
            );
        }
        let input = arguments.iter().map(number).collect::<PyResult<Vec<_>>>()?;
        validate(self.family, &input)?;
        let mut output = [Complex::new(0.0, 0.0); 3];
        self.machine(backend.resolve(false)?)?
            .evaluate(&input, &mut output)?;
        Ok(coefficients(py, &output))
    }

    /// Evaluate a sequence of kinematic rows, returning one coefficient tuple per row.
    ///
    /// All rows have ``arity`` entries and use the same precision/backend choice.
    /// An empty batch returns an empty list. No NumPy array is required.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hepkit as hep
    /// >>> from symbolica.community.hepkit import oneloop
    /// >>> evaluator = oneloop.Evaluator(oneloop.A0)
    /// >>> rows = evaluator.evaluate_batch([[1.0, 1.0], [2.0, 1.0]])
    /// >>> assert rows[0] == evaluator.evaluate([1.0, 1.0])
    /// >>> assert evaluator.evaluate_batch([]) == []
    ///
    /// Parameters
    /// ----------
    /// rows : sequence[sequence[Number]]
    ///     Ordered physical arguments, including the squared scale, for each point.
    /// prec : int or None, optional
    ///     Decimal significant digits; None uses the constructor default.
    /// backend : str or None, optional
    ///     Backend override; None uses the constructor default.
    #[pyo3(signature = (rows, *, prec=None, backend=None))]
    fn evaluate_batch(
        &mut self,
        py: Python<'_>,
        rows: Vec<Vec<Bound<'_, PyAny>>>,
        prec: Option<&Bound<'_, PyAny>>,
        backend: Option<&str>,
    ) -> PyResult<Vec<Coefficients>> {
        let digits = precision::requested_precision(prec, self.prec)?;
        let backend = backend
            .map(BackendChoice::parse)
            .transpose()?
            .unwrap_or(self.backend);
        let count = rows.len();
        if count == 0 {
            backend.resolve(digits != 16)?;
            return Ok(Vec::new());
        }
        let input_count = count
            .checked_mul(self.family.arity())
            .ok_or_else(|| PyValueError::new_err("input batch size overflow"))?;
        let output_count = count
            .checked_mul(3)
            .ok_or_else(|| PyValueError::new_err("output batch size overflow"))?;
        let mut arbitrary = digits != 16;
        for row in &rows {
            arbitrary |= arbitrary_required(row, digits)?;
        }
        if arbitrary {
            let family = self.family;
            let evaluator = self.arbitrary(digits, backend.resolve(true)?)?;
            let bits = evaluator.binary_precision();
            let mut input = Vec::with_capacity(input_count);
            for row in &rows {
                let start = input.len();
                for value in row {
                    input.push(precision::arbitrary_number(value, bits)?);
                }
                validate_arbitrary(family, &input[start..])?;
            }
            let mut output = (0..output_count)
                .map(|_| Complex::new(Float::new(bits), Float::new(bits)))
                .collect::<Vec<_>>();
            evaluator
                .evaluate_batch(&input, &mut output, count)
                .map_err(PyRuntimeError::new_err)?;
            return output
                .chunks_exact(3)
                .map(|row| arbitrary_coefficients(py, row, digits))
                .collect();
        }
        let mut input = Vec::with_capacity(input_count);
        for row in rows {
            let start = input.len();
            for value in row {
                input.push(number(&value)?);
            }
            validate(self.family, &input[start..])?;
        }
        let mut output = vec![Complex::new(0.0, 0.0); output_count];
        let arity = self.family.arity();
        self.machine(backend.resolve(false)?)?
            .evaluate_batch(&input, &mut output, count, arity)?;
        Ok(output
            .chunks_exact(3)
            .map(|row| coefficients(py, row))
            .collect())
    }

    /// Recreate cached evaluator workspaces while preserving the configured family and defaults.
    ///
    /// Native evaluation refreshes constants/workspaces; a SymJIT evaluator is
    /// recompiled. This is usually unnecessary between evaluations at new points.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hepkit as hep
    /// >>> from symbolica.community.hepkit import oneloop
    /// >>> evaluator = oneloop.Evaluator(oneloop.A0)
    /// >>> before = evaluator.evaluate([1.0, 1.0])
    /// >>> evaluator.rebuild()
    /// >>> assert evaluator.evaluate([1.0, 1.0]) == before
    fn rebuild(&mut self) -> PyResult<()> {
        let evaluator = self
            .evaluator
            .as_ref()
            .map(|(backend, _)| {
                MachineEvaluator::new(self.family, *backend, true).map(|value| (*backend, value))
            })
            .transpose()?;
        let arbitrary = self
            .arbitrary
            .as_ref()
            .map(|(digits, value)| {
                PrecisionEvaluator::with_backend(self.family, *digits, value.backend())
                    .map(|value| (*digits, value))
            })
            .transpose()
            .map_err(PyRuntimeError::new_err)?;
        self.evaluator = evaluator;
        self.arbitrary = arbitrary;
        Ok(())
    }

    /// Summarize the selected family, precision and requested backend.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.hepkit import oneloop
    /// >>> summary = repr(oneloop.Evaluator(oneloop.A0))
    fn __repr__(&self) -> String {
        format!(
            "Evaluator('{}', arity={}, prec={}, backend='{}', coefficient_order=(0,-1,-2))",
            self.family.name(),
            self.family.arity(),
            self.prec,
            self.backend.name(),
        )
    }
}

macro_rules! scalar {
    ($rust:ident, $python:literal, $kind:ident, [$($argument:ident),+], $documentation:literal) => {
        #[pyfunction(name=$python, signature=($($argument,)+ mu_squared=None, *, rebuild=false, prec=16, backend="auto"))]
        #[doc = $documentation]
        // Preserve the conventional ordered scalar-integral Python signature.
        #[allow(clippy::too_many_arguments)]
        fn $rust(py: Python<'_>, $($argument: &Bound<'_, PyAny>,)+ mu_squared: Option<&Bound<'_, PyAny>>, rebuild: bool, #[pyo3(from_py_with = precision::parse_precision)] prec: u32, backend: &str) -> PyResult<Coefficients> {
            let default_scale = 1u8.into_pyobject(py)?.into_any();
            let mut arguments = vec![$($argument.clone()),+];
            arguments.push(mu_squared.cloned().unwrap_or(default_scale));
            single_inputs(py, ScalarIntegral::$kind, &arguments, prec, rebuild, BackendChoice::parse(backend)?)
        }
    };
}

scalar!(
    a0,
    "A0",
    A0,
    [mass_squared],
    r#"Evaluate a scalar tadpole A0(mass_squared, mu_squared).

Results are ordered (finite, 1/eps, 1/eps**2). Invariants and masses are
squared quantities; omitted mu_squared is exactly one. ``prec`` is the
positive decimal significant-digit count (default 16). Decimal inputs or
higher precision produce DecimalComplex values. ``backend`` selects "auto",
"native", "symjit" (binary64 only), "expression", or "symbolica".
``rebuild=True`` refreshes the evaluator workspace.

Examples
--------
>>> from symbolica import S, E
>>> from symbolica.community import hepkit as hep
>>> from symbolica.community.hepkit import oneloop
>>> finite, pole, double_pole = oneloop.a0(1.0, 1.0)

Parameters
----------
mass_squared : Number
    Squared mass of the tadpole propagator.
mu_squared : Number or None, optional
    Squared renormalization scale; None uses exactly one.
rebuild : bool, optional
    Refresh the evaluator workspace before evaluation; default False.
prec : int, optional
    Positive number of decimal significant digits; default 16. Use Decimal
    inputs to retain input digits beyond binary64 precision.
backend : {"auto", "native", "symjit", "expression", "symbolica"}, optional
    Evaluation backend; default "auto" selects a supported backend for the
    requested precision. "symjit" supports binary64 only; "symbolica" is
    an alias for "expression"."#
);
scalar!(
    b0,
    "B0",
    B0,
    [momentum_squared, mass_0_squared, mass_1_squared],
    r#"Evaluate a scalar bubble B0(momentum_squared, mass_0_squared, mass_1_squared, mu_squared).

Results are ordered (finite, 1/eps, 1/eps**2). Invariants and masses are
squared quantities; omitted mu_squared is exactly one. ``prec`` is the
positive decimal significant-digit count (default 16). Decimal inputs or
higher precision produce DecimalComplex values. ``backend`` selects "auto",
"native", "symjit" (binary64 only), "expression", or "symbolica".
``rebuild=True`` refreshes the evaluator workspace.

Examples
--------
>>> from symbolica import S, E
>>> from symbolica.community import hepkit as hep
>>> from symbolica.community.hepkit import oneloop
>>> finite, pole, double_pole = oneloop.b0(-1.0, 1.0, 1.0, 1.0)

Parameters
----------
momentum_squared : Number
    External momentum squared; negative values describe spacelike momentum.
mass_0_squared, mass_1_squared : Number
    Squared masses of the two propagators, in B0 argument order.
mu_squared : Number or None, optional
    Squared renormalization scale; None uses exactly one.
rebuild : bool, optional
    Refresh the evaluator workspace before evaluation; default False.
prec : int, optional
    Positive number of decimal significant digits; default 16. Use Decimal
    inputs to retain input digits beyond binary64 precision.
backend : {"auto", "native", "symjit", "expression", "symbolica"}, optional
    Evaluation backend; default "auto" selects a supported backend for the
    requested precision. "symjit" supports binary64 only; "symbolica" is
    an alias for "expression"."#
);
scalar!(
    db0,
    "dB0",
    DB0,
    [momentum_squared, mass_0_squared, mass_1_squared],
    r#"Evaluate the derivative of B0 with respect to its external momentum squared.

Results are ordered (finite, 1/eps, 1/eps**2). Invariants and masses are
squared quantities; omitted mu_squared is exactly one. ``prec`` is the
positive decimal significant-digit count (default 16). Decimal inputs or
higher precision produce DecimalComplex values. ``backend`` selects "auto",
"native", "symjit" (binary64 only), "expression", or "symbolica".
``rebuild=True`` refreshes the evaluator workspace.

Examples
--------
>>> from symbolica import S, E
>>> from symbolica.community import hepkit as hep
>>> from symbolica.community.hepkit import oneloop
>>> finite, pole, double_pole = oneloop.db0(-1.0, 1.0, 1.0, 1.0)

Parameters
----------
momentum_squared : Number
    External momentum squared at which the derivative of B0 is evaluated.
mass_0_squared, mass_1_squared : Number
    Squared masses held fixed while differentiating B0.
mu_squared : Number or None, optional
    Squared renormalization scale; None uses exactly one.
rebuild : bool, optional
    Refresh the evaluator workspace before evaluation; default False.
prec : int, optional
    Positive number of decimal significant digits; default 16. Use Decimal
    inputs to retain input digits beyond binary64 precision.
backend : {"auto", "native", "symjit", "expression", "symbolica"}, optional
    Evaluation backend; default "auto" selects a supported backend for the
    requested precision. "symjit" supports binary64 only; "symbolica" is
    an alias for "expression"."#
);
scalar!(
    c0,
    "C0",
    C0,
    [
        p1_squared,
        p2_squared,
        p3_squared,
        mass_0_squared,
        mass_1_squared,
        mass_2_squared
    ],
    r#"Evaluate C0(p1_squared, p2_squared, p3_squared, mass_0_squared, mass_1_squared, mass_2_squared, mu_squared).

Results are ordered (finite, 1/eps, 1/eps**2). Invariants and masses are
squared quantities; omitted mu_squared is exactly one. ``prec`` is the
positive decimal significant-digit count (default 16). Decimal inputs or
higher precision produce DecimalComplex values. ``backend`` selects "auto",
"native", "symjit" (binary64 only), "expression", or "symbolica".
``rebuild=True`` refreshes the evaluator workspace.

Examples
--------
>>> from symbolica import S, E
>>> from symbolica.community import hepkit as hep
>>> from symbolica.community.hepkit import oneloop
>>> finite, pole, double_pole = oneloop.c0(-1.0, -2.0, -3.0, 1.0, 1.0, 1.0, 1.0)

Parameters
----------
p1_squared, p2_squared, p3_squared : Number
    Three external momentum invariants, in C0 argument order.
mass_0_squared, mass_1_squared, mass_2_squared : Number
    Squared masses of the three propagators, in C0 argument order.
mu_squared : Number or None, optional
    Squared renormalization scale; None uses exactly one.
rebuild : bool, optional
    Refresh the evaluator workspace before evaluation; default False.
prec : int, optional
    Positive number of decimal significant digits; default 16. Use Decimal
    inputs to retain input digits beyond binary64 precision.
backend : {"auto", "native", "symjit", "expression", "symbolica"}, optional
    Evaluation backend; default "auto" selects a supported backend for the
    requested precision. "symjit" supports binary64 only; "symbolica" is
    an alias for "expression"."#
);
scalar!(
    d0,
    "D0",
    D0,
    [
        p1_squared,
        p2_squared,
        p3_squared,
        p4_squared,
        s12,
        s23,
        mass_0_squared,
        mass_1_squared,
        mass_2_squared,
        mass_3_squared
    ],
    r#"Evaluate D0 with four external squared momenta, s12, s23, four squared masses, and the squared scale.

Results are ordered (finite, 1/eps, 1/eps**2). Invariants and masses are
squared quantities; omitted mu_squared is exactly one. ``prec`` is the
positive decimal significant-digit count (default 16). Decimal inputs or
higher precision produce DecimalComplex values. ``backend`` selects "auto",
"native", "symjit" (binary64 only), "expression", or "symbolica".
``rebuild=True`` refreshes the evaluator workspace.

Examples
--------
>>> from symbolica import S, E
>>> from symbolica.community import hepkit as hep
>>> from symbolica.community.hepkit import oneloop
>>> finite, pole, double_pole = oneloop.d0(-1.0, -1.0, -1.0, -1.0, -3.0, -4.0, 1.0, 1.0, 1.0, 1.0, 1.0)

Parameters
----------
p1_squared, p2_squared, p3_squared, p4_squared : Number
    Four external momenta squared, in D0 argument order.
s12 : Number
    Channel invariant (p1 + p2)**2.
s23 : Number
    Channel invariant (p2 + p3)**2.
mass_0_squared, mass_1_squared, mass_2_squared, mass_3_squared : Number
    Squared masses of the four propagators, in D0 argument order.
mu_squared : Number or None, optional
    Squared renormalization scale; None uses exactly one.
rebuild : bool, optional
    Refresh the evaluator workspace before evaluation; default False.
prec : int, optional
    Positive number of decimal significant digits; default 16. Use Decimal
    inputs to retain input digits beyond binary64 precision.
backend : {"auto", "native", "symjit", "expression", "symbolica"}, optional
    Evaluation backend; default "auto" selects a supported backend for the
    requested precision. "symjit" supports binary64 only; "symbolica" is
    an alias for "expression"."#
);

mod inspection;

fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    inspection::register(module)?;
    // Registering Python types/functions does not initialize Symbolica symbols.
    module.add_class::<Evaluator>()?;
    module.add_class::<precision::DecimalComplex>()?;
    module.add_function(wrap_pyfunction!(a0, module)?)?;
    module.add_function(wrap_pyfunction!(b0, module)?)?;
    module.add_function(wrap_pyfunction!(db0, module)?)?;
    module.add_function(wrap_pyfunction!(c0, module)?)?;
    module.add_function(wrap_pyfunction!(d0, module)?)?;
    module.add_function(wrap_pyfunction!(is_initialized, module)?)?;
    for (lower, upper) in [
        ("a0", "A0"),
        ("b0", "B0"),
        ("db0", "dB0"),
        ("c0", "C0"),
        ("d0", "D0"),
    ] {
        module.add(lower, module.getattr(upper)?)?;
    }
    module.add("COEFFICIENT_ORDER", (0, -1, -2))?;
    module.add("DEFAULT_BACKEND", backend_name(oneloop::DEFAULT_BACKEND))?;
    module.add(
        "SYMBOLICA_REVISION",
        "a19c760dd567c239f30d87e4e924ca2f8b8457ab",
    )?;
    module.add("EXPRESSION_INTEROP", cfg!(feature = "community"))?;
    #[cfg(feature = "community")]
    community::register(module)?;
    Ok(())
}

#[cfg(all(feature = "extension-module", not(feature = "community")))]
#[pymodule]
fn oneloop_native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    register(module)
}

#[cfg(feature = "community")]
mod community {
    use super::*;
    use symbolica::{
        api::python::{
            Citation, ConvertibleToExpression, PythonExpression, PythonExpressionEvaluator,
            SymbolicaCommunityModule,
        },
        atom::Atom,
    };

    pub struct CommunityModule;

    impl SymbolicaCommunityModule for CommunityModule {
        fn get_citations() -> Vec<Citation> {
            if !oneloop::was_used() {
                return Vec::new();
            }
            vec![
Citation {
                id: "https://github.com/alphal00p/oneloopmaster".into(),
                url: "https://github.com/alphal00p/oneloopmaster".into(),
                reference: "OneLoopMaster contributors. OneLoopMaster (2026).".into(),
                bibtex: r#"@software{oneloopmaster,
  author = {{OneLoopMaster contributors}},
  title = {{OneLoopMaster}},
  year = {2026},
  url = {https://github.com/alphal00p/oneloopmaster}
}"#.into(),
                reasons: vec!["Rust and Symbolica implementation of the scalar master integrals.".into()],
                description: String::new(),
                relevance: None,
            },
Citation {
                id: "arXiv:1007.4716".into(),
                url: "https://arxiv.org/abs/1007.4716".into(),
                reference: "A. van Hameren. OneLOop: for the evaluation of one-loop scalar functions. Computer Physics Communications 182 (2011) 2427–2438.".into(),
                bibtex: r#"@article{vanHameren2011OneLOop,
  author = {van Hameren, Andreas},
  title = {{OneLOop: for the evaluation of one-loop scalar functions}},
  year = {2011},
  url = {https://arxiv.org/abs/1007.4716},
  doi = {10.1016/j.cpc.2011.06.011},
  eprint = {1007.4716},
  archivePrefix = {arXiv}
}"#.into(),
                reasons: vec!["Original OneLOop algorithms; requested in the OneLoopMaster README.".into()],
                description: String::new(),
                relevance: None,
            },
Citation {
                id: "arXiv:0903.4665".into(),
                url: "https://arxiv.org/abs/0903.4665".into(),
                reference: "A. van Hameren, C. G. Papadopoulos and R. Pittau. Automated one-loop calculations: a proof of concept. JHEP 09 (2009) 106.".into(),
                bibtex: r#"@article{vanHameren2009Automated,
  author = {van Hameren, Andreas and Papadopoulos, Costas G. and Pittau, Roberto},
  title = {{Automated one-loop calculations: a proof of concept}},
  year = {2009},
  url = {https://arxiv.org/abs/0903.4665},
  doi = {10.1088/1126-6708/2009/09/106},
  eprint = {0903.4665},
  archivePrefix = {arXiv}
}"#.into(),
                reasons: vec!["Original OneLOop work; requested in the OneLoopMaster README.".into()],
                description: String::new(),
                relevance: None,
            }
            ]
        }

        fn get_name() -> String {
            "oneloop".to_owned()
        }
        fn register_module(module: &Bound<'_, PyModule>) -> PyResult<()> {
            super::register(module)
        }
        fn initialize(_py: Python<'_>) -> PyResult<()> {
            Ok(())
        }
    }

    /// Return the three Laurent coefficients of a complete primitive master call.
    ///
    /// The input includes physical arguments and squared scale, without a Laurent
    /// tag. The returned symbolic expressions carry tags 0, -1 and -2 for finite,
    /// simple-pole and double-pole coefficients and retain native evaluation hooks.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hepkit as hep
    /// >>> from symbolica.community.hepkit import oneloop
    /// >>> s = S("s")
    /// >>> coefficients = oneloop.master_coefficients(oneloop.B0(s, 0, 0, 1))
    /// >>> assert len(coefficients) == 3
    ///
    /// Parameters
    /// ----------
    /// master : Expression
    ///     Complete untagged A0/B0/dB0/C0/D0 call, including squared scale.
    #[pyfunction]
    fn master_coefficients(master: PythonExpression) -> PyResult<Vec<PythonExpression>> {
        let (family, input) =
            oneloop::master_arguments(master.expr.as_view()).map_err(PyValueError::new_err)?;
        let master = match family {
            ScalarIntegral::A0 => oneloop::A0(),
            ScalarIntegral::B0 => oneloop::B0(),
            ScalarIntegral::DB0 => oneloop::dB0(),
            ScalarIntegral::C0 => oneloop::C0(),
            ScalarIntegral::D0 => oneloop::D0(),
        };
        Ok([0, -1, -2]
            .into_iter()
            .map(|tag| {
                let arguments = std::iter::once(Atom::num(tag))
                    .chain(input.iter().cloned())
                    .collect::<Vec<_>>();
                master.call(arguments.as_slice()).into()
            })
            .collect())
    }

    const MASTER_NAMES: [&str; 5] = ["A0", "B0", "dB0", "C0", "D0"];

    /// Resolve and cache the core's primitive master Symbol on first access.
    ///
    /// Symbolica's community registration phase must not initialize symbols.
    /// Module attribute lookup defers this lightweight registration to use;
    /// it never builds formulas, numerical evaluators, or warmup caches.
    #[pyfunction(pass_module, name = "__getattr__")]
    fn primitive_symbol(module: &Bound<'_, PyModule>, name: &str) -> PyResult<Py<PyAny>> {
        if name == "__all__" {
            // Hosts can attach reduction types/functions after registration.
            // Discover their exports when star import actually asks for them.
            let names = module_dir(module)?
                .into_iter()
                .filter(|name| !name.starts_with('_'))
                .collect::<Vec<_>>();
            return Ok(names.into_pyobject(module.py())?.into_any().unbind());
        }
        let symbol = match name {
            "A0" => oneloop::A0(),
            "B0" => oneloop::B0(),
            "dB0" => oneloop::dB0(),
            "C0" => oneloop::C0(),
            "D0" => oneloop::D0(),
            _ => {
                return Err(pyo3::exceptions::PyAttributeError::new_err(format!(
                    "module '{}' has no attribute '{name}'",
                    module.name()?
                )));
            }
        };
        let expression = Py::new(module.py(), PythonExpression::from(Atom::var(symbol)))?;
        module.add(name, expression.clone_ref(module.py()))?;
        Ok(expression.into_any())
    }

    #[pyfunction(pass_module, name = "__dir__")]
    fn module_dir(module: &Bound<'_, PyModule>) -> PyResult<Vec<String>> {
        let mut names = module.dict().keys().extract::<Vec<String>>()?;
        names.extend(MASTER_NAMES.map(str::to_owned));
        names.sort();
        names.dedup();
        Ok(names)
    }

    /// Compile symbolic combinations of master coefficients for repeated evaluation.
    ///
    /// Returns a Symbolica evaluator retaining the primitive master definitions.
    /// Use its complex evaluation method and pass parameter values in the supplied
    /// order. This operation builds an evaluator; it does not evaluate a point.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hepkit as hep
    /// >>> from symbolica.community.hepkit import oneloop
    /// >>> m2 = S("m2")
    /// >>> coefficients = oneloop.master_coefficients(oneloop.A0(m2, 1))
    /// >>> compiled = oneloop.compile_native(coefficients, [m2])
    /// >>> values = compiled.evaluate_complex([1+0j])
    ///
    /// Parameters
    /// ----------
    /// expressions : sequence[Expression | int | float | complex]
    ///     Outputs to compile, usually master or reduction coefficients.
    /// parameters : sequence[Expression]
    ///     Ordered independent symbols receiving numerical values.
    #[pyfunction]
    fn compile_native(
        expressions: Vec<ConvertibleToExpression>,
        parameters: Vec<PythonExpression>,
    ) -> PyResult<PythonExpressionEvaluator> {
        oneloop::record_usage();
        let expressions = expressions
            .into_iter()
            .map(|e| e.to_expression().expr)
            .collect::<Vec<_>>();
        let parameters = parameters.into_iter().map(|p| p.expr).collect::<Vec<_>>();
        let exact = oneloop::OneLoopExpressions::new()
            .evaluator(&expressions, &parameters)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(PythonExpressionEvaluator {
            rational_constants: exact.get_constants().to_vec(),
            eval_complex: exact.map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64())),
            eval_real: None,
            jit_real: None,
            jit_complex: None,
            eval_double_float: None,
            eval_double_float_complex: None,
            eval_arb_prec: None,
            eval_arb_prec_complex: None,
            jit_compile: false,
            jit_settings: oneloop::jit_settings(),
        })
    }

    pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
        module.add_function(wrap_pyfunction!(master_coefficients, module)?)?;
        module.add_function(wrap_pyfunction!(compile_native, module)?)?;
        // super::register has already retained the numeric wrappers under
        // lowercase names. Uppercase names now resolve to the core Symbols.
        for name in MASTER_NAMES {
            module.delattr(name)?;
        }
        module.add_function(wrap_pyfunction!(primitive_symbol, module)?)?;
        module.add_function(wrap_pyfunction!(module_dir, module)?)?;
        Ok(())
    }
}

#[cfg(feature = "community")]
pub use community::CommunityModule;

/// Attach OneLOop to an existing community HEP module in the same extension.
/// The host initializes Symbolica as usual. OneLOop evaluator construction
/// remains lazy; the community initialize hook performs no numerical work.
#[cfg(feature = "community")]
pub fn register_hep_module(hep: &Bound<'_, PyModule>) -> PyResult<()> {
    let module = PyModule::new(hep.py(), "symbolica.community.hepkit_oneloop_native")?;
    register(&module)?;
    hep.add("_oneloop_native", &module)?;
    hep.py()
        .import("sys")?
        .getattr("modules")?
        .set_item("symbolica.community.hepkit_oneloop_native", &module)?;
    Ok(())
}
