//! Portable evaluator source caches and row-major batched numerical evaluation.
use crate::{ScalarIntegral, masters};
use std::{
    collections::HashSet,
    sync::{Arc, Mutex, OnceLock},
};
use symbolica::{
    domains::{float::Complex, rational::Rational},
    evaluate::{
        ExportedInstructions, ExpressionEvaluator, Instruction, JITCompilationSettings,
        JITCompiledEvaluator,
    },
};

type C = Complex<f64>;

pub(crate) fn portable_environment() -> Result<(), String> {
    // Config::default accepts architecture-pinned compiler types through these
    // overrides. Require the ordinary automatic backend both when constructing
    // a portable source cache and when recompiling it on the receiving host.
    if std::env::var_os("SYMJIT_TOML").is_some() {
        return Err(
            "unset SYMJIT_TOML when building or loading portable OneLOop evaluators".into(),
        );
    }
    if std::fs::exists("symjit.toml")
        .map_err(|error| format!("cannot check the default symjit.toml override: {error}"))?
    {
        return Err(
            "remove the working-directory symjit.toml override when building or loading portable OneLOop evaluators".into(),
        );
    }
    Ok(())
}

// The payload is Symbolica ExpressionEvaluator IR, not serialized SymJIT code.
// Keep its schema and compilation policy versioned; loading always recompiles
// with the backend selected by the consuming application's Cargo.lock.
const FORMAT: &[u8] = b"oneloop-evaluator-v3:a19c760dd567c239f30d87e4e924ca2f8b8457ab:symbolica-ir:strict-o2-complex1-simd0\0";
// These source-IR formats have identical payloads. The old backend tag recorded
// the producer, but never made its native code part of the cache.
const LEGACY_FORMATS: [&[u8]; 2] = [
    b"oneloop-evaluator-v2:a19c760dd567c239f30d87e4e924ca2f8b8457ab:symjit-2.26.0:strict-o2-complex1-simd0\0",
    b"oneloop-evaluator-v2:a19c760dd567c239f30d87e4e924ca2f8b8457ab:symjit-2.26.4:strict-o2-complex1-simd0\0",
];

fn cache_payload(data: &[u8]) -> Result<(&[u8], usize, usize), String> {
    let payload = data
        .strip_prefix(FORMAT)
        .or_else(|| {
            LEGACY_FORMATS
                .iter()
                .find_map(|prefix| data.strip_prefix(*prefix))
        })
        .ok_or("incompatible evaluator cache")?;
    let dimensions = payload.get(..16).ok_or("truncated evaluator cache")?;
    let dimension = |offset| {
        usize::try_from(u64::from_le_bytes(
            dimensions[offset..offset + 8].try_into().unwrap(),
        ))
        .map_err(|_| "evaluator dimension does not fit this platform".to_owned())
    };
    Ok((&payload[16..], dimension(0)?, dimension(8)?))
}

/// SymJIT O2 settings for scalar code and row-major batches.
/// SIMD remains disabled: SymJIT 2.26.0 mispacks scalar complex callbacks and
/// does not propagate nested conditional fallback. Ordinary translation and
/// packed complex arithmetic retain the validated production policy; fastmath
/// remains disabled. The earlier direct-power and strict-FMA defects are fixed.
pub fn jit_settings() -> JITCompilationSettings {
    JITCompilationSettings::new()
        .optimization_level(2)
        .direct_translation(false)
        .with_option("fastmath", "false")
        .with_option("fast_complex", "true")
        .with_option("use_threads", "false")
        .with_option("use_simd", "false")
        .with_option("simd_branch", "false")
        .with_option("enable_simd512", "false")
}

fn max_call_arity<T>(program: &ExportedInstructions<T>) -> usize {
    let mut pending = vec![program];
    let mut visited = HashSet::new();
    let mut maximum = 0;
    while let Some(program) = pending.pop() {
        for instruction in &program.instructions {
            if let Instruction::Fun(_, function, _) = instruction {
                maximum = maximum.max(function.2.len());
            }
        }
        for body in &program.sub_evaluators {
            if visited.insert(Arc::as_ptr(&body.instructions)) {
                pending.push(body.instructions.as_ref());
            }
        }
    }
    maximum
}

fn checked_argument_slots(complex_arguments: usize) -> Result<u32, String> {
    // SymJIT lowers every complex call argument to two real __Arg slots,
    // including arguments captured by shared evaluator bodies. The top-level
    // input count does not bound calls inside those bodies.
    let slots = complex_arguments
        .checked_mul(2)
        .and_then(|slots| u32::try_from(slots).ok())
        .ok_or("JIT function argument slot count overflow")?;
    let mut config = symjit::Config::default();
    config.set_num_args(slots);
    let accepted = config.num_args();
    if accepted < slots {
        return Err(format!(
            "JIT function requires {slots} real argument slots, but SymJIT accepts only {accepted}"
        ));
    }
    Ok(accepted)
}

fn jit_settings_for(source: &ExpressionEvaluator<C>) -> Result<JITCompilationSettings, String> {
    let slots = checked_argument_slots(max_call_arity(&source.export_instructions()))?;
    // Symbolica applies these options to the main evaluator and all shared
    // bodies. Keep every other part of the portable compilation policy intact.
    Ok(jit_settings().with_option("num_args", slots.to_string()))
}

/// A reusable evaluator for arbitrary combinations of native expressions.
/// Inputs and outputs of a batch are row-major, one complete point per row.
#[derive(Clone)]
pub struct JitEvaluator {
    inner: JITCompiledEvaluator<C>,
    source: Arc<ExpressionEvaluator<C>>,
    inputs: usize,
    outputs: usize,
}

impl JitEvaluator {
    pub(crate) fn compile(
        exact: &ExpressionEvaluator<Complex<Rational>>,
        inputs: usize,
        outputs: usize,
    ) -> Result<Self, String> {
        portable_environment()?;
        let source = exact
            .clone()
            .map_coeff(&|c| C::new(c.re.to_f64(), c.im.to_f64()));
        if source.get_input_len() != inputs || source.get_output_len() != outputs {
            return Err(
                "compiled evaluator dimensions disagree with requested expression shape".into(),
            );
        }
        let inner = source.jit_compile(jit_settings_for(&source)?)?;
        Ok(Self {
            inner,
            source: Arc::new(source),
            inputs,
            outputs,
        })
    }

    /// Number of scalar complex arguments per point.
    pub fn input_count(&self) -> usize {
        self.inputs
    }

    /// Number of complex outputs per point.
    pub fn output_count(&self) -> usize {
        self.outputs
    }

    /// Evaluate one point, checking the exact input and output dimensions.
    pub fn evaluate(&mut self, args: &[C], output: &mut [C]) -> Result<(), String> {
        crate::record_usage();
        self.validate(args.len(), output.len(), 1)?;
        self.inner.evaluate(args, output);
        Ok(())
    }

    /// Evaluate `rows` points from a flat row-major input array. The output is
    /// also row-major. Empty batches are valid and do not enter the backend.
    pub fn evaluate_batch(
        &mut self,
        args: &[C],
        output: &mut [C],
        rows: usize,
    ) -> Result<(), String> {
        crate::record_usage();
        self.validate(args.len(), output.len(), rows)?;
        if rows != 0 {
            self.inner.batch_evaluate(args, output, rows);
        }
        Ok(())
    }

    fn validate(&self, args: usize, output: usize, rows: usize) -> Result<(), String> {
        let expected_in = rows.checked_mul(self.inputs).ok_or("input size overflow")?;
        let expected_out = rows
            .checked_mul(self.outputs)
            .ok_or("output size overflow")?;
        if args != expected_in || output != expected_out {
            return Err(format!(
                "expected {expected_in} inputs and {expected_out} outputs for {rows} rows; received {args} and {output}"
            ));
        }
        Ok(())
    }

    /// Serialize Symbolica evaluator IR and every nested function definition.
    /// No SymJIT application or native machine code is serialized.
    /// Treat this as trusted executable content, not an untrusted data format.
    pub fn to_bytes(&self) -> Result<Vec<u8>, String> {
        let mut data = FORMAT.to_vec();
        data.extend_from_slice(&(self.inputs as u64).to_le_bytes());
        data.extend_from_slice(&(self.outputs as u64).to_le_bytes());
        data.extend(
            bincode::encode_to_vec(self.source.as_ref(), bincode::config::standard())
                .map_err(|error| error.to_string())?,
        );
        Ok(data)
    }

    /// Load a schema-matched trusted source cache. Native machine code is regenerated
    /// for the host with the same strict settings at every nesting level; no
    /// symbolic expression construction is required.
    /// Register any custom numerical Symbol callbacks before loading. The
    /// `SYMJIT_TOML` must be unset and the working directory must not contain
    /// `symjit.toml`, since these can override the portable compiler settings.
    pub fn from_bytes(data: &[u8]) -> Result<Self, String> {
        crate::initialization::ensure_symbolica_state();
        Self::from_bytes_raw(data)
    }

    fn from_bytes_raw(data: &[u8]) -> Result<Self, String> {
        portable_environment()?;
        let (payload, inputs, outputs) = cache_payload(data)?;
        let (source, consumed): (ExpressionEvaluator<C>, usize) =
            bincode::decode_from_slice(payload, bincode::config::standard())
                .map_err(|error| error.to_string())?;
        if consumed != payload.len() {
            return Err("trailing bytes in evaluator cache".into());
        }
        // Dimensions in the outer header cannot authorize shorter slices than
        // the executable consumes. SymJIT enters raw generated machine code.
        if source.get_input_len() != inputs || source.get_output_len() != outputs {
            return Err("evaluator cache dimensions disagree with executable".into());
        }
        let inner = source.jit_compile(jit_settings_for(&source)?)?;
        Ok(Self {
            inputs,
            outputs,
            inner,
            source: Arc::new(source),
        })
    }
}

/// Three-coefficient evaluator for one scalar family. The last argument of
/// every point is the positive real squared renormalization scale `mu_squared`.
#[derive(Clone)]
pub struct ScalarEvaluator {
    family: ScalarIntegral,
    evaluator: JitEvaluator,
}

impl ScalarEvaluator {
    /// Load the embedded backend on first use, then clone its compiled code.
    /// The returned instance has its own workspace.
    pub fn prebuilt(family: ScalarIntegral) -> Result<Self, String> {
        let _ = embedded(family)?;
        Self::cached(family)
    }

    /// Prepare this family on first use, then clone its shared backend.
    /// Without `prebuilt`, build this family from native expressions.
    pub fn cached(family: ScalarIntegral) -> Result<Self, String> {
        cache(family)?
            .lock()
            .map_err(|_| "OneLOop JIT cache poisoned".into())
            .map(|backend| backend.evaluator.clone())
    }

    /// Rebuild O2 code from the current transparent native expressions.
    /// As with manual expression construction, use an adequately sized stack.
    pub fn rebuild(family: ScalarIntegral) -> Result<Self, String> {
        crate::initialization::ensure_symbolica_state();
        Self::rebuild_raw(family)
    }

    fn rebuild_raw(family: ScalarIntegral) -> Result<Self, String> {
        portable_environment()?;
        Ok(Self {
            family,
            evaluator: masters::jit_evaluator(family)?,
        })
    }

    /// Integral family represented by this evaluator.
    pub fn family(&self) -> ScalarIntegral {
        self.family
    }

    /// Evaluate a single point into `[finite, simple pole, double pole]`.
    pub fn evaluate(&mut self, args: &[C], output: &mut [C]) -> Result<(), String> {
        crate::record_usage();
        let mut relabelled = Vec::new();
        let args = if self.family == ScalarIntegral::D0 {
            crate::labelling::regular_point(args, &mut relabelled)
        } else {
            args
        };
        self.evaluator.evaluate(args, output)
    }

    /// Evaluate row-major points, writing three consecutive coefficients per row.
    pub fn evaluate_batch(
        &mut self,
        args: &[C],
        output: &mut [C],
        rows: usize,
    ) -> Result<(), String> {
        crate::record_usage();
        let mut relabelled = Vec::new();
        let args = if self.family == ScalarIntegral::D0 {
            crate::labelling::regular_rows(args, &mut relabelled)
        } else {
            args
        };
        self.evaluator.evaluate_batch(args, output, rows)
    }

    /// Serialize this family cache, including its family and dimension checks.
    pub fn to_bytes(&self) -> Result<Vec<u8>, String> {
        let mut bytes = vec![self.family as u8];
        bytes.extend(self.evaluator.to_bytes()?);
        Ok(bytes)
    }

    /// Load trusted bytes for a specific family; rejects mismatched metadata.
    pub fn from_bytes(family: ScalarIntegral, bytes: &[u8]) -> Result<Self, String> {
        crate::initialization::ensure_symbolica_state();
        Self::from_bytes_raw(family, bytes)
    }

    fn from_bytes_raw(family: ScalarIntegral, bytes: &[u8]) -> Result<Self, String> {
        if bytes.first().copied() != Some(family as u8) {
            return Err("wrong scalar family in evaluator cache".into());
        }
        let (_, inputs, outputs) = cache_payload(&bytes[1..])?;
        if inputs != family.arity() || outputs != 3 {
            return Err("wrong dimensions in scalar evaluator cache".into());
        }
        let evaluator = JitEvaluator::from_bytes_raw(&bytes[1..])?;
        if evaluator.input_count() != family.arity() || evaluator.output_count() != 3 {
            return Err("wrong dimensions in scalar evaluator cache".into());
        }
        Ok(Self { family, evaluator })
    }
}

#[cfg(feature = "prebuilt")]
fn embedded(family: ScalarIntegral) -> Result<&'static [u8], String> {
    Ok(match family {
        ScalarIntegral::A0 => include_bytes!("../assets/evaluators/a0.bin"),
        ScalarIntegral::B0 => include_bytes!("../assets/evaluators/b0.bin"),
        ScalarIntegral::DB0 => include_bytes!("../assets/evaluators/db0.bin"),
        ScalarIntegral::C0 => include_bytes!("../assets/evaluators/c0.bin"),
        ScalarIntegral::D0 => include_bytes!("../assets/evaluators/d0.bin"),
    })
}

#[cfg(not(feature = "prebuilt"))]
fn embedded(_family: ScalarIntegral) -> Result<&'static [u8], String> {
    Err("embedded evaluators are disabled; use ScalarEvaluator::rebuild".into())
}

// A master exposes scalar Laurent tags, while its backend computes all three
// coefficients together. Keep four FIFO groups for tag-major SIMD callbacks:
// all lanes request one tag before requesting the next. Each tag is consumed
// once per group, including when several lanes have identical arguments.
// This avoids three backend calls per point without memoizing repeated points
// in throughput surveys. Bitwise keys preserve signed branch-axis zeros.
struct CoefficientGroup {
    arguments: [C; 11],
    values: [C; 3],
    used_tags: u8,
}

impl CoefficientGroup {
    fn new() -> Self {
        Self {
            arguments: [C::new(0., 0.); 11],
            values: [C::new(0., 0.); 3],
            used_tags: 0,
        }
    }

    fn reusable(&self, arguments: &[C], tag: usize) -> bool {
        self.used_tags != 0
            && self.used_tags & (1 << tag) == 0
            && arguments
                .iter()
                .zip(&self.arguments)
                .all(|(a, b)| a.re.to_bits() == b.re.to_bits() && a.im.to_bits() == b.im.to_bits())
    }

    fn invalidate(&mut self) {
        self.used_tags = 0;
    }
}

pub(crate) struct CoefficientGroups {
    groups: [CoefficientGroup; 4],
    oldest: usize,
}

impl CoefficientGroups {
    pub(crate) fn new() -> Self {
        Self {
            groups: core::array::from_fn(|_| CoefficientGroup::new()),
            oldest: 0,
        }
    }

    pub(crate) fn invalidate(&mut self) {
        for group in &mut self.groups {
            group.invalidate();
        }
        self.oldest = 0;
    }

    pub(crate) fn evaluate(
        &mut self,
        arguments: &[C],
        tag: usize,
        evaluate: impl FnOnce(&mut [C; 3]),
    ) -> C {
        crate::record_usage();
        for offset in 0..self.groups.len() {
            let index = (self.oldest + offset) % self.groups.len();
            let group = &mut self.groups[index];
            if group.reusable(arguments, tag) {
                group.used_tags |= 1 << tag;
                return group.values[tag];
            }
        }

        let index = self.oldest;
        self.oldest = (self.oldest + 1) % self.groups.len();
        let group = &mut self.groups[index];
        // Invalidate the evicted entry before calling the backend. If it
        // panics, this entry stays invalid and the shared mutex is poisoned;
        // no partially updated output is published. Other entries are intact.
        group.invalidate();
        let mut output = [C::new(0., 0.); 3];
        evaluate(&mut output);
        group.arguments[..arguments.len()].copy_from_slice(arguments);
        group.values = output;
        group.used_tags = 1 << tag;
        group.values[tag]
    }
}

struct PreparedBackend {
    evaluator: ScalarEvaluator,
}

static CACHED: [OnceLock<Result<Mutex<PreparedBackend>, String>>; 5] =
    [const { OnceLock::new() }; 5];

pub(crate) fn initialize_all() -> Result<(), String> {
    for family in [
        ScalarIntegral::A0,
        ScalarIntegral::B0,
        ScalarIntegral::DB0,
        ScalarIntegral::C0,
        ScalarIntegral::D0,
    ] {
        let _ = cache(family)?;
    }
    Ok(())
}

fn cache(family: ScalarIntegral) -> Result<&'static Mutex<PreparedBackend>, String> {
    crate::initialization::ensure_symbolica_state();
    portable_environment()?;
    CACHED[family as usize]
        .get_or_init(|| {
            #[cfg(feature = "prebuilt")]
            let evaluator = ScalarEvaluator::from_bytes_raw(family, embedded(family)?);
            #[cfg(not(feature = "prebuilt"))]
            let evaluator = ScalarEvaluator::rebuild_raw(family);
            evaluator
                .map(|evaluator| Mutex::new(PreparedBackend { evaluator }))
                .map_err(|error| {
                    format!("could not initialize {} evaluator: {error}", family.name())
                })
        })
        .as_ref()
        .map_err(Clone::clone)
}

/// Evaluate flat row-major points, preparing this family on first use.
pub fn evaluate_batch(
    family: ScalarIntegral,
    args: &[C],
    output: &mut [C],
    rows: usize,
) -> Result<(), String> {
    crate::record_usage();
    let mut backend = cache(family)?
        .lock()
        .map_err(|_| "OneLOop JIT cache poisoned")?;
    backend.evaluator.evaluate_batch(args, output, rows)
}

/// Rebuild and replace a family's shared manual SymJIT evaluator.
/// Direct Rust Symbol hooks and other existing manual evaluators are unaffected.
pub fn rebuild_cached_evaluator(family: ScalarIntegral) -> Result<(), String> {
    let evaluator = ScalarEvaluator::rebuild(family)?;
    let mut backend = cache(family)?
        .lock()
        .map_err(|_| "OneLOop JIT cache poisoned")?;
    backend.evaluator = evaluator;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jit_argument_slots_reject_overflow_and_backend_clamping() {
        assert!(checked_argument_slots(usize::MAX).is_err());
        assert!(checked_argument_slots(u32::MAX as usize).is_err());
        let mut backend = symjit::Config::default();
        backend.set_num_args(u32::MAX);
        let too_many = backend.num_args() as usize / 2 + 1;
        assert!(
            checked_argument_slots(too_many)
                .unwrap_err()
                .contains("SymJIT accepts only")
        );
        assert!(checked_argument_slots(35).unwrap() >= 70);
    }

    #[test]
    fn jit_nested_wide_calls_survive_source_cache_and_mixed_batches() {
        use symbolica::{
            evaluate::{FunctionMap, FunctionRegistrationOptions, InliningPolicy},
            prelude::*,
        };

        let x = symbol!("jit_wide_outer_x");
        let leaf = symbol!("jit_wide_leaf");
        let first = symbol!("jit_wide_first");
        let second = symbol!("jit_wide_second");
        let parameters: Vec<_> = (0..35)
            .map(|index| symbol!(format!("jit_wide_argument_{index}")))
            .collect();
        let body = parameters.iter().fold(Atom::num(0), |sum, &argument| {
            sum + argument.to_atom() * argument
        });
        let options = FunctionRegistrationOptions::new().inlining(InliningPolicy::Never);
        let mut map = FunctionMap::new();
        map.add_function_with_options(leaf, parameters, body, options.clone())
            .unwrap();
        // Repeated argument slots introduce no captured scalar constants. Thus
        // the public expression calls only narrow helpers; only their nested
        // calls require the enlarged argument storage.
        let arguments = vec![x.to_atom(); 35];
        let call = leaf.call(arguments.as_slice());
        map.add_function_with_options(first, vec![x], call.clone(), options.clone())
            .unwrap();
        map.add_function_with_options(second, vec![x], call + x, options)
            .unwrap();
        let expression = first.call((x,)) + second.call((2 * x.to_atom(),));
        let exact = expression
            .evaluator(&[x.to_atom()])
            .function_map(map)
            .direct_translation(true)
            .horner_iterations(0)
            .build()
            .unwrap();
        assert_eq!(exact.get_input_len(), 1);
        let exported = exact.export_instructions();
        assert!(
            exported
                .sub_evaluators
                .iter()
                .all(|body| body.input_count < 35)
        );
        assert!(max_call_arity(&exported) >= 35);
        let mut interpreted = exact
            .clone()
            .map_coeff(&|value| C::new(value.re.to_f64(), value.im.to_f64()));
        let fresh = JitEvaluator::compile(&exact, 1, 1).unwrap();
        let bytes = fresh.to_bytes().unwrap();
        let restored = JitEvaluator::from_bytes(&bytes).unwrap();
        let arguments = [
            C::new(0., 0.),
            C::new(1., 0.),
            C::new(-2., 3.),
            C::new(0., -0.125),
            C::new(0.03125, 8.),
            C::new(-0.25, -0.5),
            C::new(17., 0.),
        ];
        let expected: Vec<_> = arguments
            .iter()
            .map(|argument| interpreted.evaluate_single(std::slice::from_ref(argument)))
            .collect();
        for mut evaluator in [fresh, restored] {
            evaluator.evaluate_batch(&[], &mut [], 0).unwrap();
            for batch in [1, 3, 8] {
                let mut actual = [C::new(f64::NAN, f64::NAN); 7];
                for (input, output) in arguments.chunks(batch).zip(actual.chunks_mut(batch)) {
                    evaluator
                        .evaluate_batch(input, output, input.len())
                        .unwrap();
                }
                for (actual, expected) in actual.iter().zip(&expected) {
                    let scale = 1. + expected.re.abs() + expected.im.abs();
                    assert!((actual.re - expected.re).abs() <= 1e-12 * scale);
                    assert!((actual.im - expected.im).abs() <= 1e-12 * scale);
                }
            }
        }
    }

    fn counted_evaluate(
        group: &mut CoefficientGroups,
        arguments: &[C],
        tag: usize,
        calls: &mut usize,
    ) -> C {
        group.evaluate(arguments, tag, |output| {
            *calls += 1;
            for (index, value) in output.iter_mut().enumerate() {
                *value = C::new(*calls as f64, index as f64);
            }
        })
    }

    #[test]
    fn coefficient_groups_execute_once_per_point_in_every_tag_order() {
        let orders = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];
        let args = [C::new(2., 0.), C::new(1., 0.)];
        for first in orders {
            for second in orders {
                let mut group = CoefficientGroups::new();
                let mut calls = 0;
                for (point, order) in [first, second].into_iter().enumerate() {
                    for tag in order {
                        assert_eq!(
                            counted_evaluate(&mut group, &args, tag, &mut calls),
                            C::new((point + 1) as f64, tag as f64),
                            "orders={first:?}/{second:?}, point={point}, tag={tag}"
                        );
                        assert_eq!(calls, point + 1);
                    }
                }
            }
        }
    }

    #[test]
    fn coefficient_groups_repeated_tags_always_start_another_execution() {
        let args = [C::new(2., 0.), C::new(1., 0.)];
        for tag in 0..3 {
            let mut group = CoefficientGroups::new();
            let mut calls = 0;
            for expected_calls in 1..=4 {
                assert_eq!(
                    counted_evaluate(&mut group, &args, tag, &mut calls),
                    C::new(expected_calls as f64, tag as f64)
                );
                assert_eq!(calls, expected_calls);
            }
            for other_tag in (0..3).filter(|other| *other != tag) {
                assert_eq!(
                    counted_evaluate(&mut group, &args, other_tag, &mut calls),
                    C::new(1., other_tag as f64)
                );
                assert_eq!(calls, 4);
            }
        }
    }

    #[test]
    fn coefficient_groups_execute_again_for_signed_zeros_and_changed_points() {
        let args = [C::new(0., 0.); 11];
        // Exercise both component keys at every argument position, including
        // the last slot used by the largest scalar family.
        for index in 0..args.len() {
            for imaginary in [false, true] {
                let mut group = CoefficientGroups::new();
                let mut calls = 0;
                counted_evaluate(&mut group, &args, 0, &mut calls);
                let mut signed_axis = args;
                if imaginary {
                    signed_axis[index].im = -0.;
                } else {
                    signed_axis[index].re = -0.;
                }
                assert_eq!(
                    counted_evaluate(&mut group, &signed_axis, 1, &mut calls),
                    C::new(2., 1.)
                );
                assert_eq!(calls, 2);
                assert_eq!(
                    counted_evaluate(&mut group, &args, 2, &mut calls),
                    C::new(1., 2.)
                );
                assert_eq!(calls, 2);
                let mut other_point = args;
                other_point[index].re = 4.;
                assert_eq!(
                    counted_evaluate(&mut group, &other_point, 0, &mut calls),
                    C::new(3., 0.)
                );
                assert_eq!(calls, 3);
            }
        }
    }

    #[test]
    fn coefficient_groups_invalidate_explicitly_and_before_backend_failure() {
        let args = [C::new(2., 0.), C::new(1., 0.)];
        let mut group = CoefficientGroups::new();
        let mut calls = 0;
        counted_evaluate(&mut group, &args, 0, &mut calls);
        group.invalidate();
        assert_eq!(
            counted_evaluate(&mut group, &args, 1, &mut calls),
            C::new(2., 1.)
        );
        assert_eq!(calls, 2);

        // Fill the ring, leaving the previous point in the next victim slot.
        for index in 1..4 {
            let mut lane = args;
            lane[0].re += index as f64;
            counted_evaluate(&mut group, &lane, 0, &mut calls);
        }
        assert_eq!(calls, 5);

        let mut other_point = args;
        other_point[1].re = 4.;
        let failure = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            group.evaluate(&other_point, 2, |output| {
                calls += 1;
                output[0] = C::new(123., 456.);
                panic!("simulated backend failure after partial output");
            });
        }));
        assert!(failure.is_err());
        assert_eq!(calls, 6);
        assert_eq!(group.groups[0].used_tags, 0);
        assert!(group.groups[1..].iter().all(|entry| entry.used_tags == 1));
        // Tag 2 of the old point was unused: reusing it here would conceal a
        // missing invalidation on the failed transition.
        assert_eq!(
            counted_evaluate(&mut group, &args, 2, &mut calls),
            C::new(7., 2.)
        );
        assert_eq!(calls, 7);
    }

    #[test]
    fn coefficient_groups_tag_major_lanes_reuse_oldest_matching_point() {
        let orders = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];
        for lanes in [1, 2, 4] {
            for identical in [false, true] {
                for first in orders {
                    for second in orders {
                        let mut group = CoefficientGroups::new();
                        let mut calls = 0;
                        for (batch, order) in [first, second].into_iter().enumerate() {
                            for tag in order {
                                for lane in 0..lanes {
                                    let args = [
                                        C::new(if identical { 2. } else { 2. + lane as f64 }, 0.),
                                        C::new(1., 0.),
                                    ];
                                    assert_eq!(
                                        counted_evaluate(&mut group, &args, tag, &mut calls),
                                        C::new((batch * lanes + lane + 1) as f64, tag as f64),
                                        "lanes={lanes}, identical={identical}, orders={first:?}/{second:?}, batch={batch}, lane={lane}, tag={tag}"
                                    );
                                }
                                assert_eq!(calls, (batch + 1) * lanes);
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn coefficient_groups_evict_fifo_without_reusing_consumed_tags() {
        let mut group = CoefficientGroups::new();
        let mut calls = 0;
        for point in 0..5 {
            counted_evaluate(&mut group, &[C::new(point as f64, 0.)], 0, &mut calls);
        }
        assert_eq!(calls, 5);
        // Point zero was evicted; its missing tag must execute again.
        assert_eq!(
            counted_evaluate(&mut group, &[C::new(0., 0.)], 1, &mut calls),
            C::new(6., 1.)
        );
        assert_eq!(calls, 6);
        // Point two is still buffered and its tag remains unused.
        assert_eq!(
            counted_evaluate(&mut group, &[C::new(2., 0.)], 1, &mut calls),
            C::new(3., 1.)
        );
        assert_eq!(calls, 6);
        assert_eq!(
            counted_evaluate(&mut group, &[C::new(2., 0.)], 1, &mut calls),
            C::new(7., 1.)
        );
        assert_eq!(calls, 7);
    }

    #[test]
    fn coefficient_groups_manual_invalidation_clears_all_lanes() {
        for identical in [false, true] {
            let mut group = CoefficientGroups::new();
            let mut calls = 0;
            let args = core::array::from_fn::<_, 4, _>(|lane| {
                [C::new(if identical { 2. } else { lane as f64 }, 0.)]
            });
            for lane in &args {
                counted_evaluate(&mut group, lane, 0, &mut calls);
            }
            group.invalidate();
            assert!(group.groups.iter().all(|entry| entry.used_tags == 0));
            for (lane, arguments) in args.iter().enumerate() {
                assert_eq!(
                    counted_evaluate(&mut group, arguments, 1, &mut calls),
                    C::new((5 + lane) as f64, 1.)
                );
            }
            assert_eq!(calls, 8);
        }
    }
}
