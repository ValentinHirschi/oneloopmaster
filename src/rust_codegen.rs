//! Generic Rust emission from exact Symbolica evaluator instructions.
//!
//! This is the emitter used by `examples/generate_native.rs`. It preserves
//! Symbolica's shared bodies and lazy branches; it performs no symbolic algebra.
//! Entry points may have any output count. Shared bodies must have one output.
//!
//! Emitted modules import `C`, `Context`, `NativeFloat` and `primitives` from
//! their parent module. Re-export [`runtime`] there to use the existing native
//! numeric implementation. The parent also supplies `ConstantSpec` for the
//! generated constant table. Build a `Context` once from that table, then call
//! each child module's `evaluate` with its declared number of input components.
//! Input length/precision validation belongs to the caller; emitted arithmetic
//! does not convert or increase input precision. [`RustEmitter::emit_batch`]
//! additionally emits instruction-major kernels with lazy active-lane branches
//! and shared helpers. Vectorized Dual graphs use flattened scalar components after
//! applying Symbolica's `Dualizer`; component order is left unchanged.

mod batch;

use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fmt::Write,
};
use symbolica::{
    atom::Symbol,
    domains::{
        float::{Complex, Float},
        rational::Rational,
    },
    evaluate::{ExportedInstructions, ExpressionEvaluator, Instruction, Slot},
};

/// Existing numerical runtime for emitted modules, with exact cached constants.
pub mod runtime {
    pub use crate::native::batch::{self, BatchExecutionPolicy, BatchWorkspace};
    pub use crate::native::{ConstantContext as Context, ConstantSpec, primitives};
    pub use primitives::{C, NativeFloat};
}

/// Generated Rust source and statistics for its scalar entry point.
/// `emit_batch` adds batch code to `source`; these counts still describe the
/// included scalar graph and its canonical helpers, rather than counting both.
#[derive(Clone, Debug)]
pub struct GeneratedModule {
    pub source: String,
    pub input_count: usize,
    pub output_count: usize,
    pub instruction_count: usize,
    pub helper_count: usize,
    pub outlined_branch_count: usize,
}

/// Emission settings; small graphs can opt into inlined complex arithmetic.
#[derive(Clone, Copy, Debug, Default)]
pub struct EmitOptions {
    pub inline_arithmetic: bool,
}

/// A generator sharing one exact constant table across emitted modules.
///
/// Failed exports leave the table unchanged. Previously generated indices stay
/// valid when more modules are emitted. Emit the table after the last module.
#[derive(Default)]
pub struct RustEmitter {
    constants: Constants,
}

impl RustEmitter {
    /// Export an exact evaluator, checking that coefficient mapping did not
    /// silently resolve missing constant metadata, then emit generic Rust.
    pub fn emit(
        &mut self,
        evaluator: ExpressionEvaluator<Complex<Rational>>,
        options: EmitOptions,
    ) -> Result<GeneratedModule, String> {
        let graph = export_checked(evaluator)?;
        self.emit_graph(&graph, options)
    }

    /// Emit scalar and instruction-major batch entry points together.
    ///
    /// Conditional branches run only their active lanes, including shared bodies.
    /// Exported labels, joins and forward control flow retain Symbolica's lazy
    /// semantics. Malformed control flow is rejected before constants are committed.
    /// The generated `evaluate_batch(input, output, rows, ctx, workspace)` accepts
    /// row-major buffers, validates dimensions before writing, and uses the scalar
    /// entry point when `rows == 1`. Larger batches run each instruction over its
    /// active lanes before continuing. Raw inputs retain their precision and uncertainty.
    /// Shared bodies currently require one output, as in the scalar emitter.
    pub fn emit_batch(
        &mut self,
        evaluator: ExpressionEvaluator<Complex<Rational>>,
        options: EmitOptions,
    ) -> Result<GeneratedModule, String> {
        let graph = export_checked(evaluator)?;
        self.emit_batch_graph(&graph, options)
    }

    fn emit_batch_graph(
        &mut self,
        graph: &Graph,
        options: EmitOptions,
    ) -> Result<GeneratedModule, String> {
        // Include both emissions in one transaction. Neither rejection nor a
        // failed helper export can alter previously assigned constant indices.
        let mut staged = Self {
            constants: self.constants.clone(),
        };
        batch::check_control_flow(graph)?;
        let mut module = staged.emit_graph(graph, options)?;
        module
            .source
            .push_str(&batch::emit(graph, &mut staged.constants)?);
        self.constants = staged.constants;
        Ok(module)
    }

    fn emit_graph(
        &mut self,
        graph: &Graph,
        options: EmitOptions,
    ) -> Result<GeneratedModule, String> {
        let mut constants = self.constants.clone();
        let mut generator = Generator {
            constants: &mut constants,
            helpers: Vec::new(),
            normalized: HashMap::new(),
            outlines: Outlines::default(),
        };
        let source = generator.module(graph, options.inline_arithmetic)?;
        let result = GeneratedModule {
            source,
            input_count: graph.input_count,
            output_count: graph.output_count,
            instruction_count: graph.instructions.len(),
            helper_count: generator.helpers.len(),
            outlined_branch_count: generator.outlines.functions.len(),
        };
        self.constants = constants;
        Ok(result)
    }

    /// Exact descriptors in the order referenced by emitted `ctx.constant(i)`.
    pub fn constants(&self) -> &[ExactConstant] {
        &self.constants.values
    }

    /// A child Rust module defining `CONSTANTS: &[ConstantSpec]`.
    pub fn constants_source(&self) -> String {
        self.constants.source()
    }
}

// Preserve Symbolica's exact registered-constant metadata when exporting.
// Numeric placeholders must never become literal zero constants in native code.
#[derive(Clone)]
struct Graph {
    input_count: usize,
    output_count: usize,
    temporary_count: usize,
    instructions: Vec<Instruction>,
    constants: Vec<Complex<Rational>>,
    constant_functions: Vec<ConstantFunction>,
    sub_evaluators: Vec<SubEvaluator>,
}
#[derive(Clone)]
struct ConstantFunction {
    index: usize,
    symbol: Symbol,
    tags: Vec<String>,
    fixed_args: Vec<Complex<Rational>>,
}
#[derive(Clone)]
struct SubEvaluator {
    symbol: Symbol,
    tags: Vec<String>,
    instructions: Graph,
}

impl From<ExportedInstructions<Complex<Rational>>> for Graph {
    fn from(graph: ExportedInstructions<Complex<Rational>>) -> Self {
        Self {
            input_count: graph.input_count,
            output_count: graph.output_count,
            temporary_count: graph.temporary_count,
            instructions: graph.instructions,
            constants: graph.constants,
            constant_functions: graph
                .constant_functions
                .into_iter()
                .map(|constant| ConstantFunction {
                    index: constant.index,
                    symbol: constant.symbol,
                    tags: constant.tags,
                    fixed_args: constant.fixed_args,
                })
                .collect(),
            sub_evaluators: graph
                .sub_evaluators
                .into_iter()
                .map(|child| SubEvaluator {
                    symbol: child.symbol,
                    tags: child.tags,
                    instructions: std::sync::Arc::unwrap_or_clone(child.instructions).into(),
                })
                .collect(),
        }
    }
}

fn export_checked(
    evaluator: symbolica::evaluate::ExpressionEvaluator<Complex<Rational>>,
) -> Result<Graph, String> {
    let exact = evaluator.export_instructions();
    let numeric = evaluator
        .map_coeff_with_prec(
            &|c| Complex::new(c.re.to_multi_prec_float(128), c.im.to_multi_prec_float(128)),
            128,
        )
        .export_instructions();
    fn check(
        exact: &ExportedInstructions<Complex<Rational>>,
        numeric: &ExportedInstructions<Complex<Float>>,
    ) -> Result<(), String> {
        if exact.constants.len() != numeric.constants.len()
            || exact.sub_evaluators.len() != numeric.sub_evaluators.len()
        {
            return Err("coefficient mapping changed exported graph shape".into());
        }
        for (index, (value, numeric)) in exact.constants.iter().zip(&numeric.constants).enumerate()
        {
            if exact
                .constant_functions
                .iter()
                .any(|constant| constant.index == index)
            {
                continue;
            }
            let exact = value;
            if exact.re.to_multi_prec_float(128) != numeric.re
                || exact.im.to_multi_prec_float(128) != numeric.im
            {
                return Err(
                    "unresolved symbolic constant: upstream export omitted its definition".into(),
                );
            }
        }
        for (exact, numeric) in exact.sub_evaluators.iter().zip(&numeric.sub_evaluators) {
            check(&exact.instructions, &numeric.instructions)?;
        }
        Ok(())
    }
    check(&exact, &numeric)?;
    Ok(exact.into())
}

type Environment = BTreeMap<Slot, String>;
type Signature = (Symbol, Vec<String>);

/// Exact, owned constant metadata; rational strings never pass through f64.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum ExactConstant {
    Rational(String, String),
    Pi,
    Polylog2(String, String),
}

impl ExactConstant {
    /// Borrow initialization data for the existing native constant context.
    pub fn as_spec(&self) -> runtime::ConstantSpec<'_> {
        match self {
            Self::Rational(re, im) => runtime::ConstantSpec::Rational { re, im },
            Self::Pi => runtime::ConstantSpec::Pi,
            Self::Polylog2(re, im) => runtime::ConstantSpec::Polylog2 { re, im },
        }
    }
}

#[derive(Clone, Default)]
struct Constants {
    values: Vec<ExactConstant>,
    indices: HashMap<ExactConstant, usize>,
}

impl Constants {
    fn insert(&mut self, value: ExactConstant) -> usize {
        if let Some(index) = self.indices.get(&value) {
            return *index;
        }
        let index = self.values.len();
        self.indices.insert(value.clone(), index);
        self.values.push(value);
        index
    }

    fn graph_constants(&mut self, graph: &Graph) -> Result<Vec<usize>, String> {
        let mut values = graph
            .constants
            .iter()
            .map(|value| ExactConstant::Rational(value.re.to_string(), value.im.to_string()))
            .collect::<Vec<_>>();
        for definition in &graph.constant_functions {
            let value = if definition.symbol == Symbol::PI
                && definition.tags.is_empty()
                && definition.fixed_args.is_empty()
            {
                ExactConstant::Pi
            } else if definition.symbol == symbolica::transcendental::polylog()
                && definition.tags == ["2"]
                && definition.fixed_args.len() == 1
            {
                let value = &definition.fixed_args[0];
                ExactConstant::Polylog2(value.re.to_string(), value.im.to_string())
            } else {
                return Err(format!(
                    "unsupported exact constant {} {:?} {:?}",
                    definition.symbol, definition.tags, definition.fixed_args
                ));
            };
            *values
                .get_mut(definition.index)
                .ok_or("constant definition index out of range")? = value;
        }
        Ok(values.into_iter().map(|value| self.insert(value)).collect())
    }

    fn source(&self) -> String {
        let mut output = header();
        output.push_str(
            "use super::ConstantSpec;\npub(super) const CONSTANTS: &[ConstantSpec] = &[\n",
        );
        for value in &self.values {
            match value {
                ExactConstant::Rational(re, im) => writeln!(
                    output,
                    "    ConstantSpec::Rational {{ re: {re:?}, im: {im:?} }},"
                )
                .unwrap(),
                ExactConstant::Pi => output.push_str("    ConstantSpec::Pi,\n"),
                ExactConstant::Polylog2(re, im) => writeln!(
                    output,
                    "    ConstantSpec::Polylog2 {{ re: {re:?}, im: {im:?} }},"
                )
                .unwrap(),
            }
        }
        output.push_str("];\n");
        output
    }
}

fn header() -> String {
    "// Generated by examples/generate_native.rs from repaired exact expressions.\n\
     // Regenerate; do not hand-edit. No bytecode, numerical regulator, or Fortran callback.\n"
        .into()
}

struct Generator<'a> {
    constants: &'a mut Constants,
    helpers: Vec<String>,
    normalized: HashMap<String, String>,
    outlines: Outlines,
}

impl<'a> Generator<'a> {
    fn module(&mut self, graph: &Graph, inline_arithmetic: bool) -> Result<String, String> {
        let body = self.body(graph, true, &[])?;
        let mut output = header();
        output.push_str("#![allow(unused_imports, unused_variables)]\nuse super::{C, Context, NativeFloat};\nuse super::primitives::*;\n");
        if inline_arithmetic {
            output.push_str("use super::primitives::{inline_add as add, inline_mul as mul};\n");
        }
        output.push('\n');
        writeln!(output, "pub(super) fn evaluate<T: NativeFloat>(input: &[C<T>], ctx: &Context<T>) -> [C<T>; {}] {{", graph.output_count).unwrap();
        output.push_str(&body);
        output.push_str("}\n\n");
        for helper in &self.helpers {
            output.push_str(helper);
        }
        for function in &self.outlines.functions {
            output.push_str(function);
        }
        Ok(output)
    }

    /// `bound[i]` is the global constant every call passes as parameter `i`;
    /// such parameters are read from `ctx` and left out of the signature.
    fn helper(&mut self, graph: &Graph, bound: &[Option<usize>]) -> Result<Callee, String> {
        if graph.output_count != 1 {
            return Err("native helper must have one output".into());
        }
        let body = self.body(graph, false, bound)?;
        let constant = (0..graph.input_count)
            .map(|index| bound.get(index).copied().flatten().is_some())
            .collect::<Vec<_>>();
        // Child names are already canonical, so body+signature identify the
        // implementation even when Symbolica assigned different capture tags.
        let normalized = format!(
            "{}:{}:{constant:?}:{body}",
            graph.input_count, graph.output_count
        );
        if let Some(name) = self.normalized.get(&normalized) {
            return Ok(Callee {
                name: name.clone(),
                constant,
            });
        }
        let name = format!("helper_{}", self.helpers.len());
        let parameters = (0..graph.input_count)
            .filter(|index| !constant[*index])
            .map(|index| format!("a{index}: &C<T>"))
            .chain(std::iter::once("ctx: &Context<T>".into()))
            .collect::<Vec<_>>()
            .join(", ");
        self.helpers.push(format!(
            "fn {name}<T: NativeFloat>({parameters}) -> C<T> {{\n{body}}}\n\n"
        ));
        self.normalized.insert(normalized, name.clone());
        Ok(Callee { name, constant })
    }

    fn body(
        &mut self,
        graph: &Graph,
        entry: bool,
        bound: &[Option<usize>],
    ) -> Result<String, String> {
        let constants = self.constants.graph_constants(graph)?;
        let mut callees = HashMap::<Signature, Callee>::new();
        for child in &graph.sub_evaluators {
            let callee = self.helper(
                &child.instructions,
                &constant_arguments(graph, &constants, bound, child),
            )?;
            let signature = (child.symbol, child.tags.clone());
            if let Some(previous) = callees.insert(signature, callee.clone())
                && previous != callee
            {
                return Err("ambiguous exported helper signature".into());
            }
        }
        let labels = labels(&graph.instructions)?;
        let live = liveness(&graph.instructions, graph.output_count, &labels)?;
        let mut emitter = Emitter {
            instructions: &graph.instructions,
            labels,
            live,
            callees,
            next: 0,
            computed: HashMap::new(),
            computed_log: Vec::new(),
            outlines: &mut self.outlines,
        };
        let mut environment = Environment::new();
        let mut output = String::new();
        for index in 0..graph.input_count {
            if entry {
                writeln!(output, "    let a{index} = &input[{index}];").unwrap();
            }
            let value = match bound.get(index).copied().flatten() {
                Some(global) => format!("ctx.constant({global})"),
                None => format!("a{index}"),
            };
            environment.insert(Slot::Param(index), value);
        }
        for (index, global) in constants.into_iter().enumerate() {
            environment.insert(Slot::Const(index), format!("ctx.constant({global})"));
        }
        emitter.range(
            0,
            graph.instructions.len(),
            &mut environment,
            &mut output,
            1,
        )?;
        let results = (0..graph.output_count)
            .map(|index| {
                read(&environment, Slot::Out(index)).map(|value| format!("(*({value})).clone()"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        if entry {
            writeln!(output, "    [{}]", results.join(", ")).unwrap();
        } else {
            writeln!(output, "    {}", results[0]).unwrap();
        }
        Ok(output)
    }
}

/// A shared body called from a graph, as its call sites see it.
#[derive(Clone, Debug, PartialEq)]
struct Callee {
    name: String,
    /// Parameters bound to one constant at every call; calls omit them.
    constant: Vec<bool>,
}

impl From<&str> for Callee {
    fn from(name: &str) -> Self {
        Self {
            name: name.into(),
            constant: Vec::new(),
        }
    }
}

/// Symbolica passes the constants a shared body uses as arguments, through
/// every enclosing body. For each parameter of `child`, the global constant
/// that all of its calls in `graph` pass (directly or through a parameter of
/// `graph` that is itself bound), so the body can read it from `ctx`.
fn constant_arguments(
    graph: &Graph,
    constants: &[usize],
    bound: &[Option<usize>],
    child: &SubEvaluator,
) -> Vec<Option<usize>> {
    let mut common: Option<Vec<Option<usize>>> = None;
    for instruction in &graph.instructions {
        let Instruction::Fun(_, function, _) = instruction else {
            continue;
        };
        let (symbol, tags, arguments) = &**function;
        if *symbol != child.symbol || *tags != child.tags {
            continue;
        }
        let values = arguments
            .iter()
            .map(|slot| match slot {
                Slot::Const(index) => constants.get(*index).copied(),
                Slot::Param(index) => bound.get(*index).copied().flatten(),
                _ => None,
            })
            .collect::<Vec<_>>();
        common = Some(match common {
            None => values,
            Some(previous) if previous.len() == values.len() => previous
                .into_iter()
                .zip(values)
                .map(|(first, second)| if first == second { first } else { None })
                .collect(),
            Some(_) => return Vec::new(),
        });
    }
    match common {
        Some(values) if values.len() == child.instructions.input_count => values,
        _ => Vec::new(),
    }
}

/// Outline complete branch bodies, after their children have been outlined.
/// This is ordinary function extraction: captured values stay borrowed, and
/// calls remain inside the original selected arm. No arithmetic is reordered.
struct Outlines {
    minimum_lines: usize,
    functions: Vec<String>,
    normalized: HashMap<String, String>,
}

impl Default for Outlines {
    fn default() -> Self {
        Self {
            minimum_lines: 128,
            functions: Vec::new(),
            normalized: HashMap::new(),
        }
    }
}

fn value_identifier(value: &str) -> bool {
    matches!(value.as_bytes().first(), Some(b'a' | b'v'))
        && value.len() > 1
        && value.as_bytes()[1..].iter().all(u8::is_ascii_digit)
}

fn lexical_tokens(source: &str) -> Vec<&str> {
    let bytes = source.as_bytes();
    let mut tokens = Vec::new();
    let mut start = 0;
    while start < bytes.len() {
        let mut end = start + 1;
        if bytes[start].is_ascii_alphanumeric() || bytes[start] == b'_' {
            while end < bytes.len() && (bytes[end].is_ascii_alphanumeric() || bytes[end] == b'_') {
                end += 1;
            }
        } else if bytes[start].is_ascii_whitespace() {
            while end < bytes.len() && bytes[end].is_ascii_whitespace() {
                end += 1;
            }
        }
        tokens.push(&source[start..end]);
        start = end;
    }
    tokens
}

/// Alpha-normalize only values DEFINED inside this conditional. Captured SSA
/// names remain literal identities: two equal-shaped arms with different
/// inputs must never be confused. The local namespace cannot collide with
/// emitted `aN` parameters or captured `vN` values.
fn conditional_key(source: &str, outputs: usize) -> String {
    let tokens = lexical_tokens(source);
    let mut definitions = HashMap::new();
    let mut binding = false;
    for token in &tokens {
        if *token == "let" {
            binding = true;
        } else if *token == "=" {
            binding = false;
        } else if binding && value_identifier(token) {
            let index = definitions.len();
            definitions
                .entry(*token)
                .or_insert_with(|| format!("local_{index}"));
        }
    }
    let mut key = format!("conditional:{outputs}:");
    for token in tokens {
        if token.trim().is_empty() {
            key.push(' ');
        } else if let Some(local) = definitions.get(token) {
            key.push_str(local);
        } else {
            key.push_str(token);
        }
    }
    key
}

impl Outlines {
    fn factor(&mut self, body: &str, outputs: usize) -> Result<Option<String>, String> {
        if body.lines().count() < self.minimum_lines {
            return Ok(None);
        }
        let tokens = lexical_tokens(body);
        let mut definitions = BTreeMap::new();
        let mut binding = false;
        for token in &tokens {
            if *token == "let" {
                binding = true;
            } else if *token == "=" {
                binding = false;
            } else if binding && value_identifier(token) {
                let next = definitions.len();
                definitions
                    .entry((*token).to_owned())
                    .or_insert_with(|| format!("v{next}"));
            }
        }
        let mut captures = Vec::<String>::new();
        let mut capture_names = HashMap::<String, usize>::new();
        for token in &tokens {
            if value_identifier(token)
                && !definitions.contains_key(*token)
                && !capture_names.contains_key(*token)
            {
                capture_names.insert((*token).to_owned(), captures.len());
                captures.push((*token).to_owned());
            }
        }
        let mut normalized = String::new();
        let mut index = 0;
        while index < tokens.len() {
            let token = tokens[index];
            if token == "&" {
                let mut next = index + 1;
                while next < tokens.len() && tokens[next].trim().is_empty() {
                    next += 1;
                }
                if let Some(capture) = tokens.get(next).and_then(|token| capture_names.get(*token))
                    && tokens[next].starts_with('v')
                {
                    // An outer SSA local is owned (`&v42`), whereas a formal
                    // parameter is already borrowed (`a3`). Both become one
                    // borrowed outline parameter, with no extra clone.
                    write!(normalized, "a{capture}").unwrap();
                    index = next + 1;
                    continue;
                }
            }
            if let Some(name) = definitions.get(token) {
                normalized.push_str(name);
            } else if let Some(capture) = capture_names.get(token) {
                if token.starts_with('v') {
                    return Err(format!("unborrowed captured SSA local {token}"));
                }
                write!(normalized, "a{capture}").unwrap();
            } else {
                normalized.push_str(token);
            }
            index += 1;
        }
        // Whitespace has no semantics in this generated subset. Removing the
        // original nesting indentation makes alpha-equivalent arms identical.
        let normalized = normalized
            .lines()
            .map(str::trim)
            .collect::<Vec<_>>()
            .join("\n");
        let key = format!("{}:{outputs}:{normalized}", captures.len());
        let name = if let Some(name) = self.normalized.get(&key) {
            name.clone()
        } else {
            let name = format!("branch_{}", self.functions.len());
            let arguments = (0..captures.len())
                .map(|index| format!("a{index}: &C<T>"))
                .chain(std::iter::once("ctx: &Context<T>".into()))
                .collect::<Vec<_>>()
                .join(", ");
            let result = if outputs == 1 {
                "C<T>".into()
            } else {
                format!("({})", vec!["C<T>"; outputs].join(", "))
            };
            let body = normalized
                .lines()
                .map(|line| format!("    {line}\n"))
                .collect::<String>();
            self.functions.push(format!(
                "fn {name}<T: NativeFloat>({arguments}) -> {result} {{\n{body}}}\n\n"
            ));
            self.normalized.insert(key, name.clone());
            name
        };
        let arguments = captures
            .into_iter()
            .map(|name| {
                if name.starts_with('v') {
                    format!("&{name}")
                } else {
                    name
                }
            })
            .chain(std::iter::once("ctx".into()))
            .collect::<Vec<_>>()
            .join(", ");
        Ok(Some(format!("{name}({arguments})")))
    }
}

fn labels(instructions: &[Instruction]) -> Result<HashMap<usize, usize>, String> {
    let mut result = HashMap::new();
    for (index, instruction) in instructions.iter().enumerate() {
        if let Instruction::Label(label) = instruction
            && result.insert(*label, index).is_some()
        {
            return Err("duplicate instruction label".into());
        }
    }
    Ok(result)
}

fn uses(instruction: &Instruction) -> (Option<Slot>, Vec<Slot>) {
    match instruction {
        Instruction::Add(out, args, _) | Instruction::Mul(out, args, _) => {
            (Some(*out), args.clone())
        }
        Instruction::Pow(out, arg, _, _) | Instruction::Assign(out, arg) => {
            (Some(*out), vec![*arg])
        }
        Instruction::Powf(out, base, exponent, _) => (Some(*out), vec![*base, *exponent]),
        Instruction::Fun(out, data, _) => (Some(*out), data.2.clone()),
        Instruction::IfElse(condition, _) => (None, vec![*condition]),
        Instruction::Join(out, condition, yes, no) => (Some(*out), vec![*condition, *yes, *no]),
        Instruction::Goto(_) | Instruction::Label(_) => (None, Vec::new()),
    }
}

fn liveness(
    instructions: &[Instruction],
    outputs: usize,
    labels: &HashMap<usize, usize>,
) -> Result<Vec<BTreeSet<Slot>>, String> {
    // A phi operand is live only on its incoming branch, not throughout both
    // branches. In particular, an outer false-arm value must not become live
    // inside a nested true-arm conditional before that value even exists.
    let mut phi_edges = HashMap::new();
    let mut structured_joins = BTreeSet::new();
    for instruction in instructions {
        if let Instruction::IfElse(condition, false_label) = instruction {
            let false_start = *labels.get(false_label).ok_or("missing false label")?;
            let true_end = false_start.checked_sub(1).ok_or("invalid false label")?;
            if let Instruction::Goto(end_label) = &instructions[true_end] {
                let end = *labels.get(end_label).ok_or("missing end label")?;
                let mut phis = Vec::new();
                let mut join = end + 1;
                while let Some(Instruction::Join(out, join_condition, yes, no)) =
                    instructions.get(join)
                {
                    if join_condition != condition {
                        break;
                    }
                    phis.push((*out, *yes, *no));
                    structured_joins.insert(join);
                    join += 1;
                }
                if !phis.is_empty() {
                    phi_edges.insert(end, (true_end, phis));
                }
            }
        }
    }
    let mut live = vec![BTreeSet::new(); instructions.len() + 1];
    live[instructions.len()].extend((0..outputs).map(Slot::Out));
    for index in (0..instructions.len()).rev() {
        let edge = |destination: usize| -> Result<BTreeSet<Slot>, String> {
            let mut values = live[destination].clone();
            if let Some((true_end, phis)) = phi_edges.get(&destination) {
                let true_edge = index == *true_end;
                if !true_edge && index + 1 != destination {
                    return Err("non-structured phi predecessor".into());
                }
                // Joins execute in order and may overwrite a slot consumed by
                // a later join. Propagate each selected arm backwards through
                // all assignments instead of making sibling-arm values live.
                for (out, yes, no) in phis.iter().rev() {
                    values.remove(out);
                    values.insert(if true_edge { *yes } else { *no });
                }
            }
            Ok(values)
        };
        let mut current = match &instructions[index] {
            Instruction::Goto(label) => {
                let destination = *labels.get(label).ok_or("unknown jump label")?;
                if destination <= index {
                    return Err("loops are not part of the scalar expression DAG".into());
                }
                edge(destination)?
            }
            Instruction::IfElse(_, label) => {
                let destination = *labels.get(label).ok_or("unknown conditional label")?;
                if destination <= index {
                    return Err("backward conditional jump".into());
                }
                edge(index + 1)?
                    .union(&edge(destination)?)
                    .copied()
                    .collect()
            }
            _ => edge(index + 1)?,
        };
        let (defined, used) = if let Instruction::Join(out, condition, _, _) = &instructions[index]
            && structured_joins.contains(&index)
        {
            (Some(*out), vec![*condition])
        } else {
            uses(&instructions[index])
        };
        if let Some(slot) = defined {
            current.remove(&slot);
        }
        current.extend(used);
        live[index] = current;
    }
    Ok(live)
}

fn read(environment: &Environment, slot: Slot) -> Result<&str, String> {
    environment
        .get(&slot)
        .map(String::as_str)
        .ok_or_else(|| format!("unbound native slot {slot}"))
}

struct Emitter<'a, 'b> {
    instructions: &'a [Instruction],
    labels: HashMap<usize, usize>,
    live: Vec<BTreeSet<Slot>>,
    callees: HashMap<Signature, Callee>,
    next: usize,
    /// Only bindings that dominate the current instruction are reusable.
    /// Branch-local bindings never escape except through explicit phi tuples.
    computed: HashMap<String, Vec<String>>,
    computed_log: Vec<String>,
    outlines: &'b mut Outlines,
}

impl Emitter<'_, '_> {
    fn fresh(&mut self) -> String {
        let name = format!("v{}", self.next);
        self.next += 1;
        name
    }

    fn emit(&mut self, expression: String, output: &mut String, indent: &str) -> String {
        if let Some(references) = self.computed.get(&expression) {
            return references[0].clone();
        }
        let name = self.fresh();
        writeln!(output, "{indent}let {name} = {expression};").unwrap();
        let reference = format!("&{name}");
        self.remember(expression, vec![reference.clone()]);
        reference
    }

    fn remember(&mut self, expression: String, references: Vec<String>) {
        debug_assert!(!self.computed.contains_key(&expression));
        self.computed.insert(expression.clone(), references);
        self.computed_log.push(expression);
    }

    fn restore_computed(&mut self, checkpoint: usize) {
        while self.computed_log.len() > checkpoint {
            let expression = self.computed_log.pop().unwrap();
            self.computed.remove(&expression);
        }
    }

    fn range(
        &mut self,
        start: usize,
        end: usize,
        environment: &mut Environment,
        output: &mut String,
        depth: usize,
    ) -> Result<(), String> {
        let indent = "    ".repeat(depth);
        let mut index = start;
        while index < end {
            let instruction = &self.instructions[index];
            if let Instruction::IfElse(condition, false_label) = instruction {
                let condition = *condition;
                let false_start = *self.labels.get(false_label).ok_or("missing false branch")?;
                let Some(Instruction::Goto(end_label)) = false_start
                    .checked_sub(1)
                    .and_then(|i| self.instructions.get(i))
                else {
                    return Err("conditional has no structured true-branch terminator".into());
                };
                let join_label = *self.labels.get(end_label).ok_or("missing join label")?;
                let join = join_label + 1;
                let Some(Instruction::Join(_, join_condition, _, _)) = self.instructions.get(join)
                else {
                    return Err("structured branch does not terminate in Join".into());
                };
                if *join_condition != condition || false_start <= index || join >= end {
                    return Err("invalid structured branch boundary".into());
                }
                let condition_value = read(environment, condition)?.to_owned();
                let mut true_environment = environment.clone();
                let mut false_environment = environment.clone();
                let mut true_body = String::new();
                let mut false_body = String::new();
                let dominating = self.computed_log.len();
                self.range(
                    index + 1,
                    false_start - 1,
                    &mut true_environment,
                    &mut true_body,
                    depth + 1,
                )?;
                self.restore_computed(dominating);
                self.range(
                    false_start + 1,
                    join_label,
                    &mut false_environment,
                    &mut false_body,
                    depth + 1,
                )?;
                self.restore_computed(dominating);
                let mut merged = Vec::new();
                let mut after_join = join;
                while let Some(Instruction::Join(target, join_condition, yes, no)) =
                    self.instructions.get(after_join)
                {
                    if *join_condition != condition || after_join >= end {
                        break;
                    }
                    let yes = read(&true_environment, *yes)?.to_owned();
                    let no = read(&false_environment, *no)?.to_owned();
                    true_environment.insert(*target, yes.clone());
                    false_environment.insert(*target, no.clone());
                    merged.push((*target, yes, no));
                    after_join += 1;
                }
                // Register allocation may reuse slots across a branch. Merge
                // every additional live value rather than assuming only Join
                // escapes. Branch-private dead values never need zero padding.
                for slot in &self.live[after_join] {
                    if merged.iter().any(|(target, _, _)| target == slot) {
                        continue;
                    }
                    let yes = read(&true_environment, *slot)?;
                    let no = read(&false_environment, *slot)?;
                    if yes != no {
                        merged.push((*slot, yes.into(), no.into()));
                    } else {
                        // Assign aliases can overwrite the old slot with
                        // the same dominating value in both branches.
                        environment.insert(*slot, yes.to_owned());
                    }
                }
                let names = merged.iter().map(|_| self.fresh()).collect::<Vec<_>>();
                let binding = if names.len() == 1 {
                    names[0].clone()
                } else {
                    format!("({})", names.join(", "))
                };
                let branch_result = |yes: bool| {
                    let values = merged
                        .iter()
                        .map(|(_, first, second)| {
                            format!("(*({})).clone()", if yes { first } else { second })
                        })
                        .collect::<Vec<_>>();
                    if values.len() == 1 {
                        values[0].clone()
                    } else {
                        format!("({})", values.join(", "))
                    }
                };
                writeln!(true_body, "{indent}    {}", branch_result(true)).unwrap();
                writeln!(false_body, "{indent}    {}", branch_result(false)).unwrap();
                if let Some(call) = self.outlines.factor(&true_body, merged.len())? {
                    true_body = format!("{indent}    {call}\n");
                }
                if let Some(call) = self.outlines.factor(&false_body, merged.len())? {
                    false_body = format!("{indent}    {call}\n");
                }
                let expression = format!(
                    "if truth({condition_value}) {{\n{true_body}{indent}}} else {{\n{false_body}{indent}}}"
                );
                let key = conditional_key(&expression, merged.len());
                if let Some(references) = self.computed.get(&key) {
                    debug_assert_eq!(references.len(), merged.len());
                    for ((slot, _, _), reference) in merged.iter().zip(references) {
                        environment.insert(*slot, reference.clone());
                    }
                    index = after_join;
                    continue;
                }
                writeln!(
                    output,
                    "{indent}let {binding} = if truth({condition_value}) {{"
                )
                .unwrap();
                output.push_str(&true_body);
                writeln!(output, "{indent}}} else {{").unwrap();
                output.push_str(&false_body);
                writeln!(output, "{indent}}};").unwrap();
                let references = names
                    .iter()
                    .map(|name| format!("&{name}"))
                    .collect::<Vec<_>>();
                for ((slot, _, _), reference) in merged.into_iter().zip(&references) {
                    environment.insert(slot, reference.clone());
                }
                self.remember(key, references);
                index = after_join;
                continue;
            }
            let (target, expression) = match instruction {
                Instruction::Add(target, arguments, _) | Instruction::Mul(target, arguments, _) => {
                    let operation = if matches!(instruction, Instruction::Add(..)) {
                        "add"
                    } else {
                        "mul"
                    };
                    let Some(first) = arguments.first() else {
                        return Err("empty arithmetic instruction".into());
                    };
                    let mut value = read(environment, *first)?.to_owned();
                    for argument in &arguments[1..] {
                        value = self.emit(
                            format!("{operation}({value}, {})", read(environment, *argument)?),
                            output,
                            &indent,
                        );
                    }
                    // A slot assignment aliases an immutable SSA value. No
                    // new numerical operation or precision change is needed.
                    environment.insert(*target, value);
                    index += 1;
                    continue;
                }
                Instruction::Pow(target, base, exponent, is_real) => (
                    *target,
                    format!("powi({}, {exponent}, {is_real})", read(environment, *base)?),
                ),
                Instruction::Powf(target, base, exponent, is_real) => (
                    *target,
                    format!(
                        "powf({}, {}, {is_real})",
                        read(environment, *base)?,
                        read(environment, *exponent)?
                    ),
                ),
                Instruction::Assign(target, value) => {
                    let value = read(environment, *value)?.to_owned();
                    environment.insert(*target, value);
                    index += 1;
                    continue;
                }
                Instruction::Fun(target, function, _) => {
                    let (symbol, tags, arguments) = &**function;
                    let callee = self.callees.get(&(*symbol, tags.clone()));
                    let values = arguments
                        .iter()
                        .enumerate()
                        .filter(|(index, _)| {
                            !callee.is_some_and(|callee| {
                                callee.constant.get(*index).copied().unwrap_or(false)
                            })
                        })
                        .map(|(_, slot)| read(environment, *slot))
                        .collect::<Result<Vec<_>, _>>()?;
                    let call = if let Some(Callee { name, .. }) = callee {
                        format!(
                            "{name}({}ctx)",
                            if values.is_empty() {
                                String::new()
                            } else {
                                format!("{}, ", values.join(", "))
                            }
                        )
                    } else {
                        let name = primitive_name(*symbol, tags)?;
                        if values.len() != 1 {
                            return Err(format!("unexpected arity for {name}"));
                        }
                        format!("{name}({})", values[0])
                    };
                    (*target, call)
                }
                Instruction::Join(target, condition, yes, no) => (
                    *target,
                    format!(
                        "if truth({}) {{ (*({})).clone() }} else {{ (*({})).clone() }}",
                        read(environment, *condition)?,
                        read(environment, *yes)?,
                        read(environment, *no)?
                    ),
                ),
                Instruction::Label(_) | Instruction::Goto(_) => {
                    return Err("unstructured label outside a native if".into());
                }
                Instruction::IfElse(..) => unreachable!(),
            };
            let value = self.emit(expression, output, &indent);
            environment.insert(target, value);
            index += 1;
        }
        Ok(())
    }
}

// Keep primitive support and tagged-function matching identical in both emitters.
fn primitive_name(symbol: Symbol, tags: &[String]) -> Result<&'static str, String> {
    if symbol == Symbol::LOG {
        Ok("log")
    } else if symbol == Symbol::SQRT {
        Ok("sqrt")
    } else if symbol == Symbol::ABS {
        Ok("abs")
    } else if symbol == Symbol::CONJ {
        Ok("conj")
    } else if symbol == symbolica::transcendental::polylog() && tags == ["2"] {
        Ok("polylog2")
    } else {
        Err(format!("unsupported native function {symbol} {tags:?}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::atom::{Atom, AtomCore};
    use symbolica::transcendental::TranscendentalFunctions;

    #[test]
    fn public_emitter_supports_non_master_output_counts_and_shared_constants() {
        crate::initialization::ensure_symbolica_state();
        let x = symbolica::symbol!("rust_emitter_public_x").to_atom();
        let mut emitter = RustEmitter::default();
        for outputs in [1, 2, 7] {
            let expressions: Vec<_> = (1..=outputs).map(|power| x.pow(power as i64)).collect();
            let evaluator = Atom::evaluator_multiple(&expressions, std::slice::from_ref(&x))
                .cores(1)
                .build()
                .unwrap();
            let module = emitter.emit(evaluator, EmitOptions::default()).unwrap();
            assert_eq!(module.input_count, 1);
            assert_eq!(module.output_count, outputs);
            assert!(module.source.contains(&format!("-> [C<T>; {outputs}]")));
            let source = format!(
                "{PURE_NUMERIC_PROBE}\nmod primitives {{ pub(super) use super::{{add,mul,powi,truth}}; }}\nmod generated {{\n{}}}\nfn main() {{ let ctx = Context(std::marker::PhantomData); assert_eq!(generated::evaluate(&[2.0], &ctx), {:?}); }}",
                module.source,
                (1..=outputs)
                    .map(|power| 2_f64.powi(power as i32))
                    .collect::<Vec<_>>()
            );
            compile_probe(&source);
        }
        let constant = Atom::num((7, 13));
        for _ in 0..2 {
            let evaluator =
                Atom::evaluator_multiple(std::slice::from_ref(&constant), &[] as &[Atom])
                    .cores(1)
                    .build()
                    .unwrap();
            emitter.emit(evaluator, EmitOptions::default()).unwrap();
        }
        assert_eq!(
            emitter
                .constants()
                .iter()
                .filter(|constant| **constant == ExactConstant::Rational("7/13".into(), "0".into()))
                .count(),
            1
        );
        assert!(emitter.constants_source().contains("re: \"7/13\""));
    }

    #[test]
    fn failed_export_does_not_change_previously_generated_constant_indices() {
        let mut emitter = RustEmitter::default();
        let mut graph = Graph {
            input_count: 0,
            output_count: 1,
            temporary_count: 100,
            constants: vec![Complex::new(Rational::from(7), Rational::from(0))],
            constant_functions: Vec::new(),
            sub_evaluators: Vec::new(),
            instructions: vec![Instruction::Assign(Slot::Out(0), Slot::Const(0))],
        };
        emitter.emit_graph(&graph, EmitOptions::default()).unwrap();
        let before = emitter.constants_source();
        graph.constants[0] = Complex::new(Rational::from(11), Rational::from(0));
        // The constant table is populated before the instruction error occurs.
        graph.instructions = vec![Instruction::Assign(Slot::Out(0), Slot::Temp(42))];
        assert!(emitter.emit_graph(&graph, EmitOptions::default()).is_err());
        assert_eq!(emitter.constants_source(), before);
        graph.instructions = vec![Instruction::Assign(Slot::Out(0), Slot::Const(0))];
        let module = emitter.emit_graph(&graph, EmitOptions::default()).unwrap();
        assert!(module.source.contains("ctx.constant(1)"));
        assert_eq!(emitter.constants().len(), 2);
    }

    #[test]
    fn constant_arguments_are_read_from_the_context_through_nested_helpers() {
        crate::initialization::ensure_symbolica_state();
        let [outer, inner] = ["rust_emitter_constant_outer", "rust_emitter_constant_inner"]
            .map(|name| symbolica::symbol!(name));
        let number = |value: i64| Complex::new(Rational::from(value), Rational::from(0));
        let call = |target, symbol, arguments| {
            Instruction::Fun(target, Box::new((symbol, Vec::new(), arguments)), false)
        };
        let graph = |input_count, constants, instructions, sub_evaluators| Graph {
            input_count,
            output_count: 1,
            temporary_count: 2,
            instructions,
            constants,
            constant_functions: Vec::new(),
            sub_evaluators,
        };
        // inner(x, c) = x c; outer(x, c, d) = inner(x, c) + d. Every call
        // passes 7 as `c`, while `d` is 11 at one call and 7 at the other.
        let inner_body = graph(
            2,
            Vec::new(),
            vec![Instruction::Mul(Slot::Out(0), vec![Slot::Param(0), Slot::Param(1)], 0)],
            Vec::new(),
        );
        let outer_body = graph(
            3,
            Vec::new(),
            vec![
                call(Slot::Temp(0), inner, vec![Slot::Param(0), Slot::Param(1)]),
                Instruction::Add(Slot::Out(0), vec![Slot::Temp(0), Slot::Param(2)], 0),
            ],
            vec![SubEvaluator {
                symbol: inner,
                tags: Vec::new(),
                instructions: inner_body,
            }],
        );
        let entry = graph(
            1,
            vec![number(7), number(11)],
            vec![
                call(Slot::Temp(0), outer, vec![Slot::Param(0), Slot::Const(0), Slot::Const(1)]),
                call(Slot::Temp(1), outer, vec![Slot::Param(0), Slot::Const(0), Slot::Const(0)]),
                Instruction::Add(Slot::Out(0), vec![Slot::Temp(0), Slot::Temp(1)], 0),
            ],
            vec![SubEvaluator {
                symbol: outer,
                tags: Vec::new(),
                instructions: outer_body,
            }],
        );
        let mut emitter = RustEmitter::default();
        let module = emitter.emit_graph(&entry, EmitOptions::default()).unwrap();
        let source = &module.source;
        assert!(source.contains("fn helper_0<T: NativeFloat>(a0: &C<T>, ctx: &Context<T>)"), "{source}");
        assert!(source.contains("mul(a0, ctx.constant(0))"), "{source}");
        assert!(
            source.contains("fn helper_1<T: NativeFloat>(a0: &C<T>, a2: &C<T>, ctx: &Context<T>)"),
            "{source}"
        );
        assert!(source.contains("helper_0(a0, ctx)"), "{source}");
        assert!(source.contains("helper_1(a0, ctx.constant(1), ctx)"), "{source}");
        assert!(source.contains("helper_1(a0, ctx.constant(0), ctx)"), "{source}");
        let probe = PURE_NUMERIC_PROBE.replace(
            "struct Context<T>(std::marker::PhantomData<T>);",
            "struct Context<T>(Vec<T>);\nimpl<T> Context<T> { fn constant(&self, index: usize) -> &T { &self.0[index] } }",
        );
        compile_probe(&format!(
            "{probe}\nmod primitives {{ pub(super) use super::{{add,mul,powi,truth}}; }}\nmod generated {{\n{source}}}\nfn main() {{ let ctx = Context(vec![7.0, 11.0]); assert_eq!(generated::evaluate(&[2.0], &ctx), [46.0]); }}"
        ));
    }

    pub(super) const PURE_NUMERIC_PROBE: &str = r#"
#![allow(unused_variables, dead_code)]
trait NativeFloat: Copy + PartialEq + From<f64> + std::ops::Add<Output=Self> + std::ops::Mul<Output=Self> {
    fn integer_power(self, exponent: i32) -> Self;
}
impl NativeFloat for f64 {
    fn integer_power(self, exponent: i32) -> Self { assert!(self != 0.0 || exponent >= 0); self.powi(exponent) }
}
type C<T> = T;
struct Context<T>(std::marker::PhantomData<T>);
static POW_CALLS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
fn truth<T: NativeFloat>(a: &T) -> bool { *a != T::from(0.0) }
fn add<T: NativeFloat>(a: &T, b: &T) -> T { *a + *b }
fn mul<T: NativeFloat>(a: &T, b: &T) -> T { *a * *b }
fn powi<T: NativeFloat>(a: &T, exponent: i64, _: bool) -> T {
    POW_CALLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    a.integer_power(exponent as i32)
}
"#;

    #[test]
    fn registered_constants_are_exported_symbolically_not_as_zero_placeholders() {
        std::thread::Builder::new()
            .stack_size(128 * 1024 * 1024)
            .spawn(|| {
                crate::initialization::ensure_symbolica_state();
                let expressions = [Symbol::PI.to_atom(), Atom::num((1, 2)).polylog(2)];
                let evaluator = Atom::evaluator_multiple(&expressions, &[] as &[Atom])
                    .direct_translation(true)
                    .build()
                    .unwrap();
                let graph = export_checked(evaluator).unwrap();
                assert_eq!(graph.input_count, 0);
                assert_eq!(graph.constant_functions.len(), 2);
                let mut constants = Constants::default();
                constants.graph_constants(&graph).unwrap();
                assert!(constants.values.contains(&ExactConstant::Pi));
                assert!(
                    constants
                        .values
                        .contains(&ExactConstant::Polylog2("1/2".into(), "0".into()))
                );
            })
            .unwrap()
            .join()
            .unwrap();
    }

    #[test]
    fn consecutive_joins_keep_edge_liveness_and_sequential_slot_aliases() {
        let graph = Graph {
            input_count: 3,
            output_count: 2,
            temporary_count: 3,
            constants: vec![],
            constant_functions: vec![],
            sub_evaluators: vec![],
            instructions: vec![
                Instruction::IfElse(Slot::Param(0), 1),
                Instruction::Pow(Slot::Temp(0), Slot::Param(1), -1, false),
                Instruction::Assign(Slot::Temp(2), Slot::Param(1)),
                Instruction::Goto(2),
                Instruction::Label(1),
                Instruction::Pow(Slot::Temp(1), Slot::Param(2), -1, false),
                Instruction::Assign(Slot::Temp(2), Slot::Param(2)),
                Instruction::Label(2),
                Instruction::Join(Slot::Temp(0), Slot::Param(0), Slot::Temp(0), Slot::Temp(1)),
                Instruction::Join(Slot::Out(0), Slot::Param(0), Slot::Temp(0), Slot::Temp(0)),
                Instruction::Join(Slot::Temp(2), Slot::Param(0), Slot::Temp(2), Slot::Temp(2)),
                Instruction::Join(Slot::Out(1), Slot::Param(0), Slot::Temp(2), Slot::Temp(2)),
            ],
        };
        let labels = labels(&graph.instructions).unwrap();
        let live = liveness(&graph.instructions, graph.output_count, &labels).unwrap();
        assert!(!live[1].contains(&Slot::Temp(1)));
        assert!(!live[5].contains(&Slot::Temp(0)));
        let module = RustEmitter::default()
            .emit_graph(&graph, EmitOptions::default())
            .unwrap();
        let source = format!(
            "{PURE_NUMERIC_PROBE}\nmod primitives {{ pub(super) use super::{{add,mul,powi,truth}}; }}\nmod generated {{\n{}\n}}\nfn main() {{ let ctx = Context(std::marker::PhantomData); assert_eq!(generated::evaluate(&[1.,2.,0.],&ctx),[0.5,2.]); assert_eq!(generated::evaluate(&[0.,0.,4.],&ctx),[0.25,4.]); assert_eq!(POW_CALLS.load(std::sync::atomic::Ordering::Relaxed),2); }}",
            module.source
        );
        compile_probe(&source);
    }

    #[test]
    fn native_branches_merge_live_slots_without_eager_evaluation() {
        let instructions = vec![
            Instruction::IfElse(Slot::Param(0), 1),
            Instruction::Pow(Slot::Temp(0), Slot::Param(1), -1, false),
            Instruction::Assign(Slot::Temp(2), Slot::Param(1)),
            Instruction::Goto(2),
            Instruction::Label(1),
            Instruction::Pow(Slot::Temp(1), Slot::Param(2), -1, false),
            Instruction::Assign(Slot::Temp(2), Slot::Param(2)),
            Instruction::Label(2),
            Instruction::Join(Slot::Temp(3), Slot::Param(0), Slot::Temp(0), Slot::Temp(1)),
            Instruction::Add(Slot::Out(0), vec![Slot::Temp(3), Slot::Temp(2)], 0),
        ];
        let labels = labels(&instructions).unwrap();
        let live = liveness(&instructions, 1, &labels).unwrap();
        let mut emitter = Emitter {
            instructions: &instructions,
            labels,
            live,
            callees: HashMap::new(),
            next: 0,
            computed: HashMap::new(),
            computed_log: Vec::new(),
            outlines: &mut Outlines::default(),
        };
        let mut environment = (0..3)
            .map(|index| (Slot::Param(index), format!("a{index}")))
            .collect();
        let mut source = String::new();
        emitter
            .range(0, instructions.len(), &mut environment, &mut source, 0)
            .unwrap();
        assert!(source.contains("= if truth(a0) {"));
        assert!(source.contains("} else {"));
        assert!(source.find("if truth").unwrap() < source.find("powi(a1").unwrap());
        assert!(source.find("} else {").unwrap() < source.find("powi(a2").unwrap());
        assert!(
            source.contains("let ("),
            "both Join and the other live slot require phi values"
        );
    }

    #[test]
    fn nested_native_branch_phi_liveness_is_edge_specific() {
        let instructions = vec![
            Instruction::IfElse(Slot::Param(0), 10),
            Instruction::IfElse(Slot::Param(1), 20),
            Instruction::Pow(Slot::Temp(0), Slot::Param(2), -1, false),
            Instruction::Goto(21),
            Instruction::Label(20),
            Instruction::Pow(Slot::Temp(1), Slot::Param(3), -1, false),
            Instruction::Label(21),
            Instruction::Join(Slot::Temp(2), Slot::Param(1), Slot::Temp(0), Slot::Temp(1)),
            Instruction::Goto(11),
            Instruction::Label(10),
            Instruction::Pow(Slot::Temp(3), Slot::Param(4), -1, false),
            Instruction::Label(11),
            Instruction::Join(Slot::Out(0), Slot::Param(0), Slot::Temp(2), Slot::Temp(3)),
        ];
        let labels = labels(&instructions).unwrap();
        let live = liveness(&instructions, 1, &labels).unwrap();
        assert!(!live[8].contains(&Slot::Temp(3)));
        let mut emitter = Emitter {
            instructions: &instructions,
            labels,
            live,
            callees: HashMap::new(),
            next: 0,
            computed: HashMap::new(),
            computed_log: Vec::new(),
            outlines: &mut Outlines::default(),
        };
        let mut environment = (0..5)
            .map(|index| (Slot::Param(index), format!("a{index}")))
            .collect();
        let mut source = String::new();
        emitter
            .range(0, instructions.len(), &mut environment, &mut source, 0)
            .unwrap();
        assert_eq!(source.matches("if truth(").count(), 2);
        assert!(source.find("if truth(a1)").unwrap() < source.find("powi(a2").unwrap());
        assert!(source.rfind("} else {").unwrap() < source.find("powi(a4").unwrap());
    }

    #[test]
    fn scoped_cse_reuses_dominating_values_but_not_sibling_or_escaped_bindings() {
        let call = |target| {
            Instruction::Fun(
                target,
                Box::new((Symbol::LOG, Vec::new(), vec![Slot::Param(1)])),
                false,
            )
        };
        let instructions = vec![
            call(Slot::Temp(10)),
            Instruction::Assign(Slot::Temp(11), Slot::Param(1)),
            Instruction::Pow(Slot::Temp(12), Slot::Temp(11), -1, false),
            Instruction::IfElse(Slot::Param(0), 1),
            call(Slot::Temp(0)),
            Instruction::Pow(Slot::Temp(1), Slot::Param(1), -1, false),
            Instruction::Pow(Slot::Temp(2), Slot::Param(2), -1, false),
            Instruction::Goto(2),
            Instruction::Label(1),
            call(Slot::Temp(3)),
            Instruction::Pow(Slot::Temp(4), Slot::Param(1), -1, false),
            Instruction::Pow(Slot::Temp(5), Slot::Param(2), -1, false),
            Instruction::Label(2),
            Instruction::Join(Slot::Out(0), Slot::Param(0), Slot::Temp(2), Slot::Temp(5)),
            call(Slot::Out(1)),
            Instruction::Pow(Slot::Out(2), Slot::Param(2), -1, false),
        ];
        let labels = labels(&instructions).unwrap();
        let live = liveness(&instructions, 3, &labels).unwrap();
        let mut emitter = Emitter {
            instructions: &instructions,
            labels,
            live,
            callees: [((Symbol::LOG, Vec::new()), "helper_probe".into())].into(),
            next: 0,
            computed: HashMap::new(),
            computed_log: Vec::new(),
            outlines: &mut Outlines::default(),
        };
        let mut environment = (0..3)
            .map(|index| (Slot::Param(index), format!("a{index}")))
            .collect();
        let mut source = String::new();
        emitter
            .range(0, instructions.len(), &mut environment, &mut source, 0)
            .unwrap();
        assert_eq!(source.matches("helper_probe(").count(), 1, "{source}");
        assert_eq!(source.matches("powi(a1,").count(), 1, "{source}");
        assert_eq!(source.matches("powi(a2,").count(), 3, "{source}");
        assert!(source.find("helper_probe").unwrap() < source.find("if truth").unwrap());
        assert!(source.rfind("powi(a2,").unwrap() > source.rfind("} else {").unwrap());
    }

    #[test]
    fn common_branch_alias_overwrites_the_old_live_slot() {
        let instructions = vec![
            Instruction::Assign(Slot::Temp(2), Slot::Param(2)),
            Instruction::IfElse(Slot::Param(0), 1),
            Instruction::Pow(Slot::Temp(0), Slot::Param(1), -1, false),
            Instruction::Assign(Slot::Temp(2), Slot::Param(1)),
            Instruction::Goto(2),
            Instruction::Label(1),
            Instruction::Pow(Slot::Temp(1), Slot::Param(2), -1, false),
            Instruction::Assign(Slot::Temp(2), Slot::Param(1)),
            Instruction::Label(2),
            Instruction::Join(Slot::Out(0), Slot::Param(0), Slot::Temp(0), Slot::Temp(1)),
            Instruction::Assign(Slot::Out(1), Slot::Temp(2)),
        ];
        let labels = labels(&instructions).unwrap();
        let live = liveness(&instructions, 2, &labels).unwrap();
        let mut emitter = Emitter {
            instructions: &instructions,
            labels,
            live,
            callees: HashMap::new(),
            next: 0,
            computed: HashMap::new(),
            computed_log: Vec::new(),
            outlines: &mut Outlines::default(),
        };
        let mut environment = (0..3)
            .map(|index| (Slot::Param(index), format!("a{index}")))
            .collect();
        let mut source = String::new();
        emitter
            .range(0, instructions.len(), &mut environment, &mut source, 0)
            .unwrap();
        assert_eq!(read(&environment, Slot::Temp(2)).unwrap(), "a1");
        assert_eq!(read(&environment, Slot::Out(1)).unwrap(), "a1");
        assert!(
            !source.contains("let ("),
            "common aliases need no phi tuple"
        );
    }

    #[test]
    fn whole_if_keys_keep_capture_identity_aliases_constants_and_output_order() {
        let source = "if truth(&v5) { let v7 = add(a1, a2); (v7, a1) } else { let v8 = mul(a1, ctx.constant(0)); (v8, a2) }";
        let key = conditional_key(source, 2);
        assert_eq!(
            key,
            conditional_key(&source.replace("v7", "v700").replace("v8", "v800"), 2)
        );
        for changed in [
            source.replace("&v5", "&v6"),
            source.replace("a1", "a3"),
            source.replace("add(a1, a2)", "add(a1, a1)"),
            source.replace("constant(0)", "constant(1)"),
            source.replace("(v7, a1)", "(a1, v7)"),
            source.replace("(v8, a2)", "(a2, v8)"),
        ] {
            assert_ne!(key, conditional_key(&changed, 2), "{changed}");
        }
        assert_ne!(key, conditional_key(source, 1));
    }

    #[test]
    fn whole_if_reuse_executes_once_without_conflating_inputs_or_tuple_outputs() {
        fn conditional(
            instructions: &mut Vec<Instruction>,
            id: usize,
            condition: usize,
            yes: usize,
            no: usize,
            output: usize,
            tuple: bool,
            reversed: bool,
        ) {
            let label = id * 2;
            let first = Slot::Temp(id * 3);
            let second = Slot::Temp(id * 3 + 1);
            instructions.push(Instruction::IfElse(Slot::Param(condition), label));
            instructions.push(Instruction::Pow(first, Slot::Param(yes), -1, false));
            if tuple {
                instructions.push(Instruction::Assign(
                    Slot::Temp(99),
                    if reversed { first } else { Slot::Param(yes) },
                ));
            }
            instructions.extend([Instruction::Goto(label + 1), Instruction::Label(label)]);
            instructions.push(Instruction::Pow(second, Slot::Param(no), -1, false));
            if tuple {
                instructions.push(Instruction::Assign(
                    Slot::Temp(99),
                    if reversed { second } else { Slot::Param(no) },
                ));
            }
            instructions.extend([
                Instruction::Label(label + 1),
                Instruction::Join(
                    Slot::Out(output),
                    Slot::Param(condition),
                    if reversed { Slot::Param(yes) } else { first },
                    if reversed { Slot::Param(no) } else { second },
                ),
            ]);
            if tuple {
                instructions.push(Instruction::Assign(Slot::Out(output + 1), Slot::Temp(99)));
            }
        }
        fn body(instructions: Vec<Instruction>, outputs: usize) -> String {
            let graph = Graph {
                input_count: 4,
                output_count: outputs,
                temporary_count: 100,
                instructions,
                constants: Vec::new(),
                constant_functions: Vec::new(),
                sub_evaluators: Vec::new(),
            };
            let mut constants = Constants::default();
            let mut generator = Generator {
                constants: &mut constants,
                helpers: Vec::new(),
                normalized: HashMap::new(),
                outlines: Outlines::default(),
            };
            let body = generator.body(&graph, true, &[]).unwrap();
            assert!(generator.outlines.functions.is_empty());
            body
        }
        let mut instructions = Vec::new();
        conditional(&mut instructions, 1, 0, 1, 2, 0, true, false);
        conditional(&mut instructions, 2, 0, 1, 2, 2, true, false);
        conditional(&mut instructions, 3, 3, 1, 2, 4, true, false);
        conditional(&mut instructions, 4, 0, 2, 1, 6, true, false);
        conditional(&mut instructions, 5, 0, 1, 2, 8, true, true);
        conditional(&mut instructions, 6, 0, 1, 2, 10, false, false);
        conditional(&mut instructions, 7, 0, 1, 1, 11, false, false);
        conditional(&mut instructions, 8, 0, 1, 2, 12, true, false);
        let full = body(instructions, 14);
        assert_eq!(full.matches("if truth").count(), 6, "{full}");
        let mut instructions = Vec::new();
        conditional(&mut instructions, 1, 0, 1, 2, 0, true, false);
        conditional(&mut instructions, 2, 0, 1, 2, 2, true, false);
        let trap = body(instructions, 4);
        assert_eq!(trap.matches("if truth").count(), 1, "{trap}");
        let mut source = PURE_NUMERIC_PROBE.to_owned();
        writeln!(source, "fn evaluate<T: NativeFloat>(input: &[C<T>], ctx: &Context<T>) -> [C<T>;14] {{\n{full}}}").unwrap();
        writeln!(
            source,
            "fn trap<T: NativeFloat>(input: &[C<T>], ctx: &Context<T>) -> [C<T>;4] {{\n{trap}}}"
        )
        .unwrap();
        source.push_str(
            r#"
fn main() {
    use std::sync::atomic::Ordering::Relaxed;
    let ctx = Context(std::marker::PhantomData);
    for (input, expected) in [
        ([1.,2.,4.,0.], [0.5,2.,0.5,2.,0.25,4.,0.25,4.,2.,0.5,0.5,0.5,0.5,2.]),
        ([0.,2.,4.,1.], [0.25,4.,0.25,4.,0.5,2.,0.5,2.,4.,0.25,0.25,0.5,0.25,4.]),
    ] {
        POW_CALLS.store(0, Relaxed);
        assert_eq!(evaluate(&input, &ctx), expected);
        assert_eq!(POW_CALLS.load(Relaxed), 6);
    }
    for (input, expected) in [
        ([1.,2.,0.,0.], [0.5,2.,0.5,2.]),
        ([0.,0.,4.,0.], [0.25,4.,0.25,4.]),
    ] {
        POW_CALLS.store(0, Relaxed);
        assert_eq!(trap(&input, &ctx), expected);
        assert_eq!(POW_CALLS.load(Relaxed), 1);
    }
}
"#,
        );
        compile_probe(&source);
    }

    #[test]
    fn shared_nested_branches_deduplicate_and_keep_selected_tuple_outputs() {
        fn inner(
            instructions: &mut Vec<Instruction>,
            condition: usize,
            yes: Slot,
            no: Slot,
            temp: usize,
            label: usize,
        ) {
            instructions.extend([
                Instruction::IfElse(Slot::Param(condition), label),
                Instruction::Pow(Slot::Temp(temp), yes, -1, false),
                Instruction::Mul(
                    Slot::Temp(temp + 1),
                    vec![Slot::Temp(temp), Slot::Temp(temp)],
                    0,
                ),
                Instruction::Add(Slot::Temp(temp + 2), vec![Slot::Temp(temp + 1), yes], 0),
                Instruction::Assign(Slot::Temp(99), yes),
                Instruction::Goto(label + 1),
                Instruction::Label(label),
                Instruction::Pow(Slot::Temp(temp + 3), no, -1, false),
                Instruction::Mul(
                    Slot::Temp(temp + 4),
                    vec![Slot::Temp(temp + 3), Slot::Temp(temp + 3)],
                    0,
                ),
                Instruction::Add(Slot::Temp(temp + 5), vec![Slot::Temp(temp + 4), no], 0),
                Instruction::Assign(Slot::Temp(99), no),
                Instruction::Label(label + 1),
                Instruction::Join(
                    Slot::Temp(temp + 8),
                    Slot::Param(condition),
                    Slot::Temp(temp + 2),
                    Slot::Temp(temp + 5),
                ),
            ]);
        }
        let mut instructions = vec![
            Instruction::Assign(Slot::Temp(88), Slot::Param(2)),
            Instruction::IfElse(Slot::Param(0), 100),
        ];
        inner(&mut instructions, 1, Slot::Temp(88), Slot::Param(3), 0, 20);
        instructions.extend([Instruction::Goto(101), Instruction::Label(100)]);
        inner(&mut instructions, 4, Slot::Param(5), Slot::Param(6), 10, 30);
        instructions.extend([
            Instruction::Label(101),
            Instruction::Join(Slot::Out(0), Slot::Param(0), Slot::Temp(8), Slot::Temp(18)),
            Instruction::Assign(Slot::Out(1), Slot::Temp(99)),
        ]);
        let graph = Graph {
            input_count: 7,
            output_count: 2,
            temporary_count: 100,
            instructions,
            constants: Vec::new(),
            constant_functions: Vec::new(),
            sub_evaluators: Vec::new(),
        };
        let mut constants = Constants::default();
        let mut generator = Generator {
            constants: &mut constants,
            helpers: Vec::new(),
            normalized: HashMap::new(),
            outlines: Outlines {
                minimum_lines: 3,
                ..Outlines::default()
            },
        };
        // Exercise actual entry parameter aliases, not only an environment
        // assembled by the test: outlined functions cannot refer to `input`.
        let body = generator.body(&graph, true, &[]).unwrap();
        assert_eq!(
            generator.outlines.functions.len(),
            2,
            "four identical leaf arms and two identical nested arms must each deduplicate"
        );
        assert_eq!(body.matches("if truth").count(), 1);
        assert_eq!(body.matches("branch_1(").count(), 2);

        // Compile and execute the emitted Rust itself with a tiny f64 domain.
        // This test never initializes Symbolica. Inactive inverse(0) inputs
        // deliberately panic if outlining accidentally makes either arm eager.
        let mut source = PURE_NUMERIC_PROBE.to_owned();
        writeln!(
            source,
            "fn evaluate<T: NativeFloat>(input: &[C<T>], ctx: &Context<T>) -> [C<T>;2] {{\n{body}}}"
        )
        .unwrap();
        for function in &generator.outlines.functions {
            source.push_str(function);
        }
        source.push_str(
            r#"
fn main() {
    let context = Context(std::marker::PhantomData);
    for (input, expected) in [
        ([1.0,1.0,2.0,0.0,0.0,0.0,0.0], [2.25,2.0]),
        ([1.0,0.0,0.0,-2.0,0.0,0.0,0.0], [-1.75,-2.0]),
        ([0.0,0.0,0.0,0.0,1.0,0.5,0.0], [4.5,0.5]),
        ([0.0,0.0,0.0,0.0,0.0,0.0,4.0], [4.0625,4.0]),
    ] {
        assert_eq!(evaluate(&input,&context), expected);
    }
}
"#,
        );
        compile_probe(&source);
    }

    pub(super) fn compile_probe(source: &str) {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let directory = std::env::temp_dir().join(format!(
            "oneloop-native-outline-{}-{unique}",
            std::process::id()
        ));
        std::fs::create_dir(&directory).unwrap();
        let path = directory.join("probe.rs");
        let binary = directory.join("probe");
        std::fs::write(&path, source).unwrap();
        let compiler = std::env::var_os("RUSTC").unwrap_or_else(|| "rustc".into());
        let compiled = std::process::Command::new(compiler)
            .arg("--edition=2024")
            .arg(&path)
            .arg("-o")
            .arg(&binary)
            .output()
            .unwrap();
        assert!(
            compiled.status.success(),
            "emitted Rust at {} failed: {}",
            path.display(),
            String::from_utf8_lossy(&compiled.stderr)
        );
        let executed = std::process::Command::new(&binary).output().unwrap();
        assert!(
            executed.status.success(),
            "emitted Rust at {} failed: {}",
            path.display(),
            String::from_utf8_lossy(&executed.stderr)
        );
        std::fs::remove_file(binary).unwrap();
        std::fs::remove_file(path).unwrap();
        std::fs::remove_dir(directory).unwrap();
    }
}
