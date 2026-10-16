//! Branch-free lowering of the existing exported instruction stream.
//!
//! Function metadata below describes storage and emitted Rust, not a second IR.
//! Each canonical shared body gets disjoint slots in one reusable arena. Calls
//! copy argument columns, execute the helper once for all lanes, then copy its
//! result column. This preserves shared-body boundaries without per-row calls.

use super::*;

pub(super) fn check_control_flow(graph: &Graph) -> Result<(), String> {
    for instruction in &graph.instructions {
        if matches!(
            instruction,
            Instruction::IfElse(..)
                | Instruction::Join(..)
                | Instruction::Goto(..)
                | Instruction::Label(..)
        ) {
            return Err(
                "native batch emission requires branch-free instructions (including shared bodies)"
                    .into(),
            );
        }
    }
    for child in &graph.sub_evaluators {
        check_control_flow(&child.instructions)?;
    }
    Ok(())
}

#[derive(Clone)]
struct Function {
    name: String,
    base: usize,
    input_count: usize,
    output_offset: usize,
}

struct BatchGenerator<'a> {
    constants: &'a mut Constants,
    slots: usize,
    functions: Vec<String>,
    normalized: HashMap<String, Function>,
}

pub(super) fn emit(graph: &Graph, constants: &mut Constants) -> Result<String, String> {
    let mut generator = BatchGenerator {
        constants,
        slots: 0,
        functions: Vec::new(),
        normalized: HashMap::new(),
    };
    let root = generator.function(graph, false)?;
    let mut source = String::from(
        "\nuse super::BatchWorkspace;\nuse super::batch::{copy_column, copy_value};\n\n",
    );
    // All dimensions are checked before packing, arithmetic or output writes.
    writeln!(source, "pub(super) fn evaluate_batch<T: NativeFloat>(input: &[C<T>], output: &mut [C<T>], rows: usize, ctx: &Context<T>, workspace: &mut BatchWorkspace<T>) -> Result<(), String> {{
    let inputs = rows.checked_mul({}).ok_or(\"batch input size overflow\")?;
    let outputs = rows.checked_mul({}).ok_or(\"batch output size overflow\")?;
    if input.len() != inputs || output.len() != outputs {{ return Err(\"batch buffer shape mismatch\".into()); }}
    if rows == 0 {{ return Ok(()); }}
    if rows == 1 {{
        let values = evaluate(input, ctx);
        for (target, source) in output.iter_mut().zip(&values) {{ copy_value(target, source); }}
        return Ok(());
    }}
    let (v, stride) = workspace.prepare(rows, {})?;", graph.input_count, graph.output_count, generator.slots).unwrap();
    for index in 0..graph.input_count {
        writeln!(source, "    for lane in 0..rows {{ copy_value(&mut v[{} * stride + lane], &input[lane * {} + {index}]); }}", root.base + index, graph.input_count).unwrap();
    }
    writeln!(
        source,
        "    {}(v, stride, rows, {}, ctx);",
        root.name, root.base
    )
    .unwrap();
    for index in 0..graph.output_count {
        writeln!(source, "    for lane in 0..rows {{ copy_value(&mut output[lane * {} + {index}], &v[{} * stride + lane]); }}", graph.output_count, root.base + root.output_offset + index).unwrap();
    }
    source.push_str("    Ok(())\n}\n\n");
    for function in generator.functions {
        source.push_str(&function);
    }
    Ok(source)
}

impl BatchGenerator<'_> {
    fn function(&mut self, graph: &Graph, shared: bool) -> Result<Function, String> {
        if shared && graph.output_count != 1 {
            return Err("native batch helper must have one output".into());
        }
        let mut callees = HashMap::new();
        for child in &graph.sub_evaluators {
            let function = self.function(&child.instructions, true)?;
            if let Some(previous) =
                callees.insert((child.symbol, child.tags.clone()), function.clone())
                && previous.name != function.name
            {
                return Err("ambiguous exported helper signature".into());
            }
        }
        let constants = self.constants.graph_constants(graph)?;
        let output_offset = graph
            .input_count
            .checked_add(graph.temporary_count)
            .ok_or("batch slot count overflow")?;
        let slots = output_offset
            .checked_add(graph.output_count)
            .ok_or("batch slot count overflow")?;
        let index = |slot| -> Result<usize, String> {
            match slot {
                Slot::Param(index) if index < graph.input_count => Ok(index),
                Slot::Temp(index) if index < graph.temporary_count => Ok(graph.input_count + index),
                Slot::Out(index) if index < graph.output_count => Ok(output_offset + index),
                _ => Err(format!("invalid batch storage slot {slot:?}")),
            }
        };
        let read = |slot: Slot| -> Result<String, String> {
            match slot {
                Slot::Const(i) => constants
                    .get(i)
                    .map(|i| format!("ctx.constant({i})"))
                    .ok_or_else(|| format!("invalid batch constant slot {i}")),
                _ => Ok(format!("&v[(base + {}) * stride + lane]", index(slot)?)),
            }
        };
        let mut defined: BTreeSet<_> = (0..graph.input_count)
            .map(Slot::Param)
            .chain((0..constants.len()).map(Slot::Const))
            .collect();
        let mut body = String::new();
        for instruction in &graph.instructions {
            let (target, arguments) = uses(instruction);
            for argument in arguments {
                if !defined.contains(&argument) {
                    return Err(format!("batch slot read before definition: {argument:?}"));
                }
                read(argument)?;
            }
            let target = target.ok_or("unsupported control flow in native batch emission")?;
            let destination = index(target)?;
            let mut assign = |source: Slot| -> Result<(), String> {
                if matches!(source, Slot::Const(_)) {
                    writeln!(body, "    for lane in 0..rows {{ copy_value(&mut v[(base + {destination}) * stride + lane], {}); }}", read(source)?).unwrap();
                } else {
                    writeln!(
                        body,
                        "    copy_column(v, stride, rows, base + {destination}, base + {});",
                        index(source)?
                    )
                    .unwrap();
                }
                Ok(())
            };
            let expression = match instruction {
                Instruction::Assign(_, source) => {
                    assign(*source)?;
                    None
                }
                Instruction::Add(_, args, _) | Instruction::Mul(_, args, _) => {
                    let operation = if matches!(instruction, Instruction::Add(..)) {
                        "add"
                    } else {
                        "mul"
                    };
                    let first = *args.first().ok_or("empty arithmetic instruction")?;
                    if args.len() == 1 {
                        assign(first)?;
                        None
                    } else {
                        // Preserve exported operand order. Finish all reads before
                        // replacing the target, even when Symbolica reuses its slot.
                        let mut value = read(first)?;
                        for arg in &args[1..] {
                            value = format!("&{operation}({value}, {})", read(*arg)?);
                        }
                        Some(value[1..].to_owned())
                    }
                }
                Instruction::Pow(_, base, exponent, real) => {
                    Some(format!("powi({}, {exponent}, {real})", read(*base)?))
                }
                Instruction::Powf(_, base, exponent, real) => Some(format!(
                    "powf({}, {}, {real})",
                    read(*base)?,
                    read(*exponent)?
                )),
                Instruction::Fun(_, function, _) => {
                    let (symbol, tags, arguments) = &**function;
                    if let Some(callee) = callees.get(&(*symbol, tags.clone())) {
                        if arguments.len() != callee.input_count {
                            return Err("native batch helper argument count mismatch".into());
                        }
                        for (argument, source) in arguments.iter().enumerate() {
                            let target = callee.base + argument;
                            if matches!(source, Slot::Const(_)) {
                                writeln!(body, "    for lane in 0..rows {{ copy_value(&mut v[{target} * stride + lane], {}); }}", read(*source)?).unwrap();
                            } else {
                                writeln!(
                                    body,
                                    "    copy_column(v, stride, rows, {target}, base + {});",
                                    index(*source)?
                                )
                                .unwrap();
                            }
                        }
                        writeln!(
                            body,
                            "    {}(v, stride, rows, {}, ctx);",
                            callee.name, callee.base
                        )
                        .unwrap();
                        writeln!(
                            body,
                            "    copy_column(v, stride, rows, base + {destination}, {});",
                            callee.base + callee.output_offset
                        )
                        .unwrap();
                        None
                    } else {
                        let name = primitive_name(*symbol, tags)?;
                        if arguments.len() != 1 {
                            return Err(format!("unexpected arity for {name}"));
                        }
                        Some(format!("{name}({})", read(arguments[0])?))
                    }
                }
                _ => return Err("unsupported control flow in native batch emission".into()),
            };
            if let Some(expression) = expression {
                writeln!(body, "    for lane in 0..rows {{ v[(base + {destination}) * stride + lane] = {expression}; }}").unwrap();
            }
            defined.insert(target);
        }
        for output in 0..graph.output_count {
            if !defined.contains(&Slot::Out(output)) {
                return Err(format!("undefined batch output {output}"));
            }
        }
        // Just as in scalar emission, canonical child names and global constant
        // indices make the emitted body a safe implementation identity.
        let key = format!(
            "{}:{}:{}:{body}",
            graph.input_count, graph.temporary_count, graph.output_count
        );
        if let Some(function) = self.normalized.get(&key) {
            return Ok(function.clone());
        }
        let function = Function {
            name: format!("batch_body_{}", self.functions.len()),
            base: self.slots,
            input_count: graph.input_count,
            output_offset,
        };
        self.slots = self
            .slots
            .checked_add(slots)
            .ok_or("batch module slot count overflow")?;
        self.functions.push(format!("fn {}<T: NativeFloat>(v: &mut [C<T>], stride: usize, rows: usize, base: usize, ctx: &Context<T>) {{\n{body}}}\n\n", function.name));
        self.normalized.insert(key, function.clone());
        Ok(function)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(instructions: Vec<Instruction>) -> Graph {
        Graph {
            input_count: 1,
            output_count: 1,
            temporary_count: 1,
            instructions,
            constants: vec![Complex::new(Rational::from(7), Rational::from(0))],
            constant_functions: Vec::new(),
            sub_evaluators: Vec::new(),
        }
    }

    #[test]
    fn rejection_checks_all_nested_control_flow_and_is_transactional() {
        let mut emitter = RustEmitter::default();
        let valid = graph(vec![Instruction::Assign(Slot::Out(0), Slot::Const(0))]);
        emitter
            .emit_batch_graph(&valid, EmitOptions::default())
            .unwrap();
        let before = emitter.constants_source();
        for instruction in [
            Instruction::IfElse(Slot::Param(0), 1),
            Instruction::Join(Slot::Out(0), Slot::Param(0), Slot::Const(0), Slot::Const(0)),
            Instruction::Goto(1),
            Instruction::Label(1),
        ] {
            let mut bad = graph(vec![instruction]);
            bad.constants[0] = Complex::new(Rational::from(11), Rational::from(0));
            for depth in 0..3 {
                assert!(
                    emitter
                        .emit_batch_graph(&bad, EmitOptions::default())
                        .unwrap_err()
                        .contains("branch-free")
                );
                assert_eq!(emitter.constants_source(), before);
                let mut parent = valid.clone();
                parent.sub_evaluators.push(SubEvaluator {
                    symbol: symbolica::symbol!(format!("batch_rejection_{depth}")),
                    tags: vec![],
                    instructions: bad,
                });
                bad = parent;
            }
        }
        for instructions in [
            vec![Instruction::Assign(Slot::Out(0), Slot::Temp(0))],
            vec![
                Instruction::Assign(Slot::Temp(2), Slot::Const(0)),
                Instruction::Assign(Slot::Out(0), Slot::Temp(2)),
            ],
            vec![Instruction::Fun(
                Slot::Out(0),
                Box::new((
                    symbolica::symbol!("batch_unknown"),
                    vec![],
                    vec![Slot::Param(0)],
                )),
                false,
            )],
            vec![Instruction::Add(Slot::Out(0), vec![], 0)],
        ] {
            let mut bad = graph(instructions);
            bad.constants[0] = Complex::new(Rational::from(13), Rational::from(0));
            assert!(
                emitter
                    .emit_batch_graph(&bad, EmitOptions::default())
                    .is_err()
            );
            assert_eq!(emitter.constants_source(), before);
        }
    }

    #[test]
    fn reused_slots_and_aliases_keep_operand_order_and_bound_storage() {
        let graph = Graph {
            input_count: 3,
            output_count: 2,
            temporary_count: 1,
            instructions: vec![
                Instruction::Assign(Slot::Temp(0), Slot::Param(0)),
                Instruction::Add(
                    Slot::Temp(0),
                    vec![Slot::Temp(0), Slot::Param(1), Slot::Param(2)],
                    0,
                ),
                Instruction::Assign(Slot::Temp(0), Slot::Temp(0)),
                Instruction::Assign(Slot::Out(0), Slot::Temp(0)),
                Instruction::Mul(Slot::Out(1), vec![Slot::Param(2)], 0),
            ],
            constants: vec![],
            constant_functions: vec![],
            sub_evaluators: vec![],
        };
        let module = RustEmitter::default()
            .emit_batch_graph(&graph, EmitOptions::default())
            .unwrap();
        assert!(module.source.contains("workspace.prepare(rows, 6)"));
        let source = format!(
            "{}\n{}\nmod primitives {{ pub(super) use super::{{add,mul,powi,truth}}; }}\nmod generated {{\n{}\n}}\nfn main() {{ let ctx = Context(std::marker::PhantomData); let mut workspace = BatchWorkspace {{ values: vec![] }}; let input = [1e16,-1e16,1.0,2.0,3.0,4.0]; let mut output = [99.0;4]; generated::evaluate_batch(&input,&mut output,2,&ctx,&mut workspace).unwrap(); assert_eq!(output,[1.0,1.0,9.0,4.0]); assert_eq!(generated::evaluate(&input[..3],&ctx),[1.0,1.0]); }}",
            super::super::tests::PURE_NUMERIC_PROBE,
            r#"struct BatchWorkspace<T> { values:Vec<T> }
impl<T:NativeFloat> BatchWorkspace<T> { fn prepare(&mut self,rows:usize,slots:usize)->Result<(&mut[T],usize),String> { self.values.resize(rows*slots,T::from(0.0));Ok((&mut self.values,rows)) } }
mod batch {
 pub fn copy_value<T:Copy>(target:&mut T,source:&T) { *target=*source; }
 pub fn copy_column<T:Copy>(v:&mut[T],stride:usize,rows:usize,target:usize,source:usize) { for lane in 0..rows { v[target*stride+lane]=v[source*stride+lane]; } }
}"#,
            module.source
        );
        super::super::tests::compile_probe(&source);

        assert!(module.source.contains("add(&add(&v[(base + 3) * stride + lane], &v[(base + 1) * stride + lane]), &v[(base + 2) * stride + lane])"));
    }
}
