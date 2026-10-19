//! Size and time of building the exact expressions and their evaluators.
//!
//! Callers such as NinjaRS build interpreted evaluators from these expressions
//! for every family (Float and Dual routes, and the tensor bubbles in binary64).
//! Values that enter the finite triangle and box formulas many times (the
//! r-functions, quadratic roots and root lips) are registered helpers, referred
//! to by calls: inlining them once made the shared definitions 5.7 GB and
//! `OneLoopExpressions::new` take 80 s, with 4 s per evaluator; they are now
//! about 60 MB (largest body 19 MB), built in well under a second.
use oneloop::OneLoopExpressions;
use std::time::Instant;
use std::collections::HashMap;
use symbolica::atom::Atom;
use symbolica::evaluate::Instruction;

#[test]
fn definitions_and_evaluators_stay_compact() {
    std::thread::Builder::new()
        .stack_size(512 << 20)
        .spawn(|| {
            let started = Instant::now();
            let expressions = OneLoopExpressions::new();
            let sizes = expressions.definition_sizes();
            let total: usize = sizes.iter().map(|(_, size)| size).sum();
            assert!(
                sizes[0].1 < 32 << 20 && total < 96 << 20,
                "largest definitions {:?}, total {total} bytes",
                &sizes[..sizes.len().min(5)]
            );

            let x: Vec<Atom> = (0..11)
                .map(|i| symbolica::symbol!(format!("construction::x{i}")).to_atom())
                .collect();
            let bubble = |series: oneloop::LaurentSeries| series.coefficients().to_vec();
            let mut families: Vec<(&str, Vec<Atom>, &[Atom])> = vec![
                ("A0", bubble(oneloop::a0(&x[0], &x[1])), &x[..2]),
                (
                    "B0",
                    bubble(oneloop::b0(&x[0], &x[1], &x[2], &x[3])),
                    &x[..4],
                ),
                (
                    "dB0",
                    bubble(oneloop::db0(&x[0], &x[1], &x[2], &x[3])),
                    &x[..4],
                ),
                (
                    "B1",
                    bubble(oneloop::b1(&x[0], &x[1], &x[2], &x[3])),
                    &x[..4],
                ),
                (
                    "B11",
                    bubble(oneloop::b11(&x[0], &x[1], &x[2], &x[3])),
                    &x[..4],
                ),
                (
                    "B111",
                    bubble(oneloop::b111(&x[0], &x[1], &x[2], &x[3])),
                    &x[..4],
                ),
                (
                    "C0",
                    bubble(expressions.c0([&x[0], &x[1], &x[2]], [&x[3], &x[4], &x[5]], &x[6])),
                    &x[..7],
                ),
                (
                    "D0",
                    bubble(expressions.d0(
                        [&x[0], &x[1], &x[2], &x[3], &x[4], &x[5]],
                        [&x[6], &x[7], &x[8], &x[9]],
                        &x[10],
                    )),
                    &x[..11],
                ),
            ];
            let an = oneloop::an(4, &x[0], &x[1])
                .unwrap()
                .into_iter()
                .flat_map(|series| series.coefficients().to_vec())
                .collect();
            families.push(("An", an, &x[..2]));
            for (name, series, parameters) in families {
                let bytes: usize = series.iter().map(|a| a.as_view().get_byte_size()).sum();
                assert!(bytes < 4 << 20, "{name}: {bytes} bytes");
                let built = Instant::now();
                expressions.evaluator(&series, parameters).unwrap();
                let elapsed = built.elapsed();
                assert!(
                    elapsed.as_secs_f64() < 5.0,
                    "{name}: evaluator built in {elapsed:?}"
                );
            }
            let elapsed = started.elapsed();
            assert!(
                elapsed.as_secs_f64() < 30.0,
                "construction took {elapsed:?}"
            );
        })
        .unwrap()
        .join()
        .unwrap();
}

/// Calls of a registered function, repeated with the same arguments where the
/// first call dominates the repeat (the same or an enclosing branch).
fn repeated_calls(instructions: &[Instruction], expensive: &[&str]) -> Vec<String> {
    let mut scopes: Vec<(usize, Option<usize>)> = vec![(0, None)];
    let mut next = 1;
    let mut seen: HashMap<String, Vec<Vec<usize>>> = HashMap::new();
    let mut repeated = Vec::new();
    for instruction in instructions {
        match instruction {
            Instruction::IfElse(..) => {
                scopes.push((next, None));
                next += 1;
            }
            Instruction::Goto(end) => {
                scopes.pop();
                scopes.push((next, Some(*end)));
                next += 1;
            }
            Instruction::Label(label) => {
                if scopes.last().is_some_and(|scope| scope.1 == Some(*label)) {
                    scopes.pop();
                }
            }
            Instruction::Fun(_, function, _) => {
                let (symbol, _, arguments) = &**function;
                let name = symbol.to_string();
                if !expensive.iter().any(|expensive| name.ends_with(expensive)) {
                    continue;
                }
                let path: Vec<usize> = scopes.iter().map(|scope| scope.0).collect();
                let calls = seen.entry(format!("{name} {arguments:?}")).or_default();
                if calls.iter().any(|earlier| path.starts_with(earlier)) {
                    repeated.push(name);
                }
                calls.push(path);
            }
            _ => {}
        }
    }
    repeated
}

/// Symbolica's evaluator calls a registered function anew at every occurrence
/// of its call. The divided difference's cancellation test took each
/// dilogarithm two or three times, and NinjaRS's interpreted D0 and C0 took
/// twice as long at 700 bits; the dilogarithms are now arguments of the
/// selection (`sheet_exact::select_body`). No interpreted D0 or C0 body
/// evaluates a dilogarithm, or a divided difference of them, twice.
#[test]
fn interpreted_dilogarithms_are_evaluated_once() {
    std::thread::Builder::new()
        .stack_size(512 << 20)
        .spawn(|| {
            let expressions = OneLoopExpressions::new();
            let x: Vec<Atom> = (0..11)
                .map(|i| symbolica::symbol!(format!("construction::y{i}")).to_atom())
                .collect();
            let families = [
                (
                    "D0",
                    expressions.d0(
                        [&x[0], &x[1], &x[2], &x[3], &x[4], &x[5]],
                        [&x[6], &x[7], &x[8], &x[9]],
                        &x[10],
                    ),
                    &x[..],
                ),
                (
                    "C0",
                    expressions.c0([&x[0], &x[1], &x[2]], [&x[3], &x[4], &x[5]], &x[6]),
                    &x[..7],
                ),
            ];
            let expensive = [
                "__olo_sheet_dilog",
                "__olo_sheet_dilog_complement",
                "__olo_sheet_dilog_difference",
                "__olo_sheet_dilog_difference_complement",
                "__olo_sheet_dilog_difference_select",
            ];
            let mut failures = Vec::new();
            for (family, series, parameters) in families {
                let exported = expressions
                    .evaluator(series.coefficients(), parameters)
                    .unwrap()
                    .export_instructions();
                // Each body lists the bodies it calls.
                let mut pending = vec![("top".to_string(), &exported)];
                let mut visited = std::collections::HashSet::new();
                while let Some((body, program)) = pending.pop() {
                    let repeated = repeated_calls(&program.instructions, &expensive);
                    if !repeated.is_empty() {
                        failures.push(format!("{family} {body}: {repeated:?}"));
                    }
                    for sub in &program.sub_evaluators {
                        if visited.insert(sub.symbol.to_string()) {
                            pending.push((sub.symbol.to_string(), &*sub.instructions));
                        }
                    }
                }
            }
            assert!(failures.is_empty(), "{failures:#?}");
        })
        .unwrap()
        .join()
        .unwrap();
}
