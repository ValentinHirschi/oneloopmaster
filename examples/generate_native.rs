//! Generate generic Rust scalar statements from the repaired exact expressions.
//! Run: cargo run --features rust-codegen --example generate_native -- src/native
use oneloop::rust_codegen::{EmitOptions, RustEmitter};
use oneloop::{A0, B0, C0, D0, OneLoopExpressions, ScalarIntegral, dB0};
use std::path::PathBuf;
use symbolica::atom::Atom;

fn generate(directory: PathBuf) -> Result<(), String> {
    let context = OneLoopExpressions::new();
    let mut emitter = RustEmitter::default();
    let mut files = Vec::new();
    for (family, symbol, file) in [
        (ScalarIntegral::A0, A0(), "generated_a0.rs"),
        (ScalarIntegral::B0, B0(), "generated_b0.rs"),
        (ScalarIntegral::DB0, dB0(), "generated_db0.rs"),
        (ScalarIntegral::C0, C0(), "generated_c0.rs"),
        (ScalarIntegral::D0, D0(), "generated_d0.rs"),
    ] {
        eprintln!("exporting native {}", family.name());
        let parameters = (0..family.arity())
            .map(|index| {
                symbolica::symbol!(format!("native_codegen_{}_{index}", family.name())).to_atom()
            })
            .collect::<Vec<_>>();
        let calls = [0, -1, -2].map(|tag| {
            let mut arguments = vec![Atom::num(tag)];
            arguments.extend(parameters.iter().cloned());
            symbol.call(&arguments)
        });
        // Sharing compensated arithmetic pays off in the large triangle/box
        // graphs, but call overhead is significant for the small families.
        let inline_arithmetic = matches!(
            family,
            ScalarIntegral::A0 | ScalarIntegral::B0 | ScalarIntegral::DB0
        );
        let module = emitter.emit(
            context
                .evaluator(&calls, &parameters)
                .map_err(|error| error.to_string())?,
            EmitOptions { inline_arithmetic },
        )?;
        eprintln!(
            "{}: {} instructions, {} canonical helpers, {} outlined branches, {} bytes",
            family.name(),
            module.instruction_count,
            module.helper_count,
            module.outlined_branch_count,
            module.source.len()
        );
        files.push((file, module.source));
    }
    files.push(("generated_constants.rs", emitter.constants_source()));
    // Construct every file before replacing any generated module. Unknown
    // instructions or unresolved constants cannot leave a partial new set.
    std::fs::create_dir_all(&directory).map_err(|error| error.to_string())?;
    for (file, source) in files {
        std::fs::write(directory.join(file), source).map_err(|error| error.to_string())?;
    }
    eprintln!("{} exact cached constants", emitter.constants().len());
    Ok(())
}

fn main() {
    let mut arguments = std::env::args().skip(1);
    let directory = PathBuf::from(arguments.next().unwrap_or_else(|| "src/native".into()));
    assert!(
        arguments.next().is_none(),
        "usage: generate_native [output-directory]"
    );
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(move || generate(directory))
        .unwrap()
        .join()
        .unwrap()
        .unwrap();
}
