#![cfg(feature = "rust-codegen")]
//! Compile and execute the generated Rust against the real native runtime.
//! This also compares Symbolica's independent batch interpreter and records the
//! primitive visitation order; scalar equality alone cannot prove batching.
use oneloop::rust_codegen::{EmitOptions, RustEmitter};
use std::{path::Path, process::Command, time::Instant};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::{dual::HyperDual, float::Complex, rational::Rational},
    evaluate::{Dualizer, FunctionRegistrationOptions, InliningPolicy},
};

fn exact() -> symbolica::evaluate::ExpressionEvaluator<Complex<Rational>> {
    let x = symbolica::parse!("batch_test::x");
    let y = symbolica::parse!("batch_test::y");
    Atom::evaluator_multiple(&[
        symbolica::parse!("batch_test::f(batch_test::x,batch_test::y)+batch_test::f(batch_test::y,batch_test::x)"),
        symbolica::parse!("batch_test::x+batch_test::y+100000000000000000001/100000000000000000000"),
        x.clone(),
    ], &[x, y])
    .add_function_with_options(symbolica::symbol!("batch_test::g"),vec![symbolica::symbol!("batch_test::t")],symbolica::parse!("batch_test::t^2"),FunctionRegistrationOptions::new().inlining(InliningPolicy::Never)).unwrap()
    .add_function_with_options(symbolica::symbol!("batch_test::f"), vec![symbolica::symbol!("batch_test::a"),symbolica::symbol!("batch_test::b")], symbolica::parse!("(batch_test::g(batch_test::a)+batch_test::b)/(1+batch_test::b)"), FunctionRegistrationOptions::new().inlining(InliningPolicy::Never)).unwrap()
    .cores(1).build().unwrap()
}

fn dual() -> symbolica::evaluate::ExpressionEvaluator<Complex<Rational>> {
    exact()
        .vectorize(&Dualizer::new(
            HyperDual::<Complex<Rational>>::new(vec![
                vec![0, 0],
                vec![1, 0],
                vec![0, 1],
                vec![2, 0],
                vec![1, 1],
                vec![0, 2],
            ]),
            vec![],
        ))
        .unwrap()
}

fn library(directory: &Path, name: &str) -> std::path::PathBuf {
    // Resolve the exact artifacts linked into this test, not the newest rlib:
    // one target directory can contain different Symbolica feature/revision
    // builds. Cargo records the dependency fingerprint in the test fingerprint.
    let executable = std::env::current_exe().unwrap();
    let test_hash = executable
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .rsplit('-')
        .next()
        .unwrap();
    let fingerprints = directory.parent().unwrap().join(".fingerprint");
    let test = std::fs::read_to_string(fingerprints.join(format!(
        "oneloop-{test_hash}/test-integration-test-rust_batch_codegen.json"
    )))
    .unwrap();
    let marker = format!("\"{name}\",false,");
    let number: String = test
        .split_once(&marker)
        .unwrap()
        .1
        .chars()
        .take_while(char::is_ascii_digit)
        .collect();
    let fingerprint = format!("{:016x}", number.parse::<u64>().unwrap().swap_bytes());
    for entry in std::fs::read_dir(&fingerprints)
        .unwrap()
        .filter_map(Result::ok)
    {
        let entry_name = entry.file_name();
        let entry_name = entry_name.to_string_lossy();
        if let Some(hash) = entry_name.strip_prefix(&format!("{name}-"))
            && std::fs::read_to_string(entry.path().join(format!("lib-{name}")))
                .is_ok_and(|text| text.trim() == fingerprint)
        {
            let artifact = directory.join(format!("lib{name}-{hash}.rlib"));
            if artifact.exists() {
                return artifact;
            }
        }
    }
    panic!("cannot find Cargo's {name} dependency artifact {fingerprint}");
}

#[test]
fn generated_batch_numeric_precision_dual_trace_and_allocation_contracts() {
    let _ = symbolica::symbol!("batch_codegen_state");
    let started = Instant::now();
    let mut emitter = RustEmitter::default();
    let ordinary = emitter.emit_batch(exact(), EmitOptions::default()).unwrap();
    assert!(ordinary.helper_count > 0);
    let scalar_bytes = RustEmitter::default()
        .emit(exact(), EmitOptions::default())
        .unwrap()
        .source
        .len();
    eprintln!(
        "ordinary source: scalar_bytes={scalar_bytes} scalar_plus_batch_bytes={}",
        ordinary.source.len()
    );

    let dual = emitter
        .emit_batch(
            dual(),
            EmitOptions {
                inline_arithmetic: true,
            },
        )
        .unwrap();
    assert_eq!(dual.input_count, 12);
    assert_eq!(dual.output_count, 18);
    let identity = emitter
        .emit_batch(
            Atom::evaluator_multiple(
                &[symbolica::parse!("x"), symbolica::parse!("y")],
                &[symbolica::parse!("x"), symbolica::parse!("y")],
            )
            .build()
            .unwrap(),
            EmitOptions::default(),
        )
        .unwrap();
    let constants = emitter
        .emit_batch(
            Atom::evaluator_multiple(
                &[
                    symbolica::parse!("100000000000000000001/100000000000000000000"),
                    symbolica::parse!("pi"),
                ],
                &[] as &[Atom],
            )
            .build()
            .unwrap(),
            EmitOptions::default(),
        )
        .unwrap();
    let traced = emitter
        .emit_batch(
            symbolica::parse!("batch_trace::f(x)+abs(x)")
                .evaluator(&[symbolica::parse!("x")])
                .add_function_with_options(
                    symbolica::symbol!("batch_trace::f"),
                    vec![symbolica::symbol!("z")],
                    symbolica::parse!("log(z)+sqrt(z)"),
                    FunctionRegistrationOptions::new().inlining(InliningPolicy::Never),
                )
                .unwrap()
                .cores(1)
                .build()
                .unwrap(),
            EmitOptions::default(),
        )
        .unwrap();
    let empty = emitter
        .emit_batch(
            Atom::evaluator_multiple(&[] as &[Atom], &[] as &[Atom])
                .build()
                .unwrap(),
            EmitOptions::default(),
        )
        .unwrap();
    let mut source = String::from(PROBE);
    for (name, module) in [
        ("ordinary", ordinary),
        ("derivatives", dual),
        ("identity", identity),
        ("constants", constants),
        ("traced", traced),
        ("empty", empty),
    ] {
        source.push_str(&format!("\nmod {name} {{\n{}\n}}\n", module.source));
    }
    source.push_str(&format!(
        "\nmod table {{\n{}\n}}\n",
        emitter.constants_source()
    ));
    let export_time = started.elapsed();
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "oneloop-native-batch-{}-{unique}",
        std::process::id()
    ));
    std::fs::create_dir(&directory).unwrap();
    let path = directory.join("probe.rs");
    let binary = directory.join("probe");
    std::fs::write(&path, &source).unwrap();
    let executable = std::env::current_exe().unwrap();
    let dependencies = executable.parent().unwrap();
    let started = Instant::now();
    let compiled = Command::new(std::env::var_os("RUSTC").unwrap_or_else(|| "rustc".into()))
        .env("CARGO_CRATE_NAME", "native_batch_probe")
        .arg("--edition=2024")
        .arg("-C")
        .arg("opt-level=2")
        .arg(&path)
        .arg("-o")
        .arg(&binary)
        .arg("-L")
        .arg(format!("dependency={}", dependencies.display()))
        .arg("--extern")
        .arg(format!(
            "oneloop={}",
            library(dependencies, "oneloop").display()
        ))
        .arg("--extern")
        .arg(format!(
            "symbolica={}",
            library(dependencies, "symbolica").display()
        ))
        .output()
        .unwrap();
    let compile_time = started.elapsed();
    assert!(
        compiled.status.success(),
        "generated source {} failed: {}",
        path.display(),
        String::from_utf8_lossy(&compiled.stderr)
    );
    let executed = Command::new(&binary).output().unwrap();
    assert!(
        executed.status.success(),
        "generated source {} failed: {}\n{}",
        path.display(),
        String::from_utf8_lossy(&executed.stdout),
        String::from_utf8_lossy(&executed.stderr)
    );
    eprintln!(
        "native batch probe: source_bytes={} prepare_and_emit_seconds={:.6} optimized_rustc_seconds={:.6}; {}",
        source.len(),
        export_time.as_secs_f64(),
        compile_time.as_secs_f64(),
        String::from_utf8_lossy(&executed.stdout)
    );
    std::fs::remove_dir_all(directory).unwrap();
}

const PROBE: &str = r#"
#![allow(unused_imports, dead_code)]
use oneloop::rust_codegen::runtime::{C, Context, NativeFloat as RuntimeFloat, ConstantSpec, BatchWorkspace, batch};
trait NativeFloat: RuntimeFloat + symbolica::domains::float::RealLike {}
impl<T:RuntimeFloat + symbolica::domains::float::RealLike> NativeFloat for T {}
use symbolica::{atom::{Atom, AtomCore}, domains::{dual::HyperDual, float::{Complex, DoubleFloat, Float, FloatLike, Real, RealLike}, rational::Rational}, evaluate::{Dualizer, FunctionRegistrationOptions, InliningPolicy}};
use std::sync::Mutex;
static TRACE: Mutex<Vec<(&str,f64)>> = Mutex::new(Vec::new());
mod primitives {
    pub use oneloop::rust_codegen::runtime::primitives::*;
    pub fn log<T: super::NativeFloat>(a: &C<T>) -> C<T> { super::TRACE.lock().unwrap().push(("log", symbolica::domains::float::RealLike::to_f64(&a.re))); oneloop::rust_codegen::runtime::primitives::log(a) }
    pub fn sqrt<T: super::NativeFloat>(a: &C<T>) -> C<T> { super::TRACE.lock().unwrap().push(("sqrt", symbolica::domains::float::RealLike::to_f64(&a.re))); oneloop::rust_codegen::runtime::primitives::sqrt(a) }
    pub fn abs<T: super::NativeFloat>(a: &C<T>) -> C<T> { super::TRACE.lock().unwrap().push(("abs", symbolica::domains::float::RealLike::to_f64(&a.re))); oneloop::rust_codegen::runtime::primitives::abs(a) }
    pub fn powi<T: super::NativeFloat>(a: &C<T>, n:i64, real:bool) -> C<T> { super::TRACE.lock().unwrap().push(("powi", symbolica::domains::float::RealLike::to_f64(&a.re))); oneloop::rust_codegen::runtime::primitives::powi(a,n,real) }
}
// Allocation measurement is restricted to the current thread and excludes
// setup, oracle evaluation, output allocation, and trace growth.
struct Allocator;
thread_local! { static COUNT: std::cell::Cell<Option<usize>> = const { std::cell::Cell::new(None) }; }
unsafe impl std::alloc::GlobalAlloc for Allocator {
    unsafe fn alloc(&self, layout: std::alloc::Layout) -> *mut u8 {
        let _ = COUNT.try_with(|count| if let Some(n)=count.get() { count.set(Some(n+1)); });
        unsafe { std::alloc::System.alloc(layout) }
    }
    unsafe fn dealloc(&self, ptr:*mut u8, layout:std::alloc::Layout) { unsafe { std::alloc::System.dealloc(ptr,layout) } }
    unsafe fn realloc(&self, ptr:*mut u8, layout:std::alloc::Layout, size:usize)->*mut u8 {
        let _ = COUNT.try_with(|count| if let Some(n)=count.get() { count.set(Some(n+1)); });
        unsafe { std::alloc::System.realloc(ptr,layout,size) }
    }
}
#[global_allocator] static ALLOCATOR: Allocator = Allocator;
fn exact() -> symbolica::evaluate::ExpressionEvaluator<Complex<Rational>> {
    let x=symbolica::parse!("batch_test::x"); let y=symbolica::parse!("batch_test::y");
    Atom::evaluator_multiple(&[
        symbolica::parse!("batch_test::f(batch_test::x,batch_test::y)+batch_test::f(batch_test::y,batch_test::x)"),
        symbolica::parse!("batch_test::x+batch_test::y+100000000000000000001/100000000000000000000"), x.clone()], &[x,y])
    .add_function_with_options(symbolica::symbol!("batch_test::g"),vec![symbolica::symbol!("batch_test::t")],symbolica::parse!("batch_test::t^2"),FunctionRegistrationOptions::new().inlining(InliningPolicy::Never)).unwrap()
    .add_function_with_options(symbolica::symbol!("batch_test::f"),vec![symbolica::symbol!("batch_test::a"),symbolica::symbol!("batch_test::b")],symbolica::parse!("(batch_test::g(batch_test::a)+batch_test::b)/(1+batch_test::b)"),FunctionRegistrationOptions::new().inlining(InliningPolicy::Never)).unwrap().cores(1).build().unwrap()
}
fn check<T: NativeFloat>(bits:u32) where C<T>: symbolica::evaluate::EvaluationDomain {
    let prototype=T::zero_at(bits);
    let n=|value:i64| C::new(prototype.from_i64(value),prototype.zero());
    let ctx=Context::<T>::new(table::CONSTANTS,bits).unwrap();
    let mut workspace=BatchWorkspace::new(prototype.clone());
    let mut oracle=exact().map_coeff_with_prec(&|c| C::new(prototype.from_rational(&c.re),prototype.from_rational(&c.im)),bits).into_batch_evaluator();
    for rows in [0,1,8,3,32,128,1024,17,2,1,0] {
        let input=(0..rows).flat_map(|i| [n((rows-i) as i64+1),C::new(prototype.from_i64(i as i64+2),prototype.from_i64(1))]).collect::<Vec<_>>();
        let mut output=vec![n(99);rows*3]; let mut interpreted=output.clone();
        ordinary::evaluate_batch(&input,&mut output,rows,&ctx,&mut workspace).unwrap();
        oracle.evaluate_batch(&input,&mut interpreted,rows).unwrap();
        for row in 0..rows {
            let scalar=ordinary::evaluate(&input[row*2..row*2+2],&ctx);
            for j in 0..3 {
                assert_eq!(output[row*3+j],scalar[j],"batch/scalar mismatch");
                let error=(output[row*3+j].clone()-interpreted[row*3+j].clone()).norm().re.to_f64();
                let scale=1.0+interpreted[row*3+j].norm().re.to_f64();
                assert!(error/scale < 2.0_f64.powi(-(bits as i32)+12),"interpreter disagreement {error}");
            }
        }
    }
    assert_eq!(workspace.capacity(),1024);
    let mut output=vec![n(99);6]; let before=output.clone();
    for (input,rows) in [(vec![n(1)],2),(vec![],usize::MAX),(vec![],0)] {
        assert!(ordinary::evaluate_batch(&input,&mut output,rows,&ctx,&mut workspace).is_err()); assert_eq!(output,before);
    }
    let mut solo=BatchWorkspace::new(prototype.clone());
    ordinary::evaluate_batch(&[n(2),n(3)],&mut output[..3],1,&ctx,&mut solo).unwrap(); assert_eq!(solo.storage_len(),0);
    let mut values=vec![n(0);8]; constants::evaluate_batch(&[],&mut values,4,&ctx,&mut workspace).unwrap();
    assert_eq!(values[..2],constants::evaluate(&[],&ctx));
    if bits>53 { assert_ne!(values[0].re,prototype.from_i64(1)); }
    // Flattened higher/mixed Taylor components preserve Symbolica's exact order.
    let shape=vec![vec![0,0],vec![1,0],vec![0,1],vec![2,0],vec![1,1],vec![0,2]];
    let mut oracle=exact().vectorize(&Dualizer::new(HyperDual::<Complex<Rational>>::new(shape),vec![])).unwrap().map_coeff_with_prec(&|c| C::new(prototype.from_rational(&c.re),prototype.from_rational(&c.im)),bits).into_batch_evaluator();
    for rows in [1,8,32,3] {
        let input=(0..rows).flat_map(|i| [n(i as i64+1),n(1),n(0),n(0),n(0),n(0),n(i as i64+2),n(0),n(1),n(0),n(0),n(0)]).collect::<Vec<_>>();
        let mut output=vec![n(99);rows*18];let mut expected=output.clone();
        derivatives::evaluate_batch(&input,&mut output,rows,&ctx,&mut workspace).unwrap(); oracle.evaluate_batch(&input,&mut expected,rows).unwrap();
        for row in 0..rows { let scalar=derivatives::evaluate(&input[row*12..row*12+12],&ctx); assert_eq!(&output[row*18..row*18+18],&scalar); }
        for (actual,expected) in output.iter().zip(&expected) { let error=(actual.clone()-expected.clone()).norm().re.to_f64(); assert!(error/(1.0+expected.norm().re.to_f64()) < 2.0_f64.powi(-(bits as i32)+14)); }
        assert!(output[4].re != prototype.zero(),"mixed coefficient unexpectedly absent");
    }
}
fn main() {
    let _ = symbolica::symbol!("batch_codegen_state");
    check::<f64>(53); check::<DoubleFloat>(106);check::<Float>(192);
    // Distinct per-component input precision survives packing and output copies.
    let f=|bits,value| Float::with_val(bits,value);
    let input=vec![C::new(f(53,1)+Float::with_val(53,0.5).pow(30),f(320,2)),C::new(f(256,3),f(70,4)),C::new(f(160,5),f(80,6)),C::new(f(240,7),f(96,8))];
    let mut workspace=BatchWorkspace::new(Float::new(192));
    empty::evaluate_batch::<Float>(&[],&mut [],usize::MAX,&Context::new(table::CONSTANTS,192).unwrap(),&mut workspace).unwrap();
    assert_eq!(workspace.capacity(),0);assert_eq!(workspace.storage_len(),0);

    let ctx=Context::<Float>::new(table::CONSTANTS,192).unwrap();
    let mut output=vec![C::new(f(192,0),f(192,0));4];
    identity::evaluate_batch(&input,&mut output,2,&ctx,&mut workspace).unwrap();
    for (actual,expected) in output.iter().zip(&input) { assert_eq!(actual,expected);assert_eq!(actual.re.get_precision(),expected.re.get_precision()); assert_eq!(actual.im.get_precision(),expected.im.get_precision()); }
    let mut output=vec![C::new(f(192,0),f(192,0));6];
    ordinary::evaluate_batch(&input,&mut output,2,&ctx,&mut workspace).unwrap();
    for row in 0..2 { let scalar=ordinary::evaluate(&input[row*2..row*2+2],&ctx); for (a,b) in output[row*3..row*3+3].iter().zip(&scalar) { assert_eq!(a,b); assert_eq!(a.re.get_precision(),b.re.get_precision()); assert_eq!(a.im.get_precision(),b.im.get_precision()); } }
    // Constants come from the supplied context, never from the workspace's
    // prototype or the previous module's slots. Reuse after every stride/layout.
    for bits in [320,80,256] {
        let context=Context::<Float>::new(table::CONSTANTS,bits).unwrap();
        let mut output=vec![C::new(f(192,0),f(192,0));10];
        constants::evaluate_batch(&[],&mut output,5,&context,&mut workspace).unwrap();
        let expected=constants::evaluate(&[],&context);
        for row in output.chunks_exact(2) { for (a,b) in row.iter().zip(&expected) { assert_eq!(a,b);assert_eq!(a.re.get_precision(),b.re.get_precision());assert_eq!(a.im.get_precision(),b.im.get_precision()); } }
        ordinary::evaluate_batch(&input,&mut output[..6],2,&context,&mut workspace).unwrap();
        for row in 0..2 { assert_eq!(&output[row*3..row*3+3],&ordinary::evaluate(&input[row*2..row*2+2],&context)); }
    }
    let ctx=Context::<f64>::new(table::CONSTANTS,53).unwrap(); let mut workspace=BatchWorkspace::new(0.0);
    let input=[C::new(1.,0.),C::new(4.,0.),C::new(9.,0.)]; let mut output=[C::new(0.,0.);3];
    TRACE.lock().unwrap().clear(); traced::evaluate_batch(&input,&mut output,3,&ctx,&mut workspace).unwrap();
    let trace=TRACE.lock().unwrap().clone(); assert_eq!(trace.len(),9);
    for group in trace.chunks_exact(3) { assert!(group.iter().all(|entry|entry.0==group[0].0)); assert_eq!(group.iter().map(|entry|entry.1).collect::<Vec<_>>(),vec![1.,4.,9.]); }
    assert!(trace.iter().any(|entry|entry.0=="log"));assert!(trace.iter().any(|entry|entry.0=="sqrt"));assert!(trace.iter().any(|entry|entry.0=="abs"));
    println!("instruction-major trace={trace:?}");
    // Keep instrumentation capacity warm so it cannot confound allocator counts.
    TRACE.lock().unwrap().reserve(100_000);
    let input=[C::new(2.,0.),C::new(3.,0.)].repeat(8);let mut output=vec![C::new(0.,0.);24];
    ordinary::evaluate_batch(&input,&mut output,8,&ctx,&mut workspace).unwrap();
    COUNT.with(|count|count.set(Some(0)));
    for rows in [8,3,1,0,8] { ordinary::evaluate_batch(&input[..rows*2],&mut output[..rows*3],rows,&ctx,&mut workspace).unwrap(); }
    let allocations=COUNT.with(|count|count.replace(None).unwrap()); assert_eq!(allocations,0,"warmed f64 allocation regression");
    let mut fresh=BatchWorkspace::new(0.0); COUNT.with(|count|count.set(Some(0)));ordinary::evaluate_batch(&input[..2],&mut output[..3],1,&ctx,&mut fresh).unwrap();let allocations=COUNT.with(|count|count.replace(None).unwrap());assert_eq!(allocations,0);assert_eq!(fresh.storage_len(),0);
    println!("warmed f64 allocation count=0; fresh single-row allocation count=0; numeric/precision/Dual controls passed");
}
"#;
