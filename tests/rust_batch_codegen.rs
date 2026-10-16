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

fn lazy_exact() -> symbolica::evaluate::ExpressionEvaluator<Complex<Rational>> {
    Atom::evaluator_multiple(
        &[
            symbolica::parse!(
                "if(c,batch_lazy::f(d,x,y),x+pi+polylog(2,1/2)+1/100000000000000000000)"
            ),
            symbolica::parse!("if(c,x+y,x-y)"),
            symbolica::parse!("if(c,if(d,x,y),pi)"),
        ],
        &[
            symbolica::parse!("c"),
            symbolica::parse!("d"),
            symbolica::parse!("x"),
            symbolica::parse!("y"),
        ],
    )
    .add_function_with_options(
        symbolica::symbol!("batch_lazy::g"),
        vec![symbolica::symbol!("z")],
        symbolica::parse!("log(z)+sqrt(z)"),
        FunctionRegistrationOptions::new().inlining(InliningPolicy::Never),
    )
    .unwrap()
    .add_function_with_options(
        symbolica::symbol!("batch_lazy::f"),
        vec![
            symbolica::symbol!("a"),
            symbolica::symbol!("u"),
            symbolica::symbol!("v"),
        ],
        symbolica::parse!("if(a,batch_lazy::g(u+v),1/v)"),
        FunctionRegistrationOptions::new().inlining(InliningPolicy::Never),
    )
    .unwrap()
    .cores(1)
    .build()
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
    let lazy = emitter
        .emit_batch(lazy_exact(), EmitOptions::default())
        .unwrap();
    let lazy_dual = emitter
        .emit_batch(
            lazy_exact()
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
                .unwrap(),
            EmitOptions::default(),
        )
        .unwrap();
    assert_eq!((lazy_dual.input_count, lazy_dual.output_count), (24, 18));
    eprintln!(
        "lazy source: scalar_plus_batch_bytes={} flattened_dual_bytes={}",
        lazy.source.len(),
        lazy_dual.source.len()
    );
    // Laurent triangle-to-bubble subtraction contains c*d1^2/d0^3. For
    // c=d1=d0=1e110 the value is one but binary64 d0^3 overflows. Its
    // reciprocal must not conceal this overflow as a finite zero.
    let correction = emitter
        .emit_batch(
            Atom::evaluator_multiple(
                &[
                    symbolica::parse!("d1^2/d0^3"),
                    symbolica::parse!("c*d1^2/d0^3"),
                ],
                &[
                    symbolica::parse!("c"),
                    symbolica::parse!("d1"),
                    symbolica::parse!("d0"),
                ],
            )
            .build()
            .unwrap(),
            EmitOptions::default(),
        )
        .unwrap();
    let mut source = String::from(PROBE);
    for (name, module) in [
        ("correction", correction),
        ("ordinary", ordinary),
        ("derivatives", dual),
        ("identity", identity),
        ("constants", constants),
        ("traced", traced),
        ("empty", empty),
        ("lazy", lazy),
        ("lazy_derivatives", lazy_dual),
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
use symbolica::{atom::{Atom, AtomCore}, domains::{dual::HyperDual, float::{Complex, DoubleFloat, Float, FloatLike, Real, RealLike, SingleFloat}, rational::Rational}, evaluate::{Dualizer, FunctionRegistrationOptions, InliningPolicy}};
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
fn lazy_exact() -> symbolica::evaluate::ExpressionEvaluator<Complex<Rational>> {
    Atom::evaluator_multiple(&[
        symbolica::parse!("if(c,batch_lazy::f(d,x,y),x+pi+polylog(2,1/2)+1/100000000000000000000)"),
        symbolica::parse!("if(c,x+y,x-y)"),
        symbolica::parse!("if(c,if(d,x,y),pi)")
    ], &[symbolica::parse!("c"),symbolica::parse!("d"),symbolica::parse!("x"),symbolica::parse!("y")])
    .add_function_with_options(symbolica::symbol!("batch_lazy::g"),vec![symbolica::symbol!("z")],symbolica::parse!("log(z)+sqrt(z)"),FunctionRegistrationOptions::new().inlining(InliningPolicy::Never)).unwrap()
    .add_function_with_options(symbolica::symbol!("batch_lazy::f"),vec![symbolica::symbol!("a"),symbolica::symbol!("u"),symbolica::symbol!("v")],symbolica::parse!("if(a,batch_lazy::g(u+v),1/v)"),FunctionRegistrationOptions::new().inlining(InliningPolicy::Never)).unwrap()
    .cores(1).build().unwrap()
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
fn lazy_cases() -> [[i64;4];6] { [[1,1,4,0],[1,0,0,2],[0,1,0,0],[1,1,9,0],[1,0,0,4],[0,0,-2,0]] }
fn check_lazy<T:NativeFloat>(bits:u32) where C<T>:symbolica::evaluate::EvaluationDomain {
    let prototype=T::zero_at(bits);let n=|v:i64| C::new(prototype.from_i64(v),prototype.zero());
    let ctx=Context::<T>::new(table::CONSTANTS,bits).unwrap();let mut workspace=BatchWorkspace::new(prototype.clone());
    let mut oracle=lazy_exact().map_coeff_with_prec(&|c| C::new(prototype.from_rational(&c.re),prototype.from_rational(&c.im)),bits).into_batch_evaluator();
    for rows in [0,1,8,32,128,1024,17,3,8,2,1,0] {
        let input=(0..rows).flat_map(|i|lazy_cases()[(i*5+rows)%6].map(n)).collect::<Vec<_>>();
        let mut actual=vec![n(99);rows*3];let mut expected=actual.clone();
        lazy::evaluate_batch(&input,&mut actual,rows,&ctx,&mut workspace).unwrap();oracle.evaluate_batch(&input,&mut expected,rows).unwrap();
        for row in 0..rows { let scalar=lazy::evaluate(&input[row*4..row*4+4],&ctx);assert_eq!(&actual[row*3..row*3+3],&scalar); }
        for (a,b) in actual.iter().zip(&expected) { let error=(a.clone()-b.clone()).norm().re.to_f64();assert!(error/(1.0+b.norm().re.to_f64()) < 2.0_f64.powi(-(bits as i32)+14),"lazy interpreter discrepancy {error}"); }
    }
    assert_eq!(workspace.capacity(),1024);
    let mut output=vec![n(99);6];let before=output.clone();
    for (input,rows) in [(vec![n(1)],2),(vec![],usize::MAX),(vec![],0)] { assert!(lazy::evaluate_batch(&input,&mut output,rows,&ctx,&mut workspace).is_err());assert_eq!(output,before); }
    let shape=vec![vec![0,0],vec![1,0],vec![0,1],vec![2,0],vec![1,1],vec![0,2]];
    let mut oracle=lazy_exact().vectorize(&Dualizer::new(HyperDual::<Complex<Rational>>::new(shape),vec![])).unwrap().map_coeff_with_prec(&|c| C::new(prototype.from_rational(&c.re),prototype.from_rational(&c.im)),bits).into_batch_evaluator();
    for rows in [1,8,32,3] {
        let input=(0..rows).flat_map(|i| {let v=lazy_cases()[i%6]; [n(v[0]),n(0),n(0),n(0),n(0),n(0),n(v[1]),n(0),n(0),n(0),n(0),n(0),n(v[2]),n(1),n(0),n(0),n(0),n(0),n(v[3]),n(0),n(1),n(0),n(0),n(0)]}).collect::<Vec<_>>();
        let mut actual=vec![n(99);rows*18];let mut expected=actual.clone();
        lazy_derivatives::evaluate_batch(&input,&mut actual,rows,&ctx,&mut workspace).unwrap();oracle.evaluate_batch(&input,&mut expected,rows).unwrap();
        for row in 0..rows { assert_eq!(&actual[row*18..row*18+18],&lazy_derivatives::evaluate(&input[row*24..row*24+24],&ctx)); }
        for (a,b) in actual.iter().zip(&expected) {let error=(a.clone()-b.clone()).norm().re.to_f64();assert!(error/(1.0+b.norm().re.to_f64()) < 2.0_f64.powi(-(bits as i32)+14),"lazy Dual discrepancy {error}");}
        // Independent control: g(x+y)=log(x+y)+sqrt(x+y) at x=4,y=0.
        // Taylor xx/yy coefficients are g''/2=-3/64; the mixed one is -3/32.
        let rational=|a:i64,b:i64| prototype.from_rational(&(Rational::from(a)/Rational::from(b)));
        assert_eq!(actual[1].re,rational(1,2));assert_eq!(actual[2].re,rational(1,2));
        assert_eq!(actual[3].re,rational(-3,64));assert_eq!(actual[4].re,rational(-3,32));assert_eq!(actual[5].re,rational(-3,64));
    }
}

fn check_correction<T: NativeFloat>(bits: u32, from_f64: impl Fn(f64) -> T) {
    let context = Context::<T>::new(table::CONSTANTS,bits).unwrap();
    let mut workspace = BatchWorkspace::new(T::zero_at(bits));
    let expressions=[symbolica::parse!("d1^2/d0^3"),symbolica::parse!("c*d1^2/d0^3")];
    let mut oracle = Atom::evaluator_multiple(&expressions,&[symbolica::parse!("c"),symbolica::parse!("d1"),symbolica::parse!("d0")]).build().unwrap()
        .map_coeff_with_prec(&|c| C::new(Float::new(256).from_rational(&c.re),Float::new(256).from_rational(&c.im)),256).into_batch_evaluator();
    for scale in [1.,1e20,4e102,8e102,1e110,-1e110] {
        for rows in [1,3] {
            let input=vec![C::new(from_f64(scale),from_f64(0.));rows*3];
            let exact_input=vec![C::new(Float::with_val(256,scale),Float::new(256));rows*3];
            let mut expected=vec![C::new(Float::new(256),Float::new(256));rows*2];
            oracle.evaluate_batch(&exact_input,&mut expected,rows).unwrap();
            let mut actual=vec![C::new(from_f64(99.),from_f64(0.));rows*2];
            correction::evaluate_batch(&input,&mut actual,rows,&context,&mut workspace).unwrap();
            for row in 0..rows {
                let scalar=correction::evaluate(&input[row*3..row*3+3],&context);
                assert!(expected[row*2].is_finite() && expected[row*2+1].is_finite());
                assert!((expected[row*2+1].re.clone()-Float::with_val(256,1)).norm() < Float::with_val(256,1e-70));
                for column in 0..2 {
                    let analytic=if column==0{from_f64(1.)/from_f64(scale)}else{from_f64(1.)};
                    for result in [&actual[row*2+column],&scalar[column]] {
                        if result.is_finite() {
                            let relative=(result.clone()-C::new(analytic.clone(),from_f64(0.))).norm().re/analytic.norm();
                            let tolerance=if T::FIXED_BITS.is_none(){1e-65}else{1e-12};
                            assert!(relative.to_f64()<tolerance,"overflow hidden as finite value: bits={bits} scale={scale:e} value={result:?}");
                        } else {
                            assert!(T::FIXED_BITS.is_some() && scale.abs()>4e102,"unexpected nonfinite correction: bits={bits} scale={scale:e} value={result:?}");
                        }
                    }
                    for (a,b) in [(&actual[row*2+column].re,&scalar[column].re),(&actual[row*2+column].im,&scalar[column].im)] { assert!(a==b || (a.real_classify()==Some(std::num::FpCategory::Nan) && b.real_classify()==Some(std::num::FpCategory::Nan))); }
                }
            }
        }
    }
}

fn main() {
    let _ = symbolica::symbol!("batch_codegen_state");
    check_correction::<f64>(53,|x|x);
    check_correction::<DoubleFloat>(106,DoubleFloat::from);
    check_correction::<Float>(256,|x|Float::with_val(256,x));
    check::<f64>(53); check::<DoubleFloat>(106);check::<Float>(192);
    check_lazy::<f64>(53);check_lazy::<DoubleFloat>(106);check_lazy::<Float>(192);
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
    let varied=lazy_cases()[..3].iter().flatten().enumerate().map(|(i,v)| C::new(Float::with_val([53,70,160,256,320][i%5],*v),Float::with_val([80,192,240][i%3],0))).collect::<Vec<_>>();
    for bits in [320,80,256] {
        let context=Context::<Float>::new(table::CONSTANTS,bits).unwrap();let mut result=vec![C::new(f(192,0),f(192,0));9];
        lazy::evaluate_batch(&varied,&mut result,3,&context,&mut workspace).unwrap();
        for row in 0..3 {let scalar=lazy::evaluate(&varied[row*4..row*4+4],&context);for (a,b) in result[row*3..row*3+3].iter().zip(&scalar) { assert_eq!(a,b);assert_eq!(a.re.get_precision(),b.re.get_precision());assert_eq!(a.im.get_precision(),b.im.get_precision()); }}
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
    // Retained stride exceeds current rows; only lanes0/3 reach the nested
    // logarithm/square-root helper. Inactive log(0) and reciprocal(0) must not run.
    let input=lazy_cases().into_iter().flatten().map(|v|C::new(v as f64,0.)).collect::<Vec<_>>();let mut result=vec![C::new(0.,0.);18];
    workspace.prepare(128,1).unwrap();TRACE.lock().unwrap().clear();
    lazy::evaluate_batch(&input,&mut result,6,&ctx,&mut workspace).unwrap();
    let trace=TRACE.lock().unwrap().clone();
    let logarithms=trace.iter().filter(|(name,_)|*name=="log").copied().collect::<Vec<_>>();
    let roots=trace.iter().filter(|(name,_)|*name=="sqrt").copied().collect::<Vec<_>>();
    assert_eq!(logarithms,vec![("log",4.),("log",9.)]);assert_eq!(roots,vec![("sqrt",4.),("sqrt",9.)]);
    let grouped=trace.iter().filter(|(name,_)|*name=="log"||*name=="sqrt").copied().collect::<Vec<_>>();assert!(grouped[..2].iter().all(|x|x.0==grouped[0].0));
    assert!(!trace.iter().any(|(name,value)|(*name=="log"||*name=="powi")&&*value==0.));
    println!("lazy active-lane trace={trace:?}; retained stride={}",workspace.capacity());
    let mut fresh=BatchWorkspace::new(0.0);TRACE.lock().unwrap().reserve(100_000);
    COUNT.with(|count|count.set(Some(0)));lazy::evaluate_batch(&input[..4],&mut result[..3],1,&ctx,&mut fresh).unwrap();let allocations=COUNT.with(|count|count.replace(None).unwrap());assert_eq!(allocations,0);assert_eq!(fresh.storage_len(),0);
    COUNT.with(|count|count.set(Some(0)));for rows in [6,3,1,0,6] {lazy::evaluate_batch(&input[..rows*4],&mut result[..rows*3],rows,&ctx,&mut workspace).unwrap();}let allocations=COUNT.with(|count|count.replace(None).unwrap());assert_eq!(allocations,0,"lazy warmed allocation regression");
    println!("lazy warmed allocation count=0; inactive singular arms skipped");

}
"#;
