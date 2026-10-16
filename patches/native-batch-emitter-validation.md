# Branch-free native batch emitter validation

Validated on 2026-10-16 with Rust 1.98.1 (`48a229ceaefd4985c50990b14116b6d856af0985`, LLVM 22.1.8), Linux x86_64, and `RUST_MIN_STACK=134217728`. Cargo parallelism was limited to eight jobs.

The base is `de0114add9937a2225a3e5d89af1d497ccd31d07` (the reusable scalar emitter). Symbolica, Numerica and Graphica were supplied by the unmodified Symbolica checkout at `52ef4b8e3a10cf3a68385230d55ebe8c3d153948`; SymJIT was 2.27.3. An untracked Cargo configuration patched all three Symbolica packages to that checkout. The ignored validation lockfile SHA-256 was `1ffb9302dd6e72fe6462fe5311ceb8950f6d49479f951494099c27244cb83498`. Neither manifests nor generated scalar master files changed.

## Reproduction and results

With the pinned dependencies and normal native build tools available:

```sh
export RUST_MIN_STACK=134217728
cargo test -j 8 --features rust-codegen --lib rust_codegen:: -- --test-threads=1
cargo test -j 8 --features rust-codegen --lib native:: -- --test-threads=1
cargo test -j 8 --features rust-codegen --test native_primitives -- --test-threads=1
cargo test -j 8 --features rust-codegen --test rust_batch_codegen -- --nocapture --test-threads=1
```

All 27 selected tests passed: 12 emitter unit tests, three workspace/constants tests, 11 existing numerical primitive regressions, and one compiled generated-code probe containing the numerical, precision, trace and allocation controls.

The probe checks f64, DoubleFloat and Float192 against the emitted scalar path and Symbolica's independent batch interpreter. It covers nested shared helpers, multiple outputs, exact constants beyond binary64, batch sizes 0/1/8/32/128/1024 and tails, permutations, higher/mixed flattened Dual components, and dimension/overflow errors with unchanged output. One workspace is reused across empty, constant-only, ordinary and derivative kernels, changing layouts, strides and Float contexts. Copies retain distinct real/imaginary component precision, including inputs both above and below the scratch prototype's precision. A zero-slot kernel with `usize::MAX` rows leaves capacity unchanged and does not poison the next ordinary call.

The trace for `f(x) + abs(x)`, where the non-inlined shared body is `f(x) = log(x) + sqrt(x)`, was:

```text
log(1), log(4), log(9), sqrt(1), sqrt(4), sqrt(9), abs(1), abs(4), abs(9)
```

This distinguishes instruction-major shared-helper execution from a row loop. The trace test uses the real native primitives with instrumentation wrappers. The f64 allocation counter then reports zero allocations across warmed 8/3/1/0/8-row evaluations and zero for a fresh workspace's single-row call. These claims exclude preparation, context construction, output allocation and trace-capacity setup. Float arithmetic is not claimed to be allocation-free.

A synthetic exported instruction stream separately tests an overwritten temporary, a self-assignment and a single-operand instruction. Its left-fold addition deliberately distinguishes `(1e16 - 1e16) + 1` from reordered arithmetic. Recursive rejection tests cover every control-flow variant, including unused nested bodies, and verify that failed exports leave earlier exact constant indices unchanged.

## Code generation cost

The ordinary three-output graph with nested shared helpers emitted 1,049 bytes for the scalar module and 4,311 bytes including the batch entry and helpers. The complete six-module validation program plus test harness was 99,284 bytes. On this run preparation/emission took 0.009981 seconds and `rustc -C opt-level=2` took 23.778326 seconds. The latter includes monomorphized validation harness code and linking against already-built dependencies; it is not a standalone production-kernel compile benchmark. These are descriptive single-run costs, without a timing uncertainty estimate or throughput acceptance claim.

## Existing baseline limits

An unmodified detached checkout of `de0114add9937a2225a3e5d89af1d497ccd31d07`, with the same Symbolica pin and SymJIT 2.27.3, independently reproduced the following eager initialization error:

```rust
fn main() {
    let result = oneloop::initialize();
    eprintln!("baseline eager initialize: {result:?}");
    assert!(result.is_ok());
}
```

Running that program as a Cargo example produced:

```text
could not initialize D0 evaluator: Could not JIT-compile sub-evaluator
'oneloop___olo_boxc_box': function call with 70 args, but only 64 are
available. Use `Config.set_num_args`.
```

The generated-code probe initializes Symbolica state without requesting eager master/JIT initialization. This change does not fix or conceal the baseline initialization failure and does not claim a complete master-backend validation.

Strict Clippy reports four pre-existing errors: `chunks_exact_to_as_chunks` in `src/backend.rs:171` and `src/precision.rs:189`, `large_enum_variant` in `src/precision.rs:11`, and `clone_on_copy` in `src/masters.rs:275`. No changed production file was diagnosed. This scoped command passed, as did formatting and `git diff --check`:

```sh
cargo clippy -j 8 --features rust-codegen --lib --test rust_batch_codegen -- \
  -D warnings -A clippy::chunks_exact_to_as_chunks \
  -A clippy::large_enum_variant -A clippy::clone_on_copy
```

Control flow is intentionally unsupported by `emit_batch`: branch/join/label/goto instructions are rejected recursively before emission commits its constant table. Existing scalar emission continues to support lazy branches. Native master row loops are unchanged. Lazy active-lane native batching, SIMD specialization, and matched production performance measurements remain future work.
