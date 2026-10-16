# Shared-helper stack integration pin

This update starts from oneloopmaster `5091653086b32428c6e024ee89404b4b21214e88` and pins all three standalone build roots to Symbolica integration `fd7bced5307e94e994d5deba63c60b8380aa7d02`. The contributor branch `codex/oneloopmaster-integration` was fast-forwarded from `58b72d4`; the added implementation is [Symbolica PR 87](https://github.com/symbolica-dev/symbolica/pull/87).

The patch runs Symbolica's existing stack allocator once on freshly lowered shared-function bodies. Explicit arguments and captured inputs remain reserved. It neither introduces a register allocator nor changes integral expressions. The previous exact Dual constant-pool, vector branch-scope and polylog series fixes remain included. No oneloopmaster runtime Rust source or generated evaluator changes in this pin update.

## Dependency checks

`cargo metadata --format-version 1` for the Rust, Python and Python-host manifests resolves exactly one Symbolica, Numerica and Graphica, all at the published git revision above, and SymJIT 2.27.3. Both tracked Python lockfiles change only those three git source revisions; package versions are unchanged. The fetched git source and the clean local Symbolica integration have identical tree `ec93de75f96cc3b4d2530b3a27c2cc7e024bdac0`.

Both of these checks pass using the actual public git source without local path overrides:

```sh
cargo check --locked -j4 --manifest-path python/Cargo.toml
cargo check --locked -j4 --manifest-path python/host/Cargo.toml
```

These validate compilation, not Python extension import or loading.

## Combined Rust validation

Validation uses Rust 1.98.1, Linux x86_64, `RUST_MIN_STACK=134217728`, four Cargo jobs and serial test execution. The selected unit/integration targets are compiled together with all features and run in separate processes to preserve initialization and cache isolation.

The reproducible Cargo command is:

```sh
cargo test --locked -j4 --all-features --lib \
  --test evaluators --test master_callbacks_jit --test initialization \
  --test rust_batch_codegen --test c0_derivatives \
  --test native --test native_contract --test native_primitives \
  --test native_evaluator --test native_hooks --test native_hook_errors \
  --test native_startup --test function_definitions \
  -- --test-threads=1 --nocapture
```

All selected targets pass: **123 passed, 5 existing ignored, zero failures**. This is a combined selected-target run, not a claim that every repository test target was rerun.

| Target | Passed | Existing ignored |
|---|---:|---:|
| `oneloop` | 83 | 1 |
| `evaluators` | 2 | 0 |
| `master_callbacks_jit` | 1 | 0 |
| `initialization` | 1 | 1 |
| `rust_batch_codegen` | 1 | 0 |
| `c0_derivatives` | 2 | 0 |
| `native` | 5 | 1 |
| `native_contract` | 4 | 0 |
| `native_primitives` | 11 | 0 |
| `native_evaluator` | 2 | 1 |
| `native_hooks` | 1 | 0 |
| `native_hook_errors` | 1 | 0 |
| `native_startup` | 1 | 1 |
| `function_definitions` | 8 | 0 |

The ignored entries are existing diagnostic checks and fresh-process child helpers. The two startup parent tests execute their six and three child entry routes respectively; these routes passed. The checks cover all 293 original scalar master fixtures, restored JIT caches, mixed batches, native DoubleFloat/Float controls, output transactionality and generated instruction-major primitive/allocation traces. C0 controls verify degree-five complex scaling and seven Taylor components, including third and mixed orders, at four nonsingular massive spacelike points with equal/unequal masses and both Källén signs. Existing tolerances and input domains are unchanged.

`cargo fmt --all --check` and `git diff --check` also pass. No production code changed that would require a new Clippy run beyond the upstream helper patch's completed checks. Both Python compilation checks above pass after the pin update.

## Storage and remaining limitations

Before this pin, the isolated Symbolica patch passed 51 upstream-main evaluator tests and 67 tests on the existing integration, including interpreter, batch, JIT, C++ export, bincode and higher/mixed Dual controls. The oneloopmaster Rust batch probe and all 293 original scalar master fixtures also passed against the patched source.

A separate generation-only audit reduces the emitted D0 numerical arena from 210,393 to 2,405 complex cells per batch lane, or about 3.21 GiB to 37.6 MiB for 1,024 complex f64 lanes. The generated D0 source remains approximately 59.4 MB. That large module was not compiled or enabled in production. These figures describe storage/source size, not an accepted runtime performance improvement. Production native master batching and broader source-size work remain separate tasks.

Command logs, exact build/runtime arguments, dependency metadata, file hashes and per-target return codes are retained under `/tmp/ninjars-oneloop-helper-pin-validation/`. Reproduction uses the published pin and commands above; the temporary directory is evidence rather than a prerequisite. This update follows the user's direct-main policy for oneloopmaster and opens no oneloopmaster pull request.
