# Native batch source outlining and production integration audit

This change extracts generated Rust source at existing exported instruction boundaries. It does not alter symbolic expressions, introduce an instruction representation, change production dispatch, or replace native master fast sectors. Its base is `df99fb48c2f87a6fec9199cccefc83587b381737`, validated with Symbolica/Numerica/Graphica `8e23925ed7b5cbb91607c8633dda114cb43be8ce`, SymJIT 2.27.3 and Rust 1.98.1.

## Extraction contract

A batch body is split after at most 128 exported instructions or after its accumulated source reaches 32 KiB, always between complete instructions. A single large-arity instruction can exceed the byte threshold. Small bodies remain byte-identical. The existing emitter `Outlines` owns canonical source identities and extracted functions; no arithmetic is reordered or recomputed. Extracted functions are marked `inline(never)` to keep optimizer input bounded.

Numerical values remain in the original arena. Conditional blocks borrow the caller's active/resume arrays and pass back the active count. Labels retain their original absolute instruction positions, including when a branch begins in one extracted block and rejoins in another. Shared helper calls retain parent stride and their disjoint child control frames. The coordinator remains a sequence of block calls; sufficiently huge coordinators may themselves warrant a later extraction level.

The new compiled regression crosses chunk boundaries inside both selected arms and across labels/joins. It checks aliased consecutive joins, a retained stride larger than the current batch, reversed lanes, and inactive reciprocal-of-zero operands. The original scalar liveness regressions remain unchanged.

## Master generation-only feasibility audit

The original `examples/generate_native.rs` was copied into an untracked inspection example, changing only `emitter.emit` to `emitter.emit_batch` and adding source/layout logging. The normal generator and all production files were left unchanged. Both sets were generated into an external target directory; the large D0 source was not compiled.

| Family | Source bytes before | Source bytes after | Largest batch function before | Largest batch function after | Numerical slots |
|---|---:|---:|---:|---:|---:|
| A0 | 17,026 | 17,026 | 10,405 | 10,405 | 92 |
| B0 | 109,989 | 111,839 | 88,198 | 19,525 | 641 |
| dB0 | 473,493 | 482,029 | 425,851 | 19,813 | 2,895 |
| C0 | 2,788,973 | 2,825,781 | 974,000 | 33,946 | 11,504 |
| D0 | 59,332,917 | 60,075,734 | 34,262,205 | 97,040 | 210,393 |

This improves individual function size, not total source size or numerical storage. C0 needs about 180 MiB and D0 about 3.21 GiB of numerical workspace at 1,024 complex f64 rows, excluding control arrays and Float payloads. The D0 coordinator remains the largest function after extraction. These figures rule out blindly switching the production backend to the current general D0 batch export.

The shared-helper storage audit found that Symbolica already has `ExpressionEvaluator::optimize_stack`. Root evaluators use it, but the pinned fresh `FunctionBody` construction stores `param_count + instruction_count` stack cells without applying that existing allocator. Applying the owning allocator to fresh helper SSA is a separate upstream candidate; no register allocator or master-specific algebra is added here.

## Focused validation

With `RUST_MIN_STACK=134217728`, the same pinned local Cargo configuration as the lazy-batch validation, and at most four Cargo jobs:

```sh
cargo test --config /tmp/ninjars-oneloop-main-validation.toml --locked -j4 \
  --features rust-codegen --lib rust_codegen:: -- --nocapture --test-threads=1
cargo test --config /tmp/ninjars-oneloop-main-validation.toml --locked -j4 \
  --features rust-codegen --test rust_batch_codegen -- --nocapture --test-threads=1
cargo clippy --config /tmp/ninjars-oneloop-main-validation.toml --locked -j4 \
  --features rust-codegen --lib --test rust_batch_codegen -- \
  -D warnings -A clippy::chunks_exact_to_as_chunks \
  -A clippy::large_enum_variant -A clippy::clone_on_copy
```

All 14 emitter unit tests and the compiled generated-code integration probe passed. Scoped strict Clippy passed with the same three existing baseline lint allowances. The integration probe retains f64/DoubleFloat/Float, higher/mixed flattened Dual, primitive tracing, inactive singular branches, error transactionality, arbitrary parent stride and zero warmed f64 allocation controls. The ordinary small branch-free source remains 1,049 bytes scalar / 4,311 bytes scalar plus batch. The lazy Dual module grows from 234,897 to 236,377 bytes; its complete optimized validation program compiled in 46.001928 seconds on this run versus 68.509142 seconds in the preceding lazy-batch validation. These are descriptive, contended runs rather than a controlled compilation or throughput benchmark.

## NinjaRS large-kernel compilation and numerical control

An isolated copy of the NinjaRS generator and its existing numerical probe used this emitter with the same pinned Symbolica integration. The generated plans are `LaurentDivision::new(4, 5, 3)`, `TriangleBubbleSubtraction::new(true)` and `DenominatorExpansion::new(2)`, each in ordinary and flattened Dual form. The Dual shape is `[[0,0], [1,0], [0,1], [2,0], [1,1]]`. The unchanged exported subtraction Dual graph has 8,605 instructions, 165 inputs and 65 outputs.

All six generated scalar portions remain byte-identical to the preceding unoutlined source. The subtraction Dual module contains 1,503,903 bytes total, including its unchanged 314,109-byte scalar portion. Its 69 batch functions have a largest body of 17,411 bytes; the numerical workspace still uses 640 slots.

The full development-profile probe compiled successfully in 57.099648 seconds, with maximum resident memory 2,362,496 KiB. The preceding unoutlined probe was stopped after 1,364 seconds at 5,578,348 KiB. Compiler caches and process contention were not controlled identically, so these results establish compilation feasibility rather than a speedup ratio or performance acceptance.

The executable passed all 108 numerical groups: three plans, f64/DoubleFloat/Float256, and batch sizes `[0,1,8,32,128,1024,3,31,127,0,8,1]`. Each group compares native scalar and instruction-major batch results against the existing Symbolica evaluator, including all five Dual components. The unchanged tolerances are `1e-12`, `1e-25` and `1e-65` respectively. The probe also retains permutations, malformed/overflowing dimensions with unchanged outputs, warmed f64 allocation checks and the fresh width-one path.

The isolated generation/build commands were:

```sh
cargo run --manifest-path /shared/external/ninja/target/drafts/ninjars-batch-outlining/Cargo.toml \
  --locked -j4 -p ninjars-symbolica --features rust-codegen \
  --example generate_rust_kernels -- \
  /shared/external/ninja/target/batch-outline-inspection/ninjars
cargo build --manifest-path /shared/external/ninja/target/drafts/ninjars-batch-outlining/Cargo.toml \
  --locked -j4 -p ninjars-symbolica --features rust-codegen \
  --example outlined_kernel_probe
```

The temporary `outlined_kernel_probe` example only includes the existing probe source and calls `probe::run()`. Its source hashes, dimensions and measurements are retained in `/shared/external/ninja/target/batch-outline-inspection/ninjars/outlined-build.json`, with all numerical records in `numeric.json`. Root NinjaRS sources and dependency pins were not changed. Production dispatch and performance acceptance both remain false.

## Minimal later production integration

1. Make the retained batch workspace available under `generated-evaluators` as well as `rust-codegen`; preserve the current feature defaults. Store it in `NativeEvaluator`, with independent clone ownership.
2. Preserve whole-batch domain validation before writes and before Float rounding. `NativeEvaluator<Float>` intentionally rounds every input component to the evaluator's fixed working bits; raw emitted kernels intentionally do not. Round into a retained batch buffer at this existing public boundary.
3. Preserve A0/B0 scaleless handling and the existing compact massless C0/D0 formulas. Exclude these rows from the general generated path, retain their order, and stage outputs until the generated path succeeds. A first integration retaining row evaluation of compact sectors must disclose that limitation; it must not claim every sector's arithmetic is instruction-major.
4. Send the remaining rows through generated batch entries with reused buffers. Keep the width-one path exactly as it is, and preserve all reported output precision and signed-zero behavior.
5. In `PrecisionEvaluator::evaluate_batch`, delegate its `Native` variant to `NativeEvaluator::evaluate_batch`; otherwise the higher-level Float API would still wrap the new path in a scalar loop. Keep the expression variant's existing behavior.
6. Leave the four-entry f64 coefficient grouping cache and per-family locks unchanged. Manual batches currently invalidate coefficient groups once, while per-tag callbacks remain scalar. Float callback caches remain keyed by all input components' maximum precision. NativeEvaluator's `NativeFloat` contract does not itself admit Dual carriers; flattened Dual codegen remains a separate path.
7. Validate all 293 original fixtures in mixed batches, the 342-row expanded Float controls, scalar/batch component precision and signed-zero parity, compact/general mixtures, invalid late rows with unchanged outputs, clone independence, empty/tail/repeated batches, and hook-cache usability. Benchmark actual batches separately; `examples/numeric_performance.rs` currently times scalar calls.

Native D0 dispatch, a general source-size reduction, helper-slot reuse and performance acceptance remain outstanding. No main dependency pins or root NinjaRS sources were edited, and no upstream push or PR is performed by this worktree.
