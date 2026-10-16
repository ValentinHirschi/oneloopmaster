# Lazy native batch emitter validation

Validated on 2026-10-16 against oneloopmaster base `69761cdf8ffc0c02dbc2e36ce41ee6fa365268ae`. The isolated branch changes the reusable emitter and its workspace; it does not replace the native master backend's existing scalar row loop or modify generated master files.

The compiler was Rust 1.98.1 (`48a229ceaefd4985c50990b14116b6d856af0985`, LLVM 22.1.8), Linux x86_64. Symbolica, Numerica and Graphica all came from the unmodified Symbolica tree `8e23925ed7b5cbb91607c8633dda114cb43be8ce`; SymJIT was 2.27.3. A local Cargo configuration patched all three packages to that tree and reused the validation target directory. The ignored lockfile SHA-256 was `1ffb9302dd6e72fe6462fe5311ceb8950f6d49479f951494099c27244cb83498`. This evidence deliberately precedes the separately developed Symbolica constant-pool correction.

## Implementation contract

`RustEmitter::emit_batch` keeps its existing interface. Branch-free graphs retain their contiguous instruction-major code, and width-one evaluations retain the scalar entry point without workspace allocation. Conditional graphs use the original exported `IfElse`, `Goto`, `Label` and `Join` instructions, with one retained active/resume frame per canonical helper. Each helper receives the parent numerical stride and the caller's active lane indices. No symbolic algebra or additional instruction representation is introduced.

All helper dependencies have earlier frame indices, so a function borrows its own control frame separately from its children. Labels reactivate only caller lanes whose forward jump has arrived. Numerical primitives and helper argument/result copies run only on active lanes. Malformed or backward control flow is rejected before the emitter commits constants. Dimension and workspace-size errors leave output buffers unchanged. Copies preserve each Float component's tracked precision through `set_from`; they do not promote inputs to the scratch prototype's precision.

## Separate scalar emitter correction

Flattening higher/mixed Dual conditionals exposed an existing scalar-emitter failure, `unbound native slot o2`. Symbolica emits consecutive `Join` instructions after one branch, while the scalar emitter previously handled only the first. The correction groups consecutive joins with the same condition, propagates each selected phi source along its own incoming edge, and preserves sequential assignment when later joins consume an overwritten target.

The focused synthetic regression verifies lazy reciprocal arms, overwritten and reused temporary slots, multiple outputs, and exactly one reciprocal per selected branch. All existing scalar-emitter unit controls and the existing production native integration controls were retained and passed after this correction. This correction is necessary even for width-one evaluation of the new flattened Dual test; it is separate from active-lane batch scheduling.

## Reproduction and results

With the dependencies above, native build tools and the Symbolica license available, the exact selected commands were equivalent to:

```sh
export RUST_MIN_STACK=134217728
cargo test --config /tmp/ninjars-oneloop-main-validation.toml --locked -j4 \
  --features rust-codegen --lib rust_codegen:: -- --nocapture --test-threads=1
cargo test --config /tmp/ninjars-oneloop-main-validation.toml --locked -j4 \
  --features rust-codegen --lib native:: -- --nocapture --test-threads=1
cargo test --config /tmp/ninjars-oneloop-main-validation.toml --locked -j4 \
  --features rust-codegen --test rust_batch_codegen -- --nocapture --test-threads=1
cargo test --config /tmp/ninjars-oneloop-main-validation.toml --locked -j4 \
  --all-features --test native --test native_contract --test native_primitives \
  -- --nocapture --test-threads=1
```

The configuration's `[patch.crates-io]` entries point `symbolica` to the pinned checkout, `numerica` to its `lib/numerica`, and `graphica` to its `lib/graphica`. The selected results were **38 passed, one existing ignored test**: 13 emitter unit tests, four workspace/constants unit tests, one compiled generated-code probe, five native integration tests, four native contract tests, and 11 native primitive tests. This is a focused validation, not a rerun of the entire repository suite.

The generated-code probe preserves all branch-free controls and adds:

- Nested non-inlined shared helpers and divergent nested branches, with zeros in inactive logarithm and reciprocal arms.
- f64, DoubleFloat and Float192 comparisons against both scalar emitted code and Symbolica's independent batch interpreter, including multiple outputs and exact registered constants.
- Rows 0, 1, 8, 32, 128 and 1024, incomplete tails, permutations, and one workspace reused across changing helper layouts and parent strides.
- Flattened Taylor components `[00, 10, 01, 20, 11, 02]` from Symbolica's Dualizer, with 24 inputs and 18 outputs. An independent analytic control for `log(x+y)+sqrt(x+y)` at `x+y=4` checks first components `1/2`, second pure components `-3/64`, and the mixed component `-3/32`.
- Distinct input component precisions from 53 to 320 bits and contexts at 80, 256 and 320 bits. Result values and reported real/imaginary precisions agree exactly with the scalar emitted path.
- Malformed dimensions, multiplication overflow and unchanged output sentinels. Workspace tests also cover huge zero-work requests, unchanged high-water metadata after overflow, and reuse after ordinary contiguous preparation.

For six caller lanes in a workspace retaining stride 128, the primitive trace was:

```text
log(4), log(9), sqrt(4), sqrt(9), powi(2, -1), powi(4, -1)
```

Only lanes 0 and 3 enter the logarithm/square-root helper. The two reciprocal lanes run afterward; the remaining lanes skip that helper entirely. This trace distinguishes instruction-major helper execution from a scalar row loop and proves that the inactive singular arms in this fixture were not executed.

The f64 allocation counter reports zero allocations over warmed conditional calls with rows 6, 3, 1, 0 and 6, and zero allocations for a fresh workspace's width-one call. Preparation, context construction, output allocation and trace-buffer reservation are excluded. Float arithmetic is not claimed to be allocation-free.

## Source and compilation cost

The original branch-free fixture still emits 1,049 bytes of scalar code and 4,311 bytes including batch code. The new ordinary conditional fixture emits 10,013 bytes including scalar and batch code; its flattened higher/mixed Dual version emits 234,897 bytes. The complete generated validation program was 351,238 bytes. Preparation/emission took 0.045579 seconds and `rustc -C opt-level=2` took 68.509142 seconds, including the monomorphized test harness and linkage to existing dependencies. These are descriptive single-run costs, not matched production performance evidence.

Scoped strict Clippy passed with the same three pre-existing lint categories allowed as in the branch-free validation:

```sh
cargo clippy --config /tmp/ninjars-oneloop-main-validation.toml --locked -j4 \
  --features rust-codegen --lib --test rust_batch_codegen -- \
  -D warnings -A clippy::chunks_exact_to_as_chunks \
  -A clippy::large_enum_variant -A clippy::clone_on_copy
cargo fmt --all --check
git diff --check
```

The pinned Numerica dependency still reports its existing unused-import and unnecessary-unsafe warnings. Native master backend integration, SIMD specialization, broader performance measurements and production throughput acceptance remain separate work. No new oneloopmaster PR is opened; this reviewable commit follows the user's direct-main consolidation policy.
