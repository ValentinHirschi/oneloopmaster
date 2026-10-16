# Combined conditional-cache and C0 validation

This pin-only checkpoint starts from oneloopmaster
`c838e8ed728eb50e5a92457a1ac5520a5720e20d`, including the separately committed
native baseline regeneration and lazy C0 projective-sheet correction. All three
build roots select Symbolica, Numerica and Graphica from
`2c06f245835d57535418b5d78f711fe191db49cc` on the public
`ValentinHirschi/symbolica` integration branch. This adds the scope-local completed
conditional cache of Symbolica PR #88 to the previous `fd7bced5` integration.
No generated master file or numerical formula changes in this checkpoint.

## Checks

On 2026-10-16 the combined state passed 118 tests across these selected targets:

| Target | Passed | Ignored |
| --- | ---: | ---: |
| Library units | 83 | 1 |
| C0 derivatives | 3 | 1 |
| C0 projective control | 1 | 0 |
| Portable evaluators and JIT/cache | 2 | 0 |
| Function definitions | 8 | 0 |
| Generated native masters | 5 | 1 |
| Native contract | 4 | 0 |
| Native primitives | 11 | 0 |
| Rust batch code generation | 1 | 0 |

Build and test selection:

```sh
cargo test -j4 --features generated-evaluators,rust-codegen \
  --lib --test c0_derivatives --test c0_projective --test native \
  --test native_contract --test native_primitives --test evaluators \
  --test rust_batch_codegen --test function_definitions --no-run
# Run each resulting test binary with --test-threads=1.
```

The runtime controls used a local Cargo source override to a clean checkout of
exactly `2c06f245`, with no source differences from the public git pin. Separate
`cargo metadata` calls without overrides confirmed exactly one copy of each of
the three packages, all sourced from the public git revision, in the Rust,
Python and Python-host roots. Both Python roots passed `cargo check --locked
-j4` against the actual public git dependencies. Test runtime used
`RUST_MIN_STACK=134217728` and the externally configured Symbolica license.
Formatting and `git diff --check` also pass.

The evaluator target covers all 293 portable scalar fixtures and its batch and
cache controls. The C0 controls include the independently integrated exceptional
point, cyclic permutations, nearby rational points, Float384 higher/mixed
components and active continuation branches. The batch code-generation target
compiles emitted Rust and checks scalar/batch and lazy helper behavior.

## Limits

This is a selected regression suite, not the entire repository suite or a
performance acceptance result. The three ignored diagnostics remain unchanged.
In particular, `c0_near_projective_f64_derivative_conditioning` still records the
pre-existing failed third-order f64 derivative gate near the corrected finite
center; it is not claimed fixed. See
[c0-projective-sheet-guard.md](c0-projective-sheet-guard.md) for the unchanged
old/fixed failure and independent finite-value evidence.

The upstream update follows the user's direct-main policy for oneloopmaster;
there is no new oneloopmaster PR. Symbolica PR #88 remains a separate upstream
review. Primary NinjaRS dependency checkouts were not moved during this work.
