# Portable scalar evaluators

Regenerated on 2026-10-17 from the expressions with the continuation fixes of
`tests/permutation_regressions.rs` (Symbolica snapshot `e679e8c1` and SymJIT 2.27.3 pinned in
`Cargo.toml`). The generator validated all five families against all 293
acceptance and 31 benchmark rows before writing them; the regression rows,
`tests/data/real_mass_limits.txt` and `tests/native_d0_four_mass_defect.rs`
check the SymJIT backend as well. The C0 and D0 assets were rebuilt afterwards
from the compact definitions (registered r-function and root helpers): C0 with
the exact zero coefficient of the two-mass triangle at a vanishing invariant
(`tests/float_branch_decisions.rs`), D0 with the contour residue's sheet
windings separated (`tests/d0_vanishing_gram.rs`).

Previous generations: 2026-09-29 with released Symbolica and Numerica 3.0.1 and
SymJIT 2.26.0, and the C0 asset of 2026-09-22 (zero-root ordering of the finite
three-mass triangle, with 45 equal-mass regression points through the native,
expression and portable SymJIT backends).

These files contain bincode-serialized numerical `ExpressionEvaluator` graphs,
including nested definitions and registered callback metadata. They contain no
native machine code. Loading recompiles every level for the host with OneLOop's
explicit settings: O2, ordinary translation, packed complex arithmetic enabled,
and fastmath, SIMD batches, threads and AVX-512 disabled. Registered numerical
callbacks preserve complex square-root range/branches and exact predicates.
Packed complex arithmetic retains the validated production policy; the earlier
strict-FMA defect is fixed in SymJIT 2.26.0.

The generator restored each file and checked all three Laurent coefficients
against all 293 acceptance and 31 benchmark rows before writing it. All five
families passed. Fresh-process, batch and Python validation is tracked in the
[migration notes](../../SYMBOLICA_3_MIGRATION.md). Native x86-64 execution is tested;
other architectures have not been exercised here. These fixtures do not establish
exhaustive analytic-region coverage.

| Family | Bytes | SHA-256 |
| --- | ---: | --- |
| A0 | 1,488 | `b034aa19e48f3f62b7029f2c88ad4f14e19ff6d9580979d43868349c5e651a0d` |
| B0 | 15,375 | `0730c8a390d31c508fbece0760186627dc851b8ccf2da787072f8c2788f50354` |
| dB0 | 61,566 | `5f145ba984ffe2bdffb9355557a9272c3006d0a601e845cd48877d5d19ea4d77` |
| C0 | 86,330 | `1dd3da62d2fabbefa5aa4008181f3fbc74266d0c66dcaf7ad5e3f92f91913012` |
| D0 | 1,861,137 | `98dad062fd7a915136d52ebf4a109f315d82934473182bdc77c1ff1927863f40` |

The embedded files use the `oneloop-evaluator-v3` discriminator, which binds
the Symbolica IR schema and strict compilation policy without claiming backend
binary compatibility; loading always recompiles the Symbolica IR with the
current backend. The known `oneloop-evaluator-v2` source formats of earlier
generations are still accepted for migration. Original helper symbols and signatures are preserved. Upstream shares
non-inlined bodies, so explicit-pi aliases and padded argument lists are no
longer needed. Symbol ordering still controls floating-point operation order.

Regenerate with:

```sh
cargo run --release --no-default-features --example rebuild_evaluators -- assets/evaluators
```

An optional final family name writes only that family. Startup still prepares
all five evaluators. Validation failure leaves the affected existing file
untouched. Commit matching source and regenerated assets together, and update
this manifest's sizes and hashes. Load only trusted serialized evaluators.

The previous generation and its dependency patches are documented in the
[2026-09-08 audit](../../performance/2026-09-08-native/README.md).
