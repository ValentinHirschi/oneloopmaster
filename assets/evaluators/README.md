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
windings separated (`tests/d0_vanishing_gram.rs`). Both were rebuilt again on
2026-10-18 with the term-wise dilogarithm divided difference and
`log(x)/(1-x)` next to one (`tests/c0_coincident_roots.rs`) and the contour
box's exact coefficient sums (`tests/d0_tiny_invariant.rs`); both tests also
check the SymJIT backend. D0 was rebuilt once more with the contour
representation for every supported real-mass three-mass box
(`tests/d0_real_mass_orderings.rs`, which checks SymJIT too), and again with
the contour box's exact Källén discriminants, massless-end roots, exact line
sums and charts expanded about beta = 1 (`tests/d0_contour_thresholds.rs`,
likewise), and again with the far-root form of `r1` (same test).
C0 and D0 were rebuilt on 2026-10-18 with the continued dilogarithm's
unit-circle test on a real part, and D0 with the Denner-Dittmaier four-mass
chart (`tests/d0_four_mass_orderings.rs`, which checks SymJIT too). D0 was
rebuilt on 2026-10-19 with the one- and two-adjacent-mass boxes' `solabc`
roots and dilogarithm complements next to collinear limits
(`tests/d0_collinear_limits.rs`, likewise), and again with the massless
box's complements and the exact discriminants next to collinear limits
(same test). C0 and D0 were rebuilt again with the divided difference's
selection taking the dilogarithms as arguments
(`tests/expression_construction.rs`). D0 was rebuilt on 2026-10-19 with the
contour box's line-edge pole offsets from Vieta's product next to a
vanishing Gram determinant and the second contour where the first one's `K`
vanishes (`tests/d0_vanishing_gram.rs`, which checks SymJIT too), and with
Denner and Dittmaier's sheets for one-mass boxes with a complex mass
(`tests/d0_one_mass_width.rs`, likewise) and for two-mass boxes with complex
masses, whose terms are one registered function (`tests/d0_two_mass_width.rs`,
likewise). D0 was rebuilt again on 2026-10-19 with Symbolica `9cdc2949`
(alphal00p/symbolica `ninjars-pin`: `e679e8c1` plus single-rounded `Float`
products, symbolica-dev/symbolica#97) for the contour-box regression fixes
(line products only for a line nearly parallel to the contour direction, the
zero-external labelling that puts a massless corner at the origin and the
contour box's dilogarithm complements from the poles;
`tests/d0_contour_regressions.rs`, likewise) and the two-mass terms'
mass-pair complements taken from the channels' complements next to the
shells (`tests/d0_two_mass_width.rs`, `tests/data/two_mass_shell.txt`).
D0 and C0 were rebuilt on 2026-10-21 for the final candidate. D0 with the
exactly real contour norms and real sheet signs at mixed precision of
alphal00p/oneloopmaster#9, the infrared chart reclassification of #10, the
contour box's chart tests for a Gram determinant at the rounding level (a
vanishing K's far pole, the other contour where the selected one is flat,
the selected one kept where the other's K is a rounding residue, and NaN
where binary64 cannot resolve the chart;
`tests/data/contour_gram_rounding.txt`, `tests/data/contour_unresolved.txt`,
`tests/data/contour_fresh_sweep.txt`), the one- and two-mass terms'
mass-pair complements and their root nearer the branch cut taken across it,
in registered functions (`tests/d0_one_mass_width.rs`,
`tests/d0_two_mass_width.rs` with `tests/data/two_mass_conditioning.txt`),
no line product at an exactly double real pole, the line products probed
where they disagree with their pole's rounded offset or where `a1` does not
account for an offset, and kept in binary64 only where every cancelling
edge has one (`tests/data/contour_pole_errors.txt`,
`tests/data/contour_line_products.txt`), and the zero-external reflection
only for a line nearer its corner with a pole that the rounding resolves
poorly next to the corner's vertex, and not into a chart whose flat contour
is unresolved (`tests/d0_contour_regressions.rs`, which checks SymJIT too).
C0 for the first time since PR #9 changed its real sheet signs: the rebuilt
asset keeps the binary64 SymJIT values of the 10,668 C0 rows of the two
verifiers' sweeps and the NinjaRS master calls, and of the 9,516 C0 rows of
the oracle sweeps (seeds 1 to 4), bit for bit but for the sign of one NaN;
the previous asset (102,902 bytes, `f8c7109a...`) was built from the
expressions before PR #9.

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
| C0 | 85,547 | `7a7bd44280ba5dcf70d0ad9b54b93cbc29853665799a8288c8d7f1a648ac25d9` |
| D0 | 1,558,784 | `3dc226c8dccb84ee48b0d1eadcf9dd38770dc90bdcdb5e7ba5a5534d9efb9aaf` |

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
