# Native OneLOop status

This page preserves the September 8 development snapshot. See the
[Symbolica main migration report](SYMBOLICA_3_MIGRATION.md) for the current
dependency setup and validation.

## One-mass boxes with a complex mass (2026-10-19)

- **OneLOop's `boxf1` is wrong for a complex mass.** It normalizes the
  roots `X` of its quadratic by the complex `sqrt(m²)` and takes the
  logarithm of a dilogarithm argument `X r34/r13` as `log X + log r34 -
  log r13`. The normalization turns `X` by half the mass's phase, across the
  negative real axis for a root next to it, while the integral is analytic
  there. In `4 1.0 P 14.78 12.26 -18.32 15.39 -17.48 0 0 0 0 0 0 3.018
  -0.302` OneLOop 3.6 (double and quadruple precision) jumps between
  P = -0.06 and -0.05 from -0.0250 - 3.2e-4i to 0.0254 + 1.87i, and gives
  6.11 + 383i at P = -1e-6 for -0.0599 - 8.2e-4i. Of 3000 random one-mass
  boxes (widths 1e-3 to 0.5 of the mass) 131 are wrong by 0.3 to 7 of the
  value, and the same 131 with the widths divided by 1e6: OneLOop's D0 does
  not tend to its real-mass value there.
- **The expressions take Denner and Dittmaier's sheets** for a mass with a
  width (arXiv:1005.2076, eqs. (3.82)-(3.83), stated there for any real
  momenta and complex masses): the unnormalized roots `u = X sqrt(m²)/|sqrt(m²)|`
  and the ratios `Y23/Y02`, `Y13/Y01`, `Y03/m²` each with its principal
  logarithm (`box_integral::one_mass_complex_sheets`); real masses keep
  OneLOop's sheets and their bits. Where the two roots lie on opposite lips
  of the negative real axis (next to each other near a collinear limit),
  their principal sheets differ and the dilogarithm pairs cancel as
  `1/(X1 - X2)` (1e-16 became 6e-8 in binary64): the second root then takes
  the first one's lip, and the change of its terms, a multiple of `4 pi²`,
  is rounded and added back.
- **Validation** against COLLIER 1.2.9 in quadruple precision (its COLI and
  DD branches agree to 1e-15): all 3000 random boxes and the 24 relabellings
  of the 131 within 8e-13 in binary64 on every route (native, generated,
  interpreted, SymJIT) and 2.5e-15 at 256 bits; 1000 boxes next to a
  collinear limit and to `r13 r24 = r12 r34` at once within 5e-12 (before:
  41 wrong, the others up to 2.4e-7); with widths of 1e-11 of the mass,
  2000 boxes within 3e-9 of the real-mass expressions. A direct
  Feynman-parameter integral (primary sectors, contour deformation, nested
  tanh-sinh) reproduces the reproducer at P = -1e-3 to 2e-11 and at
  P = -1e-6 to 1e-9. `tests/d0_one_mass_width.rs` checks 15 rows
  (`tests/data/one_mass_width.txt`, references from the formula at 40
  digits) in every relabelling. Real-mass one-mass boxes and every other
  mass pattern are bitwise unchanged on every route.
- **Open: two-mass boxes with complex masses** (`boxf2`, `boxf5`) follow
  OneLOop, which is wrong in 8 to 13 of 1000 random boxes per mass pattern,
  at any width (also with the widths divided by 1e6), and right with real
  masses. Denner and Dittmaier's eqs. (3.78)-(3.80) cover them.

## Rank-four bubble tensors (2026-10-18)

`NativeTensorEvaluator` has a fourth family, `TensorIntegral::B1111`
(`[B0, B1, B11, B111, B1111]`, for numerator rank `n+2`), on every numeric
type, with a binary64 specialization; see
[NATIVE_BACKEND.md](NATIVE_BACKEND.md#rank-four-b1111). `B11` and `B111` are
bitwise unchanged on every route, including Ninja's recorded bubble calls.

- **OneLOop's `B1111` is wrong in every mass-swap branch.** OneLOop 3.6 and 3.4
  negate `B0' + 4B1' + 6B11' + 4B111' + B1111'` when they reconstruct
  `B1111(p; m0, m1)` from the swapped masses: at `p² = 0` with `|m1| >= |m0|`,
  and at `p² != 0` with `|m1| < |m0|` (massless or not). The native family and
  the exact `b1111` use the plus sign; OneLOop (and Ninja's `avholo`) must not
  serve as a `B1111` reference. `tests/native_tensor_b1111.rs` checks the swap
  identity `B_n = (-1)^n sum_k C(n,k) B_k'` in every branch.
- **The exact `b1111` loses digits at small momenta in fixed precision.** Its
  Passarino-Veltman form divides by `2p²` per rank above one: binary64 errors
  are 4e-6 at `p²/m² = 2^-8` and 0.6 at `2^-12`, DoubleFloat 4e-10 at `2^-16`
  and 2.9 at `2^-24`; `Float` needs about `4 log2(m²/p²)` extra bits.
- **Recommendation:** fixed-precision carriers (`f64`, `DoubleFloat`) take
  the bubble tensors `B1`-`B1111` from `NativeTensorEvaluator`, which keeps
  `B1111` within 2e-15 of `max_k |B_k|` in binary64 on the `2^-4..2^-24`
  ladder, and within 1e-14 down to `p² = 1e-300`.

## Current development status (2026-09-08)

Direct generic Rust evaluation is now the default for the five scalar masters.
Their Symbolica hooks use that native backend exclusively for `Complex<f64>`,
`Complex<DoubleFloat>`, and arbitrary-precision `Complex<Float>`. Manual calls
also retain explicit `symjit` and `expression` backends. Transparent native
FunctionMap definitions remain available and take precedence over bare hooks.
All master Symbols include the squared renormalization scale as their final
numeric argument.

Initialization eagerly registers the Symbols and prepares both sets of
binary64 backends: native Rust workspaces and embedded portable SymJIT
evaluators. Precision-specific Float workspaces are created on request. Python
supports decimal-digit `prec` requests, exact Decimal/large-integer conversion,
DecimalComplex inputs/outputs, and batching. The former approximately 900-bit
cap in the scalar-master Li2 path has been removed; precision-provenance fixes
also cover roots, phases, powers, zero arithmetic, and Float norms. A working
precision request is not a guarantee that every requested digit survives
cancellation or singular conditioning.

See [NATIVE_BACKEND.md](NATIVE_BACKEND.md), [PRECISION.md](PRECISION.md),
[ARBITRARY_PRECISION_AUDIT.md](ARBITRARY_PRECISION_AUDIT.md), and the
[master-Symbol audit](MASTER_SYMBOL_AUDIT.md) for implementation and audit scope.

Completed checks from this pass:

- Final standalone Python suite: **27 tests passed**, 83.755 seconds.
- Community adapter **compile-only check passed**, 25.13 seconds. Runtime
  validation inside the community host is still pending.
- Focused tests against real Symbolica: **11 primitive tests passed** in
  0.08 seconds and **3 simple-sector tests passed** in 0.06 seconds.
- The independent ArbPrec audit records a separate **5/5** fresh-Numerica
  precision-provenance gate, including 512/3456-bit controls.

The final rebuilt release suite passes **114 Rust tests**, with six ignored
entries. The explicit native binary64 diagnostic still fails on 14 finite
B0/dB0 coefficients. The matched mixed-input native/Fortran ratios are
0.85/0.84/1.08/1.22/1.51 for A0/B0/dB0/C0/D0 at batch 1024; the largest sector
miss at that batch size is timelike D0 at 3.80×. See the [current release audit](NATIVE_RELEASE_AUDIT.md) and
[matched performance report](performance/2026-09-08-native/README.md), which also
separate bare hooks, Python and manually selected SymJIT. These passing checks
do not establish global parity, guaranteed accuracy or publication readiness.

## Preserved SymJIT-default verified snapshot (2026-09-08)

The following measurements belong to the earlier SymJIT-default implementation,
before the direct Rust default and current arbitrary-precision changes.
"Current" and "default" within this preserved section describe that historical
snapshot. They do not replace the final native-default gates listed above.

The current manifest uses patched Symbolica dev `fb845d34`, with portable O2
evaluator/batch APIs and a separate PyO3 adapter. The release Rust suite passes
57 tests (four ignored test entries, including a startup subprocess helper).
The expanded 342-point / 1,026-output
audit passes at 256 bits, but fixed binary64 still fails on 14 small-momentum B0/dB0
finite coefficients. Passing the default suite is not a full-parity certificate.
Both Python adapter build modes compile. The standalone release extension passes
all 11 tests on CPython 3.13.12 in the original fresh-worker mode (19.740 s): seven
API tests and four stdlib benchmark-harness tests. The unchanged separate
two-worker regression also passes. Integrated community-host runtime validation
remains pending. Thirteen unique focused engine tests pass on the new pin,
including exact-constant JIT evaluation, portable complex/nested callback
restoration, and native special-function registration.
The published SymJIT backend required nine further generic fixes: direct power
registration, complex callback lane packing, nested SIMD fallback, exact complex
IF truth, principal square-root lips, negative integer powers, known-real
absolute-value registers, strict complex-arithmetic contractions, and zero-input
SIMD transposition. Twelve dependency regressions and the expanded eight-case
SIMD probe pass. Integrated
scalar-family strict-v3 caches have been generated and are enabled by default;
their restored scalar outputs pass all three Laurent coefficients for every
row of the 293-row acceptance and 31-row benchmark tables. The release suite
also checks mixed SIMD batches, portable restoration, clone lifetimes, malformed
cache dimensions, and constant-only zero-input batching. The bare-Symbol outer
JIT route, with no native FunctionMap, passes all 293 acceptance rows in canonical
tag order and five further tag orders on per-family controls, before/after
portable restoration and through scalar/mixed/repeated/tail layouts.
All five assets have
been regenerated and validated with the current patched dependencies; their
hashes are recorded in `assets/evaluators/README.md`.

The completed performance survey uses matched 4,097-point, seven-repetition
workloads against unchanged Fortran compiled with `-O2 -fno-fast-math`. Python
uses batches of 1,024 plus a one-row tail; Rust and Symbol routes are measured
separately. The requested 1.5× target is **not met across all workloads**.
See [the per-family report](performance/2026-09-08/README.md) for raw data,
correctness gates, conversion overhead and shared-host limitations.

A subsequent generic Symbolica binary64 Li2 specialization removes MPFR work
from finite order-2 evaluator callbacks; it does not change arbitrary-precision
or other-order evaluation. Three state-free tests pass, including 1,324
comparisons (662 selected inputs at 128/256 reference bits), with maximum
observed `|error|/(1+|Li2|) = 4.15e-16`, plus explicit endpoint/cut and tiny-width
checks. Twelve focused engine tests and the additional native special-function
registration test also pass with this source. These are sampled numerical and
integration checks, not exhaustive precision guarantees. Integral timings are
measured independently in the performance survey. The complete dependency delta
and reproduction details are recorded
in [patches/README.md](patches/README.md).

Both manifests now disable Symbolica's default features and explicitly retain
`tracing_max_level_info`, `integer-gmp`, `float-mpfr`, `native_code_generation`,
and `bincode`, excluding its optional `faster_alloc` global mimalloc allocator.
The sequential-thread Python crash was in mimalloc during PyO3 string extraction,
before cache cloning or evaluator construction. The system-allocator rebuild
passes the unchanged reproduction and full API suite; no mathematical change or
new initializing-thread lifetime rule is involved. An embedding host must audit
Cargo's unified features, since another dependency can re-enable `faster_alloc`.

The repeated-root real-box limit has been repaired and its formerly ignored
128/256-bit regression now passes on the current engine. A native vacuum
triangle sheet correction and broader permutation regression now pass on the
new pin: 324 evaluations across six orders, three real-sign sets, three widths,
three scales, and 128/256 bits. Analytic symmetry takes precedence over
inconsistent exceptional-point reference values; see
[VACUUM_SYMMETRY_AUDIT.md](VACUUM_SYMMETRY_AUDIT.md). The historical descriptions
below are historical; the separately completed current-engine checks above and
the dated audits take precedence. The width-1e-30 box check at 128 bits asserts
finiteness only; its analytic target comparison is at 256 bits.

## Preserved prior-engine audit baseline

The expression crate is independent of the unchanged Fortran wrapper and
standalone numerical Rust port. Symbolica is its only direct dependency.
It is not yet registered as a community Python extension.

The requested public master Symbols and dedicated Symbolica `EvaluationInfo`
hooks are now implemented, with native tagged map definitions in parallel.
All three interface tests pass, including all 879 acceptance coefficients through
both hooks and native maps. The three-mass box NaN exposed by an earlier run is
repaired. See [MASTER_SYMBOL_AUDIT.md](MASTER_SYMBOL_AUDIT.md).

**Global scalar parity is not established.** The expanded audit originally found
44 fixed-f64 disagreements. The repaired expressions now match all 342 accepted
Fortran points / 1,026 coefficients at 256 bits. The final fixed-f64 rerun retains
14 failures (six B0 and eight dB0 finite small-momentum coefficients); all sampled
A0/C0/D0 and pole coefficients pass. Focused regressions also pass
99 bubble points at 256/384 bits, 20 triangle points at f64/256/384, and 63 box
points at f64/128. Fixed-f64 cancellation, additional exact boundaries and
legacy-output versus analytic-value distinctions remain explicit audit items.
At that prior baseline, the exact repeated-root box
`p=[0,0,0,0,4,-3]`, `m=[1,1,1,1]` had an incorrect imaginary part at 128/256 bits.
The negative-mass vacuum triangle also had a native continuation error alongside
inconsistent original-Fortran mass permutations. The current repairs and
independent analytic targets are described above and in the dated audits.
The dated follow-up in `SCALAR_PARITY_AUDIT.md` takes precedence over the
preserved older acceptance results below.

## Current expression coverage

- A0 and tadpole tensors through rank four.
- B0, dB0 and bubble tensor coefficients through rank four (the longitudinal
  `B1`-`B1111` also as native `NativeTensorEvaluator` families).
- C0 and D0 dispatch for every mass mask, including the existing infrared sectors.
- The negative-Källén C0 representation and the all-zero-momentum massive
  triangle limit.
- Exact finite two-mass box boundary branches (including vanishing r12/r13),
  and the equal-mass, zero-external-momentum box `1/(6(m²)²)`.
- Continued root/log/Li2 definitions shared through a native FunctionMap.
  Scalar-master hooks now use the direct generic Rust transcription; the
  transparent expression bodies and explicit manual evaluator backends remain
  available in parallel. No Fortran runtime is used by these master hooks.

C0/D0 return compact mapped series; OneLoopExpressions shares one cached definition
map across an amplitude. Cyclic mass sectors reuse canonical definitions.
The implementation is split into two-point, triangle, triangle-continuation,
box, shared-definition and sheet-arithmetic modules.

## Preserved prior-engine checks and fixture provenance

The prior-engine 2026-09-08 default rerun passed **44 tests** (30 unit tests and
14 integration tests), including all 879 Fortran acceptance coefficients through
both master hooks and native maps as well as the direct-expression route.
Four diagnostics were ignored by default; the then-separate repeated-root
box diagnostic failed before the current repair. The expanded 342-point audit is
also opt-in and passes at 256 bits. `AUDIT.md` retains the earlier 31-test baseline.
Tests arrange their own stack without an external RUST_MIN_STACK setting.
Symbolica can emit normalization
warnings while constructing singular expressions in inactive native if branches;
the evaluated active coefficients are required to remain finite by the tests.

The checked-in Fortran suite contains 293 points / **879 Laurent coefficients**:
4 A0, 13 B0, 12 dB0, 125 C0, and 139 D0 points. It includes 48 points with nonzero
Laurent poles, every C0/D0 mass mask, scales between 1e-4 and 1e4, Euclidean and
selected timelike/scattering invariants, complex masses, cyclic permutations and
finite zero-boundary regressions. The comparison rejects NaN/infinity and applies
2e-10 + 2e-8*|reference| **after mass-dimension normalization**. The original
acceptance/hook checks use mu² as their characteristic squared scale. The expanded
audit, portable/batched tests and benchmarks instead use the largest absolute
kinematic input component, excluding mu², with the appropriate integral mass
dimension. That latter convention prevents independent scale changes from
relaxing the tolerance and small, dimensionful box values from
passing with large relative errors. The test enforces fixture inventory, row
shape, valid integral kinds, and mass-mask/pole coverage.

Regeneration with tests/generate_fixtures.py reproduced the checked-in table
byte-for-byte. Oracle provenance:

- OneLOopBridge revision: 402bf0d70bfa2a9d8ae66b8251cc4395b7b5ab99.
- Test driver: `tests/support/oracle.f90`, linked to the unchanged original library.
- Oracle executable SHA-256:
  5c4fee3ec2a2742fe92da699a2ed2f0e5954ef765e1d98b0c449157bb8863667.
- Original library archive SHA-256:
  8b429a58d6e9404608d447d3a7e0a542fc74ed2adb545cc71d35ebeafb8c854b.
- Fixed random seed 0x5EED; inputs and expected coefficients are stored in the table.
- Ordinary tests do not require Fortran, Python or the numerical Rust port.

Additional tests cover the existing A/B/tensor/IR references, isolated continuation
regressions, mixed C0+D0 expressions sharing one map, and 324 low-level sheet
products/quotients (with mandatory inventory checks). Those auxiliary sheet references come from the numerical Rust
transcription. Its truncated dilogarithm series has about 8e-8 error in selected
unit-circle cases: auxiliary Li2 comparisons therefore allow 2e-7, while raw
values/sheet indices/logs remain strict. Independently, the native results are
compared against 128-bit evaluation with a 2e-12 bound. This does **not** relax
the independent Fortran integral tolerance.

## Fixes and native evaluation

The wrong reciprocal root convention and expanded real/imaginary projections
could select an incorrect continuation sheet. Projections and sign selection
now use explicit native function definitions with single evaluated arguments.
This preserves their exact meaning without a floating tolerance.

For the preserved prior-engine baseline only, the sibling Symbolica dev checkout
was based on `0b57776bf911faeea7e28ea133706fb03740ffeb` and contained the following
local fixes. The current `fb845d34` dependency and its additional repairs are
described in the current-work section above and in `patches/README.md`.

1. External-function index remapping when constants are lifted from a non-inlined
   evaluator.
2. Native absolute values on the real axis, avoiding a square/square-root round trip.
3. Polylogarithm convergence near the unit circle, including avoiding premature
   termination at the zeta function's structural zeros, and exact endpoint/
   real-axis predicates that do not round arbitrary-precision inputs to f64.
4. Reuse of repeated non-inlined definitions before recursive compilation,
   retaining captures, tags and lifted constants. This prevents exponential
   recompilation of the compact native dependency graph.

These are native Symbolica changes, not replacement numerical integral helpers.
The manifest's local paths must be replaced with an agreed upstream revision or
maintained patch before distribution. The changes are supplied in
patches/symbolica-dev.patch for review.

Large generic evaluator builds need an adequate stack. Tests arrange a 128 MiB
thread for the complete symbolic workload; see README.md for application and
restricted-license threading constraints. There is no hidden worker thread in
the library, no custom adaptive-precision gate, and no new CI.

## What these checks do not establish

Passing the sampled cases is not proof of global Fortran equivalence. In particular,
arbitrary Gram/Cayley degeneracies, all threshold/real-axis limits, arbitrary
negative-real-part masses, and the complete timelike/complex kinematic domain
are not exhaustively covered.

The new one-mass C0 quadratic-degeneration formula includes the regular
coincident limit: at squared invariants `[-4,-1,-1]` and squared masses
`[0,0,1]` it simplifies to `-ln(2)`, instead of the former zero guard. The original
Fortran reports an error and returns zero at this point, so this is an independent
analytic check rather than a successful-Fortran parity fixture. It passes at
f64/256/384, with the last two compared against native `Float::log(2)` to an
absolute error below 1e-60. Generic
degeneracies in other sectors still require separate analytic validation.

The fully scaleless derivative dB0(0,0,0) is undefined in Fortran and excluded
from finite-value fixtures. Direct construction marks it indeterminate; generic
evaluation returns a nonfinite value rather than a spurious zero.

The historical expanded-table checks at 256 bits and focused 128/256/384-bit
checks remain sampled evidence, not a global stability proof. The current pass
adds 512/3456-bit precision-provenance controls and Python arbitrary-precision
checks; their scope and provenance are recorded in
[ARBITRARY_PRECISION_AUDIT.md](ARBITRARY_PRECISION_AUDIT.md) and
[PRECISION.md](PRECISION.md). The former approximately 900-bit threshold cap in
the master-reachable Li2 implementation has been removed. Neither that repair
nor requesting `prec=1000` certifies 1,000 accurate digits for all kinematics.
Independent tensor fixtures and broader higher-precision D0/threshold sweeps
remain useful release work. The final rebuilt Rust and measured performance
results are recorded in [NATIVE_RELEASE_AUDIT.md](NATIVE_RELEASE_AUDIT.md).

The native map now includes the boxc T/T1/T13/planar reductions, contour residues
and exact infinitesimal directions. Their presence does not establish every
kinematic limit: repeated-root thresholds and general Gram/Cayley degeneracies
remain explicit audit items. Higher tensor triangles/boxes and higher Laurent
orders are not promised by this API.

Community integration, dependency/version alignment, final packaging metadata and
the licensing permission being handled by the project owner remain release work.
See COMMUNITY_INTEGRATION.md for the inspected Spenso-style adapter plan.
