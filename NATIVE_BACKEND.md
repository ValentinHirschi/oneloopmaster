# Direct Rust backend

This lane adds ahead-of-time Rust functions for the same scalar masters. It is
separate from the existing native-expression/SymJIT backend. Direct Rust is now
the manual default and the sole numerical backend of bare master-Symbol hooks.
Manual `backend="symjit"` and `backend="expression"` remain available; native
compilation is not a promise of a speedup for every input or batch.

## One formula source

`examples/generate_native.rs` exports the exact repaired expressions and their
non-inlined definitions, then emits generic Rust functions under `src/native/`.
Generated functions use ordinary arithmetic, named locals, function calls and
lazy Rust `if` blocks. There is no instruction interpreter, expression builder,
function-map lookup, or JIT compiler in a warmed native scalar evaluation.

The generator tracks branch-local definitions and live outputs explicitly. It
shares repeated arithmetic and identical lazy conditional results only where
the earlier result dominates the later use. Conditional bodies are not hoisted
or evaluated speculatively. Large branches are factored into ordinary Rust
functions; there are no per-point hash tables or runtime expression caches.
Executable generator regressions cover captures, tuple ordering, branch-local
aliasing and an inactive division by zero.

Small exact sectors are written directly in [simple.rs](src/native/simple.rs):
scaleless cases, massless one/two-scale triangles, and the all-on-shell massless
box. They retain separate continued logarithms. For example, with
`L_s = log((-s-i0)/mu_squared)`, the normalized one-scale massless triangle is
`[L_s²/(2s), -L_s/s, 1/s]`; the massless on-shell box is
`[(2 L_s L_t - pi²)/(s t), -2(L_s+L_t)/(s t), 4/(s t)]`.
Exact zero dispatch is used, not a small-parameter cutoff.

## Hand-written binary64 sector kernels

For `f64`, the dominant massive sectors use hand-written kernels in
[`src/native/handwritten/`](src/native/handwritten/) instead of the generated
code: the three-mass triangle (`triangle.rs`'s `trif3` representation, nonzero
momenta with a nonnegative Källén function), and for three or four massive
propagators the contour box of `box_complex.rs` (OneLOop `boxc`) and the
four-mass Denner-Nierste-Scharf box (`boxf4`), with the same channel exchange
and contour selection as `box_series`. Each is a branch-for-branch transcription
of the expressions: the same exact zero tests (`x - |x|` sign tests included),
infinitesimal directions, root labels and continued sheets (`sheet.rs` mirrors
`sheet_exact.rs`). They evaluate one sector dispatch per call and never an
inactive branch, which makes them two to five times faster than the generated
code. Like the expression, the triangle takes the quadratic root that would
cancel in `(b ± d)/(2a)` from the product of the roots (OneLOop's `solabc`), so
that an exact zero produced by roundoff cannot reach a logarithm. Every other sector,
and every sector for `DoubleFloat` and `Float` (whose tracked precision follows
the expression's own grouping), uses the generated code. Regenerating the
expressions does not update these kernels: an expression change in these
sectors must be transcribed as well (as were the continuation fixes checked by
`tests/permutation_regressions.rs`: exact r-function roots at zero momenta and
thresholds, the roots, chart choice and root lips of `trif3`, the root lips of
`boxf4`, and the contour box's double-root split and real residue logarithms). `tests/native_handwritten.rs` checks them
against the generated expressions at 256 bits on randomized kinematics and on
every distinct C0/D0 call of the recorded Ninja benchmark streams
(`tests/data/ninja_master_rows.txt`); `src/native/handwritten/box_contour.rs`
checks each contour-box definition against its registered expression. Timings:
[performance/2026-10-17-claude-c0d0](performance/2026-10-17-claude-c0d0/README.md).

This keeps the expression implementation as the mathematical source. The older
dependency-free numerical port is not copied into production: its finite-width
regulators and then-unrepaired sectors would undo repairs already made here.
The original Fortran wrapper and that separate numerical port remain untouched.

Exact rational constants, pi and fixed-argument dilogarithms are described in
generated initialization data and converted once at working precision. The
Symbolica instruction export includes symbolic constant metadata: its numeric
zero placeholders must never be mistaken for exact zeros by another backend.

## Numeric domains

The generated functions are generic over `NativeFloat`, implemented for `f64`,
Symbolica `DoubleFloat`, and arbitrary-precision `Float`. The core still has only
Symbolica as a direct dependency. Complex masses use `Complex<T>`.

The native primitives preserve exact complex zero tests and principal square-root
and logarithm lips. For binary64 off the axes and inside `[2^-500, 2^500]`, complex
reciprocals use one division and norms `sqrt(x²+y²)` (within a few ulps; tested
against 256 bits); axes, non-finite values and extreme scales keep the scaled
algorithms. They avoid converting arbitrary-precision predicates through
f64. Dilogarithms call the same state-free native Symbolica numerical kernels:
binary64 calls `symbolica::transcendental::dilog_complex_f64` directly (a
Bernoulli-series kernel with Symbolica's lips and signed zeros, Symbolica PR89),
without the registered-callback lookup or MPFR. Its normwise relative error is at
most 1.2e-15 (about 10 units of 2^-53, reached next to Re z = 1/2 where the
reflection cancels; 99.86% of points within 4 units) on a 200,745-point review
grid against mpmath, where the previous MPFR route was correctly rounded. DoubleFloat's
dilogarithm uses a precision-preserving MPFR conversion of both
compensated components, while ordinary generated arithmetic remains DoubleFloat.
This is not a runtime scalar-master evaluator callback.

## Tensor coefficients

[`tensor.rs`](src/native/tensor.rs) is a hand-written generic port of OneLOop's
bubble tensor routines (`bub0`-`bub111`, `bnlog`, `solabc`) and `tadpn`, exposed
as `NativeTensorEvaluator<T>` with `TensorIntegral::{An, B11, B111}`. It is not
generated from the symbolic expressions: their Passarino-Veltman forms cancel
like `(m²/p²)^rank` at small momenta, while OneLOop's root-pair representation
with the logarithmic `bnlog` series stays accurate there. It keeps OneLOop's
mass ordering and reconstruction, replaces its `10*EPSN` degeneracy tolerances by
exact zero tests and its `IEPS` regulator by explicit lips, and shares
`log(1-1/x)` between ranks. A massless propagator's unit root `x = 1` is formed
exactly (`(p²-m²)/p²` is divided, not scaled by `1/p²`). `f64` evaluators run
a binary64 specialization of the same algorithm
([`tensor/binary64.rs`](src/native/tensor/binary64.rs)): real kinematics in
real arithmetic, one `bnlog` (`2 Re`) for the complex-conjugate roots of real
masses below threshold, and the `bnlog` series by Estrin's scheme on a fixed
8, 16 or 24 coefficients instead of Horner's rule; `DoubleFloat` and `Float`
keep the generic code. `tests/native_tensor_binary64.rs` compares it with the
generic binary64 code and 320 bits on 8,000 random rows over every sector, with
the symbolic `b0`-`b111`/`an` expressions at 256 bits, and with OneLOop on every
bubble/tadpole call of Ninja's recorded streams. Besides OneLOop's
values and Ninja's B111 formulas, `tests/native_tensor.rs` checks all four
bubble coefficients against direct Feynman-parameter integrals (595 rows,
`tests/generate_tensor_feynman_fixtures.py`, mpmath): 320-bit results agree
to 3e-27. Validation and timings are in
[the 2026-10-17 report](performance/2026-10-17-claude/README.md).

## Verification and performance gates

The completed gates compare all three Laurent coefficients independently:

- Original acceptance fixtures for all five families, scalar and mixed batches
  including partial tails; explicit diagnostics for known fixed-f64 failures.
- DoubleFloat against original controls and higher-precision primitives.
- Genuine complex C0/D0 at 1000 digits against higher-precision expression
  evaluation, plus independent analytic limits and tiny-width controls.
- Public Symbol hooks, Python Decimal boundaries, and explicit manual SymJIT
  evaluation checked as separate routes.
- Matched warmed Rust/native, Rust/SymJIT, hook, batched Python, and unchanged
  Fortran measurements using the existing survey protocol. The target remains
  at most 1.5 times Fortran; measured misses must be reported explicitly.

The prior binary64 report remains a historical measurement of the SymJIT
implementation at its recorded commit. It is not evidence for this new backend.

The final release suite passes **114 Rust tests**; the standalone extension
passes **27 Python tests**. The expanded 342-point table passes at native Float
256 bits, while the explicit binary64 diagnostic retains 14 B0/dB0 failures.
See the [release audit](NATIVE_RELEASE_AUDIT.md),
[precision audit](ARBITRARY_PRECISION_AUDIT.md), and
[matched timing report](performance/2026-09-08-native/README.md) for the full
inventory, retained failures and per-workload performance limits.

### Baseline correctness checks (2026-09-08)

The first generated implementation passed the 293-point acceptance set for
all five families with `f64` (scalar and mixed batches 1/3/4/5/31/256/1024) and
`DoubleFloat`. A separate test passed genuine complex inputs in all five
families at 3456 working bits against 1120-decimal-digit expression evaluation,
requiring relative error below `1e-1000` independently for real and imaginary
components. These are sampled checks, not exhaustive analytic-region coverage.

The unoptimized test times were 0.76 seconds, 2.44 seconds and 152.32 seconds
respectively, including the tests' different setup costs. They are **not**
matched Fortran performance measurements. The final factored implementation's
matched measurements are recorded separately in the report linked above.
