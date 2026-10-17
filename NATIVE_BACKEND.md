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
and logarithm lips. They avoid converting arbitrary-precision predicates through
f64. Dilogarithms call the same state-free native Symbolica numerical kernels:
binary64 calls `symbolica::transcendental::dilog_complex_f64` directly (a
Bernoulli-series kernel with Symbolica's lips and signed zeros, about 2.5 ulp,
Symbolica PR89), without the registered-callback lookup or MPFR. DoubleFloat's
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
`log(1-1/x)` between ranks. Validation and timings are in
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
