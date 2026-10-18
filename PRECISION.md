# Arbitrary-precision evaluation

The same exact Symbolica expression graph is available as a native expression
interpreter and as generated generic Rust arithmetic. `NativeEvaluator<T>`
supports `f64`, `DoubleFloat` and `Float`; `PrecisionEvaluator` (Rust) and the
Python `prec` keyword provide decimal-precision `Complex<Float>` interfaces.
The direct native route is a transcription of that graph, not an alternative
numerical OneLOop algorithm, a Fortran call, or a high-precision wrapper around
binary64 output. The core uses Symbolica's numeric types.

See [the native evaluator](src/native/mod.rs),
[numeric primitives](src/native/primitives.rs),
[precision adapter](src/precision.rs) and [source generator](examples/generate_native.rs).
`get_expression` exposes the fully expanded exact Laurent coefficients with
explicit resource limits; see [expression inspection](README.md#combining-integrals).

## Contract

- `prec`, `PrecisionEvaluator::new` and `PrecisionEvaluator::with_backend`
  specify positive **decimal digits**.
  Python defaults to 16. The Rust binary-precision constructor instead takes
  bits and does not add guards.
- Decimal requests use a conservative decimal-to-binary conversion plus 32
  fixed guard bits. This is a fixed-precision
  calculation, not an adaptive error estimate or a promise that all requested
  digits survive ill-conditioning. Working bits are exposed by
  `PrecisionEvaluator::binary_precision()`.
- Inputs and exact rational expression coefficients are converted directly to
  the working Float precision. Python Decimal and large integers are parsed
  from exact text, including the squared renormalization scale. Ordinary Python
  floats retain their already-rounded binary64 values.
- Raw `NativeEvaluator<f64>` uses 53 bits; `NativeEvaluator<DoubleFloat>` uses
  106 bits. `NativeEvaluator<Float>` defaults to 128 bits and accepts explicit
  positive working bits. Fixed-type constructors reject a different bit count.
  DoubleFloat elementary arithmetic retains its compensated components; only
  Li2 converts both components directly to Float at 160 bits and back. It never
  narrows the whole computation to f64.
- Output precision is not artificially increased after evaluation. Python
  returns DecimalComplex components with at most the requested decimal digits
  and no more than the computed Float precision supports. Python Decimal
  construction from text does not round to the ambient decimal context.
- Real external squared invariants, nonpositive imaginary squared masses,
  positive real `mu_squared`, and finite input components remain required.
  Domain checks use the supplied high-precision values, not f64 approximations.
- Rust supports scalar and flat row-major batches; Python supports scalar and
  list-of-rows batches. Each row returns finite, simple-pole, and double-pole
  coefficients. Arbitrary-precision batches evaluate scalar rows through the
  selected Native or Expression backend, not SIMD. Shape/domain validation is
  completed before any batch output is written.
- Public master Symbols always use direct native callbacks for `Complex<f64>`,
  `Complex<DoubleFloat>` and `Complex<Float>`, independently of the manual
  default. Their Float caches are keyed by the maximum component bit precision
  supplied to the callback, including `mu_squared`. Transparent function-map compilation
  remains available and avoids making the master opaque to an outer optimizer.

- Two expression branches depend on the working precision. The dilogarithm
  divided difference of nearly coincident arguments (`sheet_exact::near_difference`)
  uses a Bernoulli expansion through `s^22` only when its first omitted term,
  bounded by `8 (max|s|²/39)^12`, leaves `1` unchanged at the working precision
  (`(1 + bound) - 1 == 0`). That holds throughout in binary64, but at 256 bits
  only for `|s| < 3e-3`; otherwise the direct difference is used, whose loss
  stays in the tracked precision. DoubleFloat represents `1 + bound` exactly and
  keeps the direct difference. Likewise the contour box's `r1` takes its
  far-root form, `fy log(-z) - sum_{n<=8} M_n/(n z^n)`, only when
  `16 (max(1, |y1|², |y2|²)/|z|²)⁴` leaves 1 unchanged (in binary64 for
  `|z|` beyond about 140 `max(1, |y|)`; at 128 bits and more only for huge
  roots, where the general form costs eps `|z|`). Every other branch is an
  exact test.
- Each quadratic of the contour box's `s3` is the Feynman polynomial of a
  two-point subgraph, `p x² + (mi - mj - p) x + mj`: its discriminant is the
  Källén function `λ(p, mi, mj)` and its values at the ends of the
  integration are the masses. The callers pass `λ` formed from the invariants
  as `(p - (√mi + √mj)²)(p - (√mi - √mj)²)` and the mass at the far end. A
  pseudo-threshold next to a massless line (an external invariant equal to
  its massive neighbour's mass, `λ = (p - mi)²`) is then an exact double root,
  and next to it the roots keep their separation `sqrt(λ)/p` and the signs of
  their slopes `±sqrt(λ)`; a massless end is an exact root 0 or 1. Formed
  from rounded coefficients, the discriminant kept a residue of rounding
  size there whose square root cost binary64 8 to 10 digits; collapsing it
  within rounding cost the separation instead (up to 3e-7 at 5e-7 from the
  threshold). Binary64 now keeps 1e-14 on both sides
  (`tests/d0_contour_thresholds.rs`). Within an ulp of such a threshold the
  second root rounds onto the massless end, so both roots take their slopes
  and continuations from the exact end values. The contour line's coefficient
  sums (`g + h`, `g + j`, `g + h + j`) are formed from the invariants too: a
  vanishing one is an exactly vanishing linear coefficient or an exact root at
  an end (with `p4² = s23` in the zero-external chart the rounded `g + h` was
  1e-16, put a root at 1e15 and made one of two equivalent contours O(1)
  wrong in binary64). A small invariant `a + b + c` also puts a root `beta` of
  the chart's `b beta² + c beta + a` next to 1, where `d + beta e` and the
  like cancel; nearer 1 than 0 the chart coefficients use `beta - 1` from the
  product of the roots (`box_complex::beta_data`).
- `Float` tracks a difference with the bits left after its cancellation and
  computes later operations at that precision, even when the difference is
  exact. Coefficients that are invariants (`m2 - m3 - p2` for equal masses, or
  a coefficient sum equal to a tiny invariant) are therefore formed without
  such cancellation in the contour box (`tests/d0_tiny_invariant.rs`).

`EvaluationBackend::{Native, SymJit, Expression}` selects the manual Rust route.
`DEFAULT_BACKEND` is Native, selected from measured mixed-family throughput.
`PrecisionEvaluator::with_backend(..., EvaluationBackend::Native)` selects direct
native Float; `Expression` selects the interpreter, and `SymJit` is rejected for
Float. Without an explicit selection, `PrecisionEvaluator` defaults
to Native. `evaluate_with_backend` and `evaluate_batch_with_backend` provide
equivalent explicit choices for binary64 convenience calls.

Python `backend="auto"` likewise uses Native for ordinary
float/complex/small-integer calls at `prec=16`, returning Python complex numbers.
Decimal inputs, large integers, or any nondefault precision select the Float
native implementation and DecimalComplex outputs. `backend="native"` explicitly
selects direct native evaluation in either domain; `"expression"` selects the
interpreter, and `"symjit"` rejects AP inputs. In particular, supplying a Decimal
at the default precision preserves
its value until conversion to the requested working precision; it does not
request unlimited precision implicitly.

## Dependency scope

OneLOop pins unmodified Symbolica main at
`821b02451256a92039a0665006628bd5d91470cc` and its bundled Numerica. That revision
includes arbitrary-precision polylogarithm convergence and endpoint checks,
Cartesian complex roots, dominant-axis phase evaluation, precision-preserving
power identities, and direct half-integer powers. Numerical dilogarithms use
the public registered `polylog(2,z)` callback in each numeric domain.

The upstream scaled `hypot` and complex logarithm now construct exact constants
at the stronger component's precision. OneLOop uses these operations directly;
its temporary logarithm and magnitude adapters have been removed. Computed small
corrections retain their uncertainty. Regression coverage checks both a weak
tiny component beside a precise dominant component and an accurately known tiny
component beside a low-precision exact unit value. The earlier
[upstream Numerica patch](patches/numerica-hypot-precision.md) is retained as a
historical artifact and is not applied by any OneLOop manifest.

These checks cover arithmetic used by the scalar expressions; they do not
establish thousand-digit accuracy for every Symbolica special function or every
kinematic point. Portable JIT caches remain binary64. See the
[current migration validation](SYMBOLICA_3_MIGRATION.md) and the historical
[precision audit](ARBITRARY_PRECISION_AUDIT.md) for their respective evidence.

## Usage

For a native Float workspace with exactly the supplied working bits:

```rust,ignore
use oneloop::{NativeEvaluator, ScalarIntegral};
use symbolica::domains::float::Float;

let mut evaluator = NativeEvaluator::<Float>::with_binary_precision(
    ScalarIntegral::C0, 3456)?;
// Supply Complex<Float> arguments directly; evaluate / evaluate_batch use
// finite, simple-pole, double-pole output order, with mu_squared last in each row.
```

See the decimal-precision [Rust example](README.md#stack-and-precision) and the
[Python examples](python/README.md#arbitrary-precision), including Decimal complex
masses, scales, batches, and per-call precision/backend changes. Symbol access,
`initialize()`, Python import, and higher-level convenience construction eagerly
prepare all five binary64 Native and all five SymJIT backends. Float workspaces
are precision-specific and built on request. That eager symbolic work requires
an adequate calling-thread stack and the installed Symbolica license's thread
limits. Native evaluation executes generated Rust arithmetic; its first
dilogarithm call can initialize Symbolica state while resolving the registered
polylog callback. The adapter does not spawn a hidden worker.

The prior performance survey measures binary64 evaluation only. Its observed
slowdowns are not upper bounds for arbitrary precision. Likewise, the retained
Fortran reference tables are double-precision references: matching them cannot
establish 1000-digit accuracy. High-precision identity and convergence tests
must be considered separately from sampled Fortran parity.

## Precision checks (2026-09-08)

The following results describe the earlier patched snapshot. Current results
are recorded in the [migration report](SYMBOLICA_3_MIGRATION.md).

The final standalone Python suite passes **27 tests** (83.755 seconds), including
the Native default, explicit Native/SymJIT/Expression selection and rebuilds,
all five families at 32/1000 decimal digits, Decimal elementary vacuum controls,
context28 independence, scales from `1e-1000`, a 1050-digit tiny width,
1024-row machine/AP32 batches, exact large integers, invalid inputs and
helper-free expression inspection. Community integration passes a compile-only
check; runtime inside the community host is not established.

The rebuilt state-free native primitive gate passes **11 tests** against real
Symbolica, covering all three numeric types, signed axes, subnormal widths,
extreme magnitudes, compensated low components, unbalanced Float component
precision, and 1000-digit Li2/root controls. Three independent simple-sector
tests additionally cover 468 parameter combinations for binary64 and Float,
including signs, mass permutations, scales, and exact off-sector rejection.

Separate dependency gates include five precision-provenance regressions at
512/3456 bits, four independent MPC root regressions (including a 324-case
component-precision/magnitude matrix), sixteen existing Numerica Float tests,
and exact zero-arithmetic tests. Four state-free Symbolica Li2 checks at 3456
bits cover Catalan/reflection/inversion identities, refinement, tiny components,
both cut lips, exact predicates, checked convergence limits and reflected-series
coefficients. The binary64 Li2 grid retains a maximum scaled error near
`4.15e-16`. See [the independent audit](ARBITRARY_PRECISION_AUDIT.md) for the
distinction between direct-source, dependency-only and integration evidence.

The final rebuilt all-target Rust suite passes **114 tests**, including all five
complex native masters at 3456 bits against 1120-decimal-digit Expression
evaluation, independently per coefficient component to better than `1e-1000`.
The [release audit](NATIVE_RELEASE_AUDIT.md) records the complete test inventory,
six ignored entries and separately failing binary64 stress test. Cross-backend agreement of the
same formulas is not an independent integral oracle or a guarantee for every
kinematic point.

The repaired high-precision Li2 uses a reflected-coefficient recurrence and
MPFR's integer-zeta specialization. This avoids both premature stopping and the
very expensive general negative-zeta coefficients.
