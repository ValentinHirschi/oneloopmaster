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
  exact test. Sign tests take real values: the continued dilogarithm's
  unit-circle test signed `1 - x conj(x)` itself, and `Float` evaluation can
  leave it a rounding-level imaginary part, which a sign test reports as
  negative. Small arguments were then inverted (a `4 pi²` jump on windings
  beyond `-1, 0, 1`), at 128 to 256 bits in four-mass boxes whose invariants
  all exceed their pseudo-thresholds; the test now signs the real part.
- Exact tests rely on identities that binary64 keeps: `x - x = 0`, `(-1) x
  = -x` and `a b = b a`. Each `Float` carries its own precision, and a
  square root has one bit more than its argument. Symbolica's numerica
  used to round the more precise factor of a product to the other's
  precision first. Evaluators form `a - b` as `a + (-1) b` with a `-1` at
  the working precision. So `x - |x|` (the sign test) and `x - conj(x)` (the
  imaginary part) came out as one-bit residues for about half of all square
  roots, and `z conj(z)` kept an imaginary residue. On these residues, 62
  native and 44 interpreted rows of the 256-bit D0 sweeps took wrong
  branches (by up to 4e6) where binary64 was right. The pinned Symbolica
  (symbolica-dev/symbolica#97) rounds products and quotients once and
  keeps factors of exactly ±1 exact (`tests/float_exact_products.rs`,
  STATUS.md). Values equal only analytically, formed by different
  expressions, can still leave residues: STATUS.md lists the rows where
  they still branch.
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
- The contour box's poles `y1, y2` (`beta_data`) are where the line
  `g x + h y + j` meets the conic along the contour direction `(1, beta)`:
  the roots of `K y² + L y + N = a1 Q(y) - (h y + j) B(y)` with
  `a1 = g + beta h`. At the root `u` of the line on an edge of the triangle,
  `h u + j` is a multiple of `a1` (`-a1`, `j a1/(g + h)`, `j a1/g` on the three
  edges), so the quadratic there is `a1 (Q(u) - ratio B(u))`. Next to a
  vanishing Gram determinant `a1` is small (`a1 a1' = (b g² - c g h + a h²)/b`
  over the two betas vanishes with it): the line is nearly parallel to the
  contour direction and one pole lies next to the root of every line edge
  (within 5e-6 to 7e-5 at the NinjaRS gg → HHH box, `a1/g = 1.4e-5`). The
  `s3` terms of those edges took `y - z` from the rounded pole and root, each
  edge rounding it differently, and their sum lost 1e-12 of itself (3.6e-11
  in the box). Where the nearer pole's `y - z` cancels by more than 64
  (`4096 |y - z|² < |z|²`), `r1` now takes it as Vieta's product
  `(y1 - z)(y2 - z) = a1 (Q(u) - ratio B(u))/(K e²)` over the farther pole's
  offset (`box_complex::line_product`, `r1_line`, `r0_line`), with the poles'
  own rounded `a1` and `K`: the edges then see the line the rounded `a1`
  describes, and a consistent `a1` costs nothing (in a scan box, moving
  `beta` by 1.9e-12 with everything derived from it moves the six terms of a
  256-bit `tfun` by up to 7e-18 and their sum, 8.6e-11, by 5e-24). An exactly
  real `y - z` stays real, so its lip follows the infinitesimal directions. A
  threshold of 16 also fired on the recorded NinjaRS gg → HH stream, without
  an accuracy gain on the scans, for 0.56% more retired instructions per D0
  call there (0.47% with 64).
  `Q(u) - ratio B(u)` is the conic at the line's intersection with the edge,
  so a pole also lies next to an edge's root where the conic nearly vanishes
  there, with an `a1` of the size of its terms. Then the polynomial cancels at
  least as much as `y - z`, and the line's term cancels against the conic's
  term on the same edge, whose root lies next to the same pole; both must
  take `y - z` from the same rounded pole (with the product, three-mass boxes
  lost 1.7e-11 and 3.9e-10 where 09b2c74 kept 5.8e-13 and 2.9e-12). The
  product is therefore taken only where `a1` itself cancels by more than 64
  (`4096 |a1|² < |g|² + |beta h|²`), and not on an edge whose intersection
  lies on the conic to within the square root of the rounding; an `a1` that
  rounds to zero is taken as `2^-52 g`, whose logarithm cancels between the
  edges (exactly zero offsets gave `log 0` on each edge).
  `K` of one beta is `b a1'` of the other: as `a1'` vanishes, one pole of
  this contour goes to infinity and its terms vanish with `K`, by which they
  are divided. A `K` that rounds to zero (0/0: every binary64 route returned
  NaN in 28% to 40% of the 1-ulp perturbations of a box with
  `Gram/(p1² p2² p3²)` at rounding level) is replaced by `-L/(2^52 (1 +
  |N/L|²))`, which puts the far pole 2^52 beyond the near pole `-N/L` and the
  triangle (`box_complex::chart_k`); the terms are then their limit to
  O(2^-52). The other contour, whose own `a1'` is then at the rounding level,
  is singular there: its poles sit at the line edges' roots with lips decided
  by the rounding, and from the same rounded arguments its 256-bit value
  differs from the first contour's by O(1). Where the selected contour's
  `a1` keeps fewer than 24 bits (`2^56 |a1|² < |g|² + |beta h|²`), the
  other's does not and the other has a pole within 2^7 of the triangle,
  `tfun_chart` evaluates both and keeps the selected chart only if the two
  agree to 2^-26; with both poles far (`K` and `L` vanish together) the other
  chart's terms cancel and the selected one is kept. Switching to the other
  contour unconditionally lost accuracy on more rows than it gained (the line
  products keep the selected chart's digits down to a small `a1`). The
  thresholds are fixed powers of two on every binary64 route and in
  `DoubleFloat`, so that all of them take the same chart on binary64 inputs:
  the earlier test that followed the working precision, `1 + x` rounding to
  1, took other charts on the generated, interpreted and SymJIT routes, where
  the norm `z conj(z)` of a composite `z` kept an imaginary residue (35% of
  random vanishing-`a1` samples; v1 rows of the verifier's sweep lost up to
  7e-9), and is never true in `DoubleFloat`, which keeps `1 + x` exactly.
  `Float` takes the selected chart: at 128 bits the other chart lost two boxes
  with exactly vanishing invariant combinations to 4.2e-7 and 1.7e-6 (28 and
  22 tracked bits; 09b2c74 4.5e-9 and 1.1e-16). Where no chart is regular,
  binary64 tests the chart it uses: a flat selected chart must keep 10 bits
  when the line's `g` moves to `g (1 + 2^-52)` (its sums by the same exact
  amount), and a chart whose `K` rounds to zero (or the other chart taken,
  whose `K` is `b a1` of the flat one) 13 bits when `K` is set to its rounding
  level `2^-52 (|b a1| + |h b1|)`; otherwise it returns NaN, which the caller
  detects. Of the 26 boxes of the verifier's 32,600-row sweep that 09b2c74
  returned as NaN and f2adce4 as finite wrong values (up to 2.8e2), 11 are
  NaN again, 14 right and one within 10 times its input conditioning; 7
  boxes that 09b2c74 already got wrong by more than 1e-6 are NaN, none that
  it got right, and 23 boxes that f2adce4 got right (NaN or wrong at
  09b2c74) are NaN again (`tests/data/contour_unresolved.txt`). Nine
  boxes with the Gram determinant at the rounding level, wrong by 1.4e-2 to
  140 times their value or non-finite at 09b2c74 and f2adce4, keep 1.5e-11
  (4.9e-9 one of them; `tests/data/contour_gram_rounding.txt`). The
  NinjaRS box keeps 1.9e-14 to 4.1e-14 on the binary64 routes, and 199 of 200
  random boxes next to a vanishing Gram determinant keep 3.2e-12
  (`tests/d0_vanishing_gram.rs`, references OneLOop 3.6 in quadruple
  precision); OneLOop's double precision loses more than 1e-9 on 16 of them.
  The remaining one (1.3e-11 interpreted to 3.2e-10 hand-written, OneLOop
  5.5e-11) has a nearly double `beta` in one `tfun` and nearly coincident
  poles next to an edge's end in the others, whose quadratic `s3` terms
  cancel by 20.
- A pole next to an end of a contour-box edge puts the arguments of `r1`'s
  dilogarithms (`(y - z)/(-z)` and the like) next to one, where a sheet with
  windings carries `i pi n log(1 - x)`; from the rounded quotient `1 - x` keeps
  only its absolute precision (8% of the logarithm for a pole 7e-16 from
  x = 0). `r1` now takes the complements from the poles (`y/z`, `y/(y - z)`,
  `(1 - y)/(1 - z)`, `(y - 1)/(y - z)`, `sheet_exact::complement_difference`).
  Next to the end x = 1 the pole itself keeps only eps of its offset: with a
  zero external invariant, a massless `m2` and the line through the corner
  (1, 0) of the parameter triangle within an ulp (an invariant within an ulp
  of a squared mass), a pole lay 4e-31 from that corner and the box was wrong
  by up to 130 times its value. `zero_external` then takes the reflected
  labelling, which moves that corner to the origin
  (`tests/d0_contour_regressions.rs`: every relabelling within 7.7e-12),
  unless one of the reflected chart's conics has a double contour direction
  (its discriminant `c² - 4 a b`, a Källén function of three invariants,
  within 2^-26 of its terms): there the two contours nearly coincide, and a
  box with the Gram determinant at the rounding level kept 2.4e-6 reflected
  where it keeps 2.3e-9 unreflected.
- The finite massless, one-mass and two-adjacent-mass boxes (OneLOop's
  `boxf0`, `boxf1` and `boxf2`) take dilogarithms of `X r34/r13` and
  `X r24/r12` for the roots `X` of a quadratic that tends to
  `(r34 X - r13)(r24 X - r12)` as an invariant between two massless lines
  vanishes (`r23`, or `r14` without masses). Next to such a collinear limit
  one argument of each pair is next to one, on sheets whose dilogarithm
  carries `i pi n log(1 - x)`; formed from the rounded argument, `1 - x` lost
  `|1 - x|` relative digits (2.9e-7 at an invariant of 7e-9 and O(1) at
  1e-14, as in OneLOop's double precision). The roots follow `solabc` (the
  larger directly, the other from Vieta's product, or OneLOop's labels
  without masses), and `1 - x` comes from the roots' offsets from `r13/r34`
  and `r12/r24` (`ordered_root_offsets`), where the quadratic is the small
  invariant times a closed form. Above half the dilogarithms take it for
  `log(1 - x)` and `Li2(1 - x)`; it also decides the side of the unit
  circle, which sheets with two or more windings jump across
  (`sheet_exact::complement_difference`). Far from one the values are kept:
  formed from other expressions, the complements need not round as the
  values do. Next to `r13 r24 = r12 r34` the discriminant is small as well
  and the dilogarithm terms grow as its inverse; `b² - 4ac` cancelled by
  `b²/d²` there (3e-5 in a one-mass box), and takes the exact form
  `(r13 r24 - r12 r34)² - T (2 (r13 r24 + r12 r34) - T) + rest`
  (`collinear_square`). The two dilogarithm terms' shifted coefficients are
  `∓(r13 r24 - r12 r34)` up to collinear terms; summed in the evaluator's
  own order, each rounded that difference differently and their `1/d` parts
  no longer cancelled (1e-11 in a two-mass box at 2e-3 from it, which
  OneLOop's double precision keeps to 4e-13), so it is formed once
  (`channel_gap`). Binary64 keeps 2e-13 down to invariants of 1e-14
  (`tests/d0_collinear_limits.rs`); next to both limits the rounding of
  `r13 r24 - r12 r34` remains (6 of 300 such one- and two-adjacent-mass boxes
  beyond 1e-12, up to 1.4e-8; 61 in OneLOop's double precision).
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
