# Lazy zero sheet corrections at a regular C0 point

The massive Euclidean point `C0([-1,-2,-3], [2,3,5], 7)` is finite, but the
negative-Källén representation evaluated `0 * log((y-1)/y)` with `y=0` in two
of its three `s3` terms. Both the exact interpreter and the generated native
implementation returned NaN in f64, DoubleFloat, and Float256 before the fix.
The Källén discriminant is −8, and the Feynman-parameter denominator is strictly
positive; this is a removable intermediate limit, not an integral singularity.

The change wraps the existing correction term in the existing lazy
`if_nonzero` helper. The correction contains only discrete sheet predicates
and multiples of i*pi. Within the selected analytic neighborhood its derivative
components are zero, so a zero primal correction cannot hide a nonzero
continuous derivative. No integral formula, denominator, branch condition,
precision, or tolerance is changed. The existing continued dilogarithm
neighborhood at a unit sheet argument supplies the necessary derivatives.
This does not establish global derivative support at continuation boundaries.

## Independent controls

`tests/generate_c0_projective_fixtures.py` integrates the two-dimensional
Feynman-parameter representation with mpmath. It compares 60- and 90-digit
quadrature and writes exact rational inputs plus 85-digit finite values to
`tests/data/c0_projective.tsv`. The central value is

```
-0.135212396225436074537371827943492306170596769774993734247205...
```

The existing independent three-mass representation also yields this value in
all three cyclic routings and all three precisions. The retained quadrature
fixtures cover the center and ±1/64 perturbations of each invariant and mass.
All 13 points are tested under all three cyclic routings: 39 rows per carrier.

`tests/c0_projective.rs` checks every Laurent coefficient with the exact
interpreter and generated native scalar evaluator, warmed forward/reverse
batches of widths 1, 4, and 16 including tails and empty input, f64, DoubleFloat,
and Float256. Absolute criteria are 1e-12, 1e-27, and 1e-65 respectively.
The same finite values and zero poles are checked with the JIT evaluator.
The unchanged source fails this test on its first central interpreter value.

The new `c0_zero_hv_projective_parameter_higher_mixed_controls` test checks
all six invariant/mass directions in complex-seeded pairs. It validates first,
second, third and mixed Taylor components against converged high-precision
finite differences using the existing Symbolica-generated stencil weights.
The central components also meet the existing f64/Float384 comparison gate.
Float384 neighborhood batches retain every component under lane permutations
and tails. The existing nonzero-correction continuation fixture now additionally
verifies that all three correction factors equal −2*pi*i, ensuring the active
branch remains covered by its independent numerical values.

## Explicit remaining binary64 derivative limitation

At the nearby point `C0([-1,-2,-3], [129/64,3,5], 7)`, with seeds
`d(m0_squared)/dt = i` and `d(p2)/du = 2-i`, the `[3,0]` Taylor coefficient has
an f64 error of approximately 1.0968975723e-8 against Float384. The original
unchanged implementation produces exactly the same result and failure, so
this conditioning issue is not introduced or resolved by the lazy guard.
The criterion remains unchanged (approximately 9.32e-10 at this point).

`c0_near_projective_f64_derivative_conditioning` retains an assertion of that
criterion as an explicitly ignored known-failure diagnostic. Running it with
`--ignored` fails on both the original and fixed implementations. No full
binary64 higher-derivative neighborhood accuracy claim is made; stable
alternative-expression selection or precision rescue is separate work.

## Reproduction and generated source

Validation uses oneloop base e9d14f7b5fa546af33fced2a421c616220178643 and Symbolica,
Numerica, and Graphica integration fd7bced5307e94e994d5deba63c60b8380aa7d02,
with the usual 128 MiB Rust thread stack and serial tests:

```sh
RUST_MIN_STACK=134217728 cargo test --features generated-evaluators \
  --lib --test c0_projective --test c0_derivatives --test native \
  --test native_contract --test evaluators -- --test-threads=1
RUST_MIN_STACK=134217728 cargo test --features generated-evaluators \
  --test c0_derivatives c0_near_projective_f64_derivative_conditioning \
  -- --ignored --exact --test-threads=1
python tests/generate_c0_projective_fixtures.py
cargo run --features rust-codegen --example generate_native -- src/native
```

A fresh unmodified-source generation repeats byte-identically in two separate
processes, but differs from the formerly committed files through mechanical
constant-pool ordering and identifier changes. The baseline refresh is kept in
its own commit. Comparing that freshly generated baseline to the fixed output
changes only the generated C0 helper by adding its lazy branch; the other five
generated files are byte-identical. No generated arithmetic is edited by hand.

## Validation checkpoint

The selected fixed-source suite passes 82 tests, with three explicitly ignored
diagnostics (two pre-existing and the retained nearby derivative limitation):
67 library tests, one independent projective regression, three C0 derivative
tests, five native backend tests, four native contract tests, and two JIT/cache
and original-fixture tests. The JIT suite includes all 293 original scalar
fixtures and six batch layouts; native controls include DoubleFloat, Float256,
and 1000-digit comparisons. The active-correction continuation tests also pass
as a focused rerun. Formatting, Python syntax, diff whitespace, and targeted
Clippy pass; Clippy retains the repository's existing allowances for
`chunks_exact_to_as_chunks` and `clone_on_copy`.

The new projective regression fails on unchanged e9d14f7 with NaN at the first
central value. Explicit execution of the ignored nearby derivative diagnostic
fails on both old and fixed source with identical printed actual/expected
components and errors. Those failures are retained as evidence, not counted as
acceptance or silently converted into successes.
