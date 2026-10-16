# Massive spacelike C0 derivative controls

This audit covers the public `OneLoopExpressions::dual_evaluator` constructor at
four fixed real kinematic points. It adds regression evidence; it does not enable
a downstream derivative provider or certify all massive spacelike triangles.
No production formulas, continuations, numerical tolerances or domain gates are
changed.

## Tested points and components

Arguments are `[p1²,p2²,(p1+p2)²,m1²,m2²,m3²,mu²]`. All masses are positive and
all invariants are strictly negative. The two invariant triples are
`[-1,-5,-4]` and `[-1,-2,-9]`, with Källén discriminants `-16` and `28`.
Each uses both `[2,3,5]` and `[3,3,3]` for the squared masses and scale `1`.
This exercises the 't Hooft–Veltman and ordered three-mass formula branches.
Neither point has a zero Gram determinant or a selected projective singularity.

The seven-component dynamic Dual shape is
`[[0,0],[1,0],[0,1],[2,0],[1,1],[3,0],[2,1]]`.
Components are Taylor coefficients, including factorial normalization.
All finite, single-pole and double-pole components are checked. Pole components
vanish for this massive finite sector.

- Nonlinear simultaneous scaling of every mass and invariant has arbitrary
  complex first, second, third and mixed input coefficients. Symbolica
  differentiates `C0(center)/scale(t,u)` independently, using the homogeneity
  `C0(s*p,s*m)=C0(p,m)/s`. The Float384 comparison uses a scaled `2^-220` bound.
- Independent invariant/mass and mass/invariant directions carry seeds `i` and
  `2-i`. A seven-node stencil is obtained with Symbolica's exact rational
  `Matrix::solve`; its tensor product estimates mixed Taylor coefficients.
  Steps `2^-12`, `2^-18`, `2^-24` must improve every coefficient's error by at
  least a factor of 100 per refinement and finish within scaled `2^-80`.
- The finite integral is independent of the renormalization scale, including
  pure higher/mixed scale derivatives.
- Empty, mixed-branch, permuted and three-plus-one tail batches reproduce scalar
  coefficient evaluations within scaled `2^-240`.
- f64 evaluations agree with the Float384 components within scaled `2^-30`
  (slightly tighter than `10^-9`).
- A separate univariate control exercises degree-five complex scaling, with
  independent Symbolica derivatives of `C0(center)/(1+(1+2i)*t)`.

The tests are in [tests/c0_derivatives.rs](tests/c0_derivatives.rs). Symbolica
provides all differentiation, exact linear algebra and numeric arithmetic.
The native matrix implementation and its example/tests were checked; no reusable
finite-difference stencil generator was found in Symbolica or FeynKit. This test
does not introduce a local linear solver or symbolic manipulation routine.

## Why complex seeds use a real-axis difference oracle

At these centers the Feynman parameter denominator is strictly positive on the
integration simplex: positive masses and negative invariants prevent a pinch.
The integral is therefore holomorphic in a neighborhood of each real center.
A Taylor coefficient in independently phased directions follows from the
real-direction coefficient multiplied by the corresponding powers of the
complex seeds. This permits an independent complex-seed check using accurate
scalar evaluations with real invariants and masses.

This distinction matters: direct scalar evaluation at complex external
invariants has a continuation limitation. At
`[-1,-5,-4,2,3,5,1]`, the finite scalar C0 is approximately
`-0.12562502313360096`. Perturbing only the first invariant by
`-i*2^-24` approaches this value, with imaginary part about `1.9963e-10`.
The perturbation `+i*2^-24` instead gives approximately
`-0.3475526120194918 - 0.8154728128494999 i`.
The positive-Källén control also has a discontinuity for this positive imaginary
perturbation. Direct complex-invariant scalar differences are consequently not
used as the derivative oracle, and these controls make no claim that such
complex scalar inputs are supported. The discontinuity requires a separate
continuation audit and fix.

## Reproduction and limitations

The implementation was checked on oneloopmaster `69761cd` with the compatible
Symbolica/Numerica/Graphica integration
`58b72d438e7bcfce158d132612db7d35fdefa4ba`. This includes the shared Dual vectorizer,
branch-scope correction, native polylog-at-zero series, and bounded literal-pool
follow-up. Patch all three packages together when overriding the source.

```sh
cargo test --test c0_derivatives -- --test-threads=1
```

Both tests pass together with `--test-threads=1` (2 passed, 0 failed).
All 48 independent directional coefficient comparisons satisfy the convergence
requirements. The largest final absolute difference in the recorded Float384
run was approximately `5.34e-33`; the smallest was `6.33e-54`. The separate
degree-five scaling control also passes at all four points. These are
correctness measurements, not performance acceptance evidence.

Zero-Gram surfaces, thresholds, zero masses/invariants, projective degenerations,
wound or mismatched sheets, arbitrary complex primal inputs, and D0 derivatives
remain outside this validation. A generic scalar limit formula does not establish
its transverse derivative neighborhood. The tests do not substitute for those
neighborhoods or a mathematical classification of the full production domain.
