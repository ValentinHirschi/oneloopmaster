# Native OneLOop expressions

An independent Rust library built on Symbolica. It does not link the Fortran
wrapper or the standalone numerical Rust implementation.
All integral algorithms are exact symbolic expressions; continuation definitions
are transparent `FunctionMap` entries. The same expression graph also generates
inspectable generic Rust arithmetic for direct numerical evaluation. Public
master-Symbol hooks interpret those expressions by default, or use the generated
Rust backend when `generated-evaluators` is enabled.

Development status: the expanded audit passes at 256 bits, but retains 14
fixed-binary64 B0/dB0 disagreements. The measured 1.5×-Fortran performance target
is also not met across all workloads. See the [native release audit](NATIVE_RELEASE_AUDIT.md)
and [performance report](performance/2026-09-08-native/README.md). Full scalar parity is not
established; this is not publication-ready. The [2026-10-17 report](performance/2026-10-17-claude/README.md)
measures the current binary64 backend: all mixed scalar workloads, the gg → HHH
top-loop pentagon (1.15–1.24× per point), the gg → HH box and the H → gg triangle
are within 1.3× of OneLOop Fortran; the γγ → γγ top box is not (1.6–1.7× per
diagram: its massless-leg D0 about 2.1×, C0 1.8×), and single real timelike D0
rows remain about 2.2×. The four-mass real-mass D0 rows with `s12 · s23 < 0`
that had a wrong real part are fixed in the expressions
(`tests/native_d0_four_mass_defect.rs`), together with further continuation
defects found by a OneLOop sweep (`tests/oracle_sweep.py`,
`tests/permutation_regressions.rs`). Two of them are OneLOop 3.6's own: the
sign of its zero-momentum `B1111` when the second mass is the heavier one, and
real-mass contour boxes whose value differs from OneLOop's vanishing-width
limit (`tests/data/real_mass_limits.txt`); the expressions follow the limit.
OneLOop's one-mass box with a complex mass (`boxf1`) is wrong too, at any
width, in 131 of 3000 random such boxes (by 0.3 to 7 of the value; one row
of the regression data by a factor of 6400): it takes the logarithms of
roots normalized by the complex `sqrt(m²)`, which can turn them across the
cut. The expressions take Denner and Dittmaier's sheets there (arXiv:1005.2076,
eq. (3.82)), checked against COLLIER in quadruple precision and a direct
Feynman-parameter integral (`tests/d0_one_mass_width.rs`,
`tests/data/one_mass_width.txt`). Its two-mass boxes with complex masses
(`boxf2`, and `boxf5` for opposite masses) are wrong in the same way, in
0.4% to 2.8% of random such boxes depending on how they are sampled
(STATUS.md; by 0.6 to 17 of the value; next to a
collinear limit by up to 9e4), also with one real mass and with |Im m²| down to
1e-9 of Re m²; the expressions take Denner and Dittmaier's two-mass sheets
there (eqs. (3.78)-(3.80), `tests/d0_two_mass_width.rs`,
`tests/data/two_mass_width.txt`). Where both massive lines sit at or next to
their shells, two of those dilogarithm arguments approach one; their
complements come from the channels' complements, since formed from the
rounded arguments they lost up to 0.16 of the value in binary64
(`tests/data/two_mass_shell.txt`). Next to a collinear limit the one- and
two-mass terms take the root nearer the branch cut across it, and a small
mass-pair complement from the roots' offsets (one mass) or, in binary64,
from the product (two masses): the boxes of
`tests/data/two_mass_conditioning.txt`, which lost 1.2e-11 to 1.6e-9 in
binary64, keep at most 5.7e-14 in every relabelling, and the one-mass boxes
of `tests/d0_one_mass_width.rs` 1e-12.
Dilogarithm divided differences of nearly coincident or small arguments are
evaluated term by term, as OneLOop's `dilog2_c` (a two-mass C0 with nearly
coincident roots was off by 2.4e-13 in binary64, now 2e-15 to 9e-15 on every
route: `tests/c0_coincident_roots.rs`), and the contour box no longer loses a
tiny invariant to cancellation in `Float` (`tests/d0_tiny_invariant.rs`).
The two dilogarithms reach that choice as arguments of a registered
function: Symbolica's evaluator calls a registered function anew at every
occurrence of its call, and with the dilogarithm calls repeated inside the
cancellation test the interpreted D0 and C0 took twice as long at 700 bits
(`tests/expression_construction.rs`).
Massless, one- and two-adjacent-mass boxes next to a collinear limit (a
small invariant between two massless lines) take `1 - x` of their
dilogarithm arguments next to one from the kinematics: binary64 keeps 2e-13
down to invariants of 1e-14, where it lost up to O(1) as OneLOop does
(`tests/d0_collinear_limits.rs`).
Real-mass three-mass boxes use the contour representation wherever it
supports the kinematics: OneLOop's real-mass `boxf33` chart, which it takes
outside timelike configurations, depends on the labelling, with a wrong real
or imaginary part in some relabellings of about a quarter of three-mass boxes,
OneLOop 3.6 included (`tests/d0_real_mass_orderings.rs`). Four-mass boxes keep
OneLOop's selection, but its `boxf4` chart (Denner-Nierste-Scharf), wrong in
about 2% of the relabellings of random real four-mass boxes (OneLOop 3.6
included), is now Denner and Dittmaier's form of it (arXiv:1005.2076, eq. (3.61)):
an invariant below its pseudo-threshold in the `p23` slot, its r-function
exactly real, `x r13` principal with the explicit `eta(x, r13) Log R` term, and
exact infinitesimal decisions (`tests/d0_four_mass_orderings.rs`). The continued
dilogarithm's unit-circle test signed `1 - x conj(x)` without taking its real
part: in `Float` a rounding-level imaginary part made it invert small arguments
(a `4 pi²` jump on higher windings, seen at 128 to 256 bits in four-mass boxes
above all their pseudo-thresholds); it now signs the real part.
The contour box's edge quadratics take their
discriminants as exact Källén functions and their roots at massless ends
exactly, so a pseudo-threshold next to a massless line keeps 1e-14 in
binary64 on both sides instead of losing 8 to 10 digits at it (or 7 next to
it), also within an ulp of it (`tests/d0_contour_thresholds.rs`).
Next to a vanishing Gram determinant (two external momenta nearly parallel,
their invariant next to its threshold, as for two Higgs bosons nearly at rest
relative to each other in gg → HHH) the contour box's line edges take the
offset of a pole from their root from Vieta's product with the poles' own
`g + beta h`: the NinjaRS gg → HHH box is now within 3.4e-14 on every binary64
route (3.6e-11 before), and of 200 random such boxes all but one are within
3.9e-12, where 34 were beyond 5e-12 and 22 beyond 1e-9
(`tests/d0_vanishing_gram.rs`). With the Gram determinant at the rounding
level, a contour whose `K` rounds to zero takes its far pole 2^52 beyond the
near one, and where the selected contour's `g + beta h` keeps fewer than 24
bits the other contour's chart is evaluated too and taken unless the two agree
to 2^-26: nine such boxes, wrong by up to 140 times their value or non-finite
at 09b2c74, keep 1.6e-11 (one of them 2.2e-9, SymJIT 4.9e-9) in every
relabelling (`tests/d0_contour_regressions.rs`). Where no chart is regular,
binary64 returns NaN unless the chart it uses keeps 10 bits when the line moves
by one rounding unit (a flat chart) or 13 bits with its `K` at the rounding
level (a `K` that rounds to zero), so that such a box is non-finite, as at
09b2c74, rather than silently wrong (f2adce4: up to 1.9;
`tests/data/contour_unresolved.txt`); `DoubleFloat` takes binary64's charts
and `Float` the selected one, untested. Where the other chart's `K` is a
rounding residue, the selected chart is kept if it moves 2^10 times less
under its line's move and keeps 20 bits; and the zero-external reflection
is not taken into a chart whose flat contour has no regular alternative
(`tests/data/contour_fresh_sweep.txt`). These hold on the measured sets
(STATUS.md): on the second verifier's 36,840-row sweep, 12 relabellings of
one box still decide differently by route (finite and wrong by 2.9e-2 on the
hand-written route in 8, NaN elsewhere, as wrong at 09b2c74), and 25 binary64
route entries on 12 rows are more than ten times worse than 09b2c74 (none
wrong by more than 1e-6). The
product is taken only where `g + beta h` itself cancels: where the conic
nearly vanishes at the line's intersection with an edge, the line's and the
conic's terms keep `y - z` from the same pole (three-mass boxes lost up to
3.9e-10 otherwise). Next to a corner of the Feynman-parameter triangle, a box
with a zero external invariant and a massless line takes the labelling that
puts the corner at the origin (unless that labelling's chart has a double
contour direction), and the contour box's dilogarithms take their
complements `1 - x` from the poles: a box with an invariant one ulp from a
squared mass, wrong by up to 130 times its value, is within 7.7e-12 in every
relabelling (`tests/d0_contour_regressions.rs`); that labelling is taken only
where a pole of the original chart lies poorly resolved next to the corner.
A soft-divergent box with a massless line whose leg partners are on their
mass shells is found soft only on the external invariants, as in OneLOop:
the numerical D0 entry points (native code in every carrier, SymJIT, the
interpreted expressions and the masters' numeric callbacks) evaluate it in
a labelling of the same box that keeps those on-shell pairs off the
channels, so that all 24 labellings give the poles and the finite part
(`tests/d0_infrared_relabelling.rs`; pyAmpliCol's e+ e- → t t̄ h box with
lines (m, m, 0, m) was finite without its soft pole in 16 of them). Symbolic
uses of the D0 expressions (`OneLoopExpressions::d0`, `crate::d0`,
`d0_mapped` and evaluators built from them) must relabel the point with
`oneloop::regular_d0_labelling` first. [STATUS.md](STATUS.md) lists the
rows on which this candidate is less accurate than 09b2c74, each accepted
and explained, and the verification against COLLIER in quadruple
precision.

This work is based on **Andreas van Hameren's OneLOop package**. The original
algorithms, analytic continuation conventions and Fortran reference results are
his work; see [Credits and references](#credits-and-references).

Inputs are squared invariants, squared masses with `Im(m²) <= 0`, and a positive
real squared renormalization scale. External invariants are real. Coefficients
are ordered `[finite, 1/epsilon, 1/epsilon²]`, with OneLOop normalization. This is
a Laurent expansion through the finite term, not an unexpanded function of epsilon.

## Positive epsilon coefficients

`A0(1, m2, mu2)` and `B0(1, p2, m02, m12, mu2)` evaluate the coefficient
of epsilon. The exact constructors are `a0_epsilon` and `b0_epsilon`;
`ScalarIntegral::laurent_orders()` reports the supported range. Existing
three-coefficient APIs still return `[finite, simple pole, double pole]`.

Both positive orders use the common `1/r_Gamma` normalization, where
`r_Gamma = Gamma(1+eps)*Gamma(1-eps)^2/Gamma(1-2*eps)`. In particular,
A0's positive coefficient is `m2*(1-L+L^2/2+pi^2/6)`, with
`L=log(m2/mu2)`. Replacing `1/r_Gamma` by `exp(gamma_E*eps)` changes
that coefficient by `-pi^2*m2/12`.

B0 uses an exact logarithm/dilogarithm expression for
`pi^2/6 + integral_0^1 log((m02*(1-x)+m12*x-p2*x*(1-x)-i0)/mu2)^2 dx/2`,
with separate scaleless, zero-momentum, massless and coincident-root limits.
The public hooks and transparent function map support binary64 and arbitrary
precision. Tests compare both mass orderings against independently integrated
90-digit parameter-space fixtures, including real-axis branch cuts.

Higher positive orders, and positive orders of dB0/C0/D0, remain unsupported.
Adding these coefficients does not remove the existing numerical limitations
of the finite-order evaluators documented below.

## Development dependency

The direct Symbolica dependency remains versioned, but the Rust and both Python
build roots pin Symbolica, Numerica and Graphica to one tested integration through
`[patch.crates-io]`. The current integration is `9cdc2949` from
[`alphal00p/symbolica`](https://github.com/alphal00p/symbolica), on
`ninjars-pin`, the branch NinjaRS pins its Symbolica from (the same commit is
`claude/oneloopmaster-integration-float-products` of
[`ValentinHirschi/symbolica`](https://github.com/ValentinHirschi/symbolica)).
It is the previous integration `e679e8c1`
(`claude/oneloopmaster-integration-dilog`) plus
single-rounded `Float` products and exact factors of ±1 from
[PR97](https://github.com/symbolica-dev/symbolica/pull/97). Without them,
exact branch tests in arbitrary precision took wrong branches (STATUS.md).
`e679e8c1` is in turn `2c06f245` (`codex/oneloopmaster-integration`) plus the
binary64 dilogarithm kernel of
[PR89](https://github.com/symbolica-dev/symbolica/pull/89). Cargo fetches this
source automatically; no sibling checkout is required for the standalone build.

This integration combines the recursive shared-function Dual support and vector
branch-scoping repairs with the native polylog-at-zero series needed for the
coincident dilogarithm neighborhood. It also includes exact Dual constant-pool
interning, temporary-stack reuse for shared evaluator bodies, and caching of
completed conditionals within their existing branch scope. Conditional reuse
does not reorder arithmetic or evaluate inactive branches. The version
string `3.0.1` alone does not identify these unreleased fixes. The constructor/Dual and Rust-emitter validation
must use this integration or a tested successor containing those changes. See the [combined conditional-cache/C0 validation](patches/conditional-cache-c0-pin-validation.md) for the current checks and the
[shared-helper integration validation](patches/shared-helper-stack-pin-validation.md) for the preceding checkpoint.

An embedding workspace controls its own `[patch.crates-io]` table: patches in a
dependency's manifest do not propagate to the workspace root. Select compatible
Symbolica, Numerica and Graphica sources together to retain one shared kernel.
The [consolidated validation record](patches/consolidated-main-validation.md)
documents exact revisions and test results.

OneLOop uses the public registered `polylog(2,z)` callbacks for numerical
dilogarithms, retains its own exact function definitions for inspection, and
serializes numerical evaluators through Symbolica's public bincode support.
The [migration notes](SYMBOLICA_3_MIGRATION.md) track validation and any remaining
runtime limitations. Historical audits describe the older patched dependency
and do not establish results for this revision.

All three manifests disable Symbolica's default features and explicitly retain
`tracing_max_level_info`, `integer-gmp`, `float-mpfr`,
`native_code_generation`, and `bincode`. This excludes `faster_alloc`, its optional
global mimalloc allocator; the standalone extension uses Rust's system allocator.
The earlier Python sequential-thread crash was in mimalloc during PyO3 argument
extraction, before evaluator construction. The unchanged reproducer passes with
this configuration, without a cache-lifetime or mathematical change.
Cargo features are additive: another dependency enabling Symbolica's defaults
or `faster_alloc` can re-enable mimalloc in an embedding host. Audit the resolved
feature graph at the consuming build root, not just this library's manifest.

## Quick start

Clone this repository as `oneloopmaster`. Cargo fetches the versioned dependencies
and the pinned Symbolica integration described above.
A recent Rust toolchain supporting
edition 2024 and a working Symbolica license/setup are required (the validated
toolchain is Rust 1.98.1). No Fortran compiler is needed for ordinary use or tests.

```sh
git clone https://github.com/alphal00p/oneloopmaster.git
cd oneloopmaster
cargo run --example basic
```

[examples/basic.rs](examples/basic.rs) is a complete runnable example. Its core is:

```rust,ignore
use symbolica::prelude::*;

let p = parse!("p");
let mass_squared = Atom::num(1);
let mu_squared = Atom::num(1);
let bubble = oneloop::b0(&p, &mass_squared, &mass_squared, &mu_squared);
let views = bubble.coefficients().each_ref().map(Atom::as_view);
let mut evaluator = Atom::evaluator_multiple(&views, &[p])
    .build()
    .expect("compile B0")
    .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
let mut output = [Complex::new(0.0, 0.0); 3];
evaluator.evaluate(&[Complex::new(-1.0, 0.0)], &mut output);
println!("[finite, simple pole, double pole] = {output:?}");
```

At `p²=-1`, `m0²=m1²=mu²=1`, this gives approximately
`[-0.1520447048200202, 1, 0]`, with zero imaginary parts. The example also prints
the exact A0 expressions. To evaluate complex masses, make the mass arguments
symbolic parameters and supply `Complex::new(real_part, nonpositive_imag_part)`.
Use complex-valued evaluation even for real inputs: timelike integrals can have
an imaginary part. See the precision and stack requirements below.

The scalar entry points are `a0(m², mu²)`, `b0(p², m0², m1², mu²)`,
`db0(p², m0², m1², mu²)`, `c0([p1²,p2²,p3²], [m0²,m1²,m2²], mu²)`,
and `d0([p1²,p2²,p3²,p4²,s12,s23], [m0²,m1²,m2²,m3²], mu²)`.
Here `s12=(p1+p2)²`, `s23=(p2+p3)²`; ordering follows the original OneLOop API.
`db0` differentiates B0 with respect to `p²`, not the renormalization scale.
Arguments are references to Symbolica `Atom`s. The routines construct the three
Laurent coefficients; they do not numerically integrate Feynman parameters.

## Combining integrals

`a0`, `an`, `b0`, `db0` and bubble tensor coefficients return plain native Laurent
series. `c0` and `d0` return `MappedLaurentSeries`, including the definitions needed
to evaluate them. Their `_mapped` names are aliases. For amplitudes, construct a
single `OneLoopExpressions` context and share its definitions:

```rust,ignore
use oneloop::OneLoopExpressions;
use symbolica::prelude::*;

let context = OneLoopExpressions::new();
let p = [parse!("p1"), parse!("p2"), parse!("p3")];
let m = [parse!("m1"), parse!("m2"), parse!("m3")];
let mu2 = parse!("mu2");
let c = context.c0([&p[0], &p[1], &p[2]], [&m[0], &m[1], &m[2]], &mu2);
let combined = c.coefficients()[0].clone() + oneloop::a0(&m[0], &mu2).coefficients()[0].clone();
let parameters = [p[0].clone(), p[1].clone(), p[2].clone(),
                  m[0].clone(), m[1].clone(), m[2].clone(), mu2];
let exact_evaluator = context.evaluator(&[combined], &parameters)?;
// Map exact coefficients with map_coeff for f64, or map_coeff_with_prec for Float.
```

Keep the function map when exporting expressions or building an evaluator through
Symbolica directly. Compact function calls alone are not self-contained formulas.
The context caches the native definitions once per process; cyclic mass sectors
reuse canonical formulas instead of constructing copies.

`OneLoopExpressions::dual_evaluator(expressions, parameters, dual_shape)` and
`MappedLaurentSeries::dual_evaluator(parameters, dual_shape)` prepare local
neighborhoods before invoking Symbolica's `Dualizer`. They infer the required
total degree from a static Dual or dynamic HyperDual shape, without a fixed
derivative-order limit. The resulting evaluator takes and returns flattened
components in Symbolica's usual order; its coefficients can be mapped to the
desired precision and its instructions can be evaluated in batches.

This currently covers the mapped sheet helper `log(x)/(1-x)` at `x=1` on a
sheet whose continued logarithm vanishes there. Its local polynomial is generated
by Symbolica's `Atom::series`; a nonzero logarithm at one remains a pole.
Equal nonzero sheet quotient arguments retain independent derivative components.
The derivative constructor leaves an unresolved `0/0` quotient nonfinite instead
of applying the historical scalar convention of assigning it one.

Mapped sheet dilogarithm divided differences also retain independently seeded
arguments at coincidence. Symbolica generates the numerator series before
cancellation, including the joint `x=y=1` neighborhood; the same unit series
preserves a single dilogarithm argument's jets at one. At coincidence the two
arguments must carry matching sheet tags. Zero arguments, differing sheet tags,
and wound unit points do not receive an analytic continuation from this
constructor and remain nonfinite. The native polylogarithm series at zero
requires Symbolica [PR84](https://github.com/symbolica-dev/symbolica/pull/84).
The local branches apply at exact coincidence. Nearby points use the ordinary
continuation formulas and retain their cancellation and precision requirements.

Shared-body Dual evaluation requires Symbolica's recursive vectorizer
([PR82](https://github.com/symbolica-dev/symbolica/pull/82)) and its optimizer
scope correction for branches with multiple results
([PR83](https://github.com/symbolica-dev/symbolica/pull/83)). Massive spacelike
C0 controls cover unequal/equal positive masses, both Källén signs, nonlinear
higher/mixed complex seeds, scale independence and mixed-branch batches at four
fixed real kinematic points. Symbolica-derived scaling controls and converged
384-bit real-axis differences are independent oracles. See
[C0_DERIVATIVE_AUDIT.md](C0_DERIVATIVE_AUDIT.md) for the precise scope and a
distinction between unsupported complex scalar invariants and validated complex
derivative seeds at real centers. These tests do not establish a general C0/D0
derivative domain.

Constructing a scalar evaluator and then dualizing it cannot recover derivatives
discarded by a scalar limit. Other limits still need neighborhoods, including
the separate plain-expression log-ratio and divided-difference helpers, and
inspection formulas. Complete C0/D0 derivative-domain coverage is
not established by these constructors; callers must validate their domain.

For a self-contained formula, `oneloop::get_expression(master_call)` (also
`family.get_expression(arguments)` for a Rust `ScalarIntegral`) returns the three exact Laurent coefficients
with every OneLOop helper expanded into Symbolica primitives. It preserves the
supplied symbols' assumptions. `get_expression_with_options` sets explicit node
and depth limits; exceeding a limit returns an error, not a truncated formula.
Expansion can duplicate work and change numerical conditioning: use the compact
function map for routine numerical compilation. See the
[inspection implementation](src/inspection.rs),
[native source generator](examples/generate_native.rs), and
[generated Rust families](src/native/).

### `get_expression` example: a parametric equal-mass B0

For `B0(psq, m², m²; mu²=1)`, keep both the momentum and squared mass symbolic.
The arguments are `[p_squared, mass_0_squared, mass_1_squared, mu_squared]`;
the result is `(finite, simple_pole, double_pole)`, not a numerical evaluation:

Use the [shared-kernel Python build](python/README.md):

```python
from symbolica import E, S, N, Replacement
from symbolica.community.hepkit import oneloop as olo

psq = S("psq", is_real=True)
m2 = S("mz_masses::m2", is_positive=True)  # m², not m
B0 = olo.B0
one_master = B0(psq, m2, m2, 1)  # mu²=1; no Laurent tag here
finite, pole, double_pole = olo.get_expression(one_master)
print(finite)       # Complete expression, including all analytic branches.
print(pole)
print(double_pole)  # 0
finite_only = olo.get_expression(one_master, coefficient=0)
assert finite_only == finite

# Arguments may contain composite expressions and ordinary Symbolica symbols:
m3, mu_r = S("m3", "mu_r")
other_master = olo.B0(psq, m2 + N(4), m3, mu_r)
other_coefficients = olo.get_expression(other_master)
```

`get_expression` accepts a master **Expression**, not a family string or a list
of strings. It accepts composite kinematics and relies on Symbolica's own
real/positive attributes and inference; there are no separate assumption lists.
Inspection accepts `oneloop::` as shorthand for `oneloopmaster::`; the numerical
coefficient hooks remain registered in `oneloopmaster::`.

### Select an analytic branch without substituting the kinematics

```python
all_branches = (finite, pole, double_pole)
selected = olo.select_branch(all_branches, [
    Replacement(m2, N(1)),
    Replacement(psq, N("3.23")),
])
print(selected[0])  # Parametric in psq and m2, with no if left for this point.
assert selected[1:] == (N(1), N(0))

# Exact integer and fractional probes are also supported.
selected_integer = olo.select_branch(all_branches, [
    Replacement(m2, N(1)), Replacement(psq, N(3)),
])
selected_fraction = olo.select_branch(finite, [
    Replacement(m2, N(1)), Replacement(psq, E("1/2")),
])
```

For the sample `psq=3.23, m2=1`, the actual selected finite expression has no
`if`, `abs` or `conj`. Naming its repeated roots, abbreviating the variables,
and collecting the two constant terms gives:

```text
r_plus  = (s + sqrt(s^2 - 4*s*M))/(2*s)
r_minus = (s - sqrt(s^2 - 4*s*M))/(2*s)
finite = 2 - log(M)
         - (1-r_plus)*log(1-1/r_plus)
         - (1-r_minus)*log(1-1/r_minus)
simple_pole = 1
double_pole = 0
```

This display applies to the selected below-threshold region `M>0, 0<s<4*M`.
The [verbatim before/after output](examples/expressions/b0_equal_mass.txt) and
[runnable Rust example](examples/select_branch.rs) preserve the actual expression
without introducing root aliases; its selected finite string is 758 characters
with the example's fully qualified symbol names.

Rules are applied only to temporary copies of the first arguments of native
`if(condition, yes, no)` expressions. Resolved conditions select a branch;
its kinematics stay symbolic. Every nested conditional in a retained branch or
condition is visited, without a helper-imposed depth cutoff. Undecidable
conditions remain parametric. Use the same namespaced symbols in replacements:
`S("m2")` does not match `S("mz_masses::m2")`.

This selects a formula valid in the probe's analytic region, **not** a globally
equivalent replacement of the original expression. Boundary probes select the
boundary formula. Do not reuse it across a branch cut without selecting again.
Native `if` tests zero/nonzero (including complex nonzero), not positivity.
Probe points must respect the input symbols' assumptions and the integral's
physical domain. Finite floating-point probes preserve their exact represented
values during condition substitution; temporary native guards defer arithmetic
to numerical evaluation and do not appear in the returned expression.

For nonzero positive `M=m²`, a **compact mathematical equivalent**, not a promise
of the current raw formatter's algebraic form, is:

```text
beta = sqrt(1 - 4*M/(s+i0))
finite = 2 - log(M) - beta*log((beta+1)/(beta-1))
simple_pole = 1
double_pole = 0
```

Here `i0` denotes the analytic continuation, not a finite numerical regulator;
principal roots/logs are understood with that prescription. At `s=0`, the finite
limit is `-log(M)`; at `s=4*M`, it is `2-log(M)`. Above threshold it has imaginary
part `+pi*sqrt(1-4*M/s)`. The emitted expression implements the continuation with
native conditionals, not a symbolic `i0` variable or opaque OneLOop functions.
Before branch selection it retains the massless/scaleless cases and redundant
guards: positive attributes are not used as proof of nonzero values. Omitting
the mass attribute allows a complex squared mass, subject to `Im(m²)<=0`.
See the [inspection guide](EXPRESSION_INSPECTION.md) for Rust branch selection,
assumptions, output examples and expansion limits. The standalone legacy numeric
extension does not expose symbolic inspection: use one shared Symbolica kernel.

The massive C0 case can be generated completely **before** selecting a branch:

```python
mass2, mass2B, mu2 = S("mass2", "mass2B", "mu2")
C0 = S("oneloopmaster::C0")
master = C0(0, -mass2, mass2, mass2, mass2, mass2B, mu2)
probes = [Replacement(mass2, N(2)), Replacement(mass2B, N(1)),
          Replacement(mu2, N(1))]
all_branches = olo.get_expression(master)  # ordinary Symbolica Expressions, with ifs
selected = olo.select_branch(all_branches, probes)
print(selected[0])  # Full native Symbolica expression, masses still parametric.
print(selected[0].evaluate({mass2: 2, mass2B: 1, mu2: 1},
                           decimal_digit_precision=60))
```

This family uses an exact four-dilogarithm representation, including zero-mass
branches and the prescribed real-axis lips. Run [the C0 example](examples/c0_selected_branch.py)
to verify the actual repository output against the known value and native backend.

For **fully generic C0/D0**, avoid duplicating every repeated root and branch
subexpression: request complete shared storage, then select afterwards:

```python
p1, p2, p3, m1, m2, m3, mu2 = S("p1", "p2", "p3", "m1", "m2", "m3", "mu2")
master = C0(p1, p2, p3, m1, m2, m3, mu2)
complete = olo.get_expression(master, shared=True)
print(complete[0])  # root, binding count, storage size; not a huge formula dump
# complete[0].root and every (alias, body) in .definitions are Symbolica Expressions.
probes = [Replacement(x, N(v)) for x, v in zip(
    [p1, p2, p3, m1, m2, m3, mu2], [-1, -2, -3, 4, 5, 6, 7])]
selected = olo.select_branch(complete, probes)  # plain, parametric Expressions
```

`shared=True` stores **all branches**, with every defining expression included.
It uses Symbolica's native `AliasedAtom`, not master evaluation hooks or opaque
OneLOop helpers. Rust exposes `get_expression_shared(master)` and
`select_branch_shared(coefficient, &rules, options)`; each shared coefficient
supports native `.evaluator(&parameters)`. Python exposes `SharedExpression.root`,
`.definitions`, `.num_definitions`, `.byte_size`, and bounded `.to_expression()`.
Do not evaluate the root alone without its bindings.

Release-mode construction of the generic shared D0 takes approximately 10–13 seconds
on the development machine; subsequent complete branch selection takes about one
second. The ordinary C0 case above takes milliseconds. These are construction
timings, not numerical-backend benchmarks. Unrestricted **fully duplicated**
generic C0/D0 trees can still exceed the resource budget; the shared form is the
scalable complete-expression API. See [the inspection guide](EXPRESSION_INSPECTION.md)
and [measurements and validation scope](performance/2026-09-25-expressions/README.md).

Early selection remains available as `get_expression(master, branch_rules=probes)`
or Rust's `get_expression_on_branch`. Probe values never replace retained bodies.
Partial probes may retain conditionals; flattening a large unresolved result can
still exceed the budget.

## Compact master Symbols

`oneloop::{A0(), B0(), dB0(), C0(), D0()}` return the Symbolica
Symbols `oneloopmaster::{A0,B0,dB0,C0,D0}` with dedicated `EvaluationInfo` hooks.
The crate registers a lightweight Symbolica `initialize!` hook, ordered after
Symbolica's special-function initializer. It attaches master and JIT callback
symbols only; parsing a master needs no prior accessor call. Numerical backends
are prepared per family on first use. Transparent definitions are built only
when an expression API requests them. Importing Python does not build formulas,
load evaluator assets, or compile SymJIT code.

Call `oneloop::initialize()?` explicitly to warm up all definitions and all five
Native and SymJIT backends. `is_initialized()` reports completion of this optional
full warmup; it stays false during ordinary lazy use.

Each master call takes a leading Laurent-power tag `0`, `-1` or `-2`, then
the same numeric arguments as the lowercase constructor, flattened into one list.

```rust,ignore
use symbolica::prelude::*;

let p = parse!("p");
let finite = oneloop::B0().call((0, &p, 1, 1, 1));
let pole = oneloop::B0().call((-1, &p, 1, 1, 1));
let context = oneloop::OneLoopExpressions::new();
let exact_evaluator = context.evaluator(&[finite, pole], &[p])?;
```

The context supplies transparent tagged native definitions; Symbolica gives them
precedence over the hooks. Bare calls can instead evaluate through cached hooks
for `Complex<f64>`, `Complex<DoubleFloat>` and `Complex<Float>`. These hooks
use the expression interpreter by default, or generated Rust when
`generated-evaluators` is enabled, independently of the manual backend selection.
Both routes cache their evaluators rather than constructing them per point.
They are opaque to the *outer* evaluator's optimizer, so prefer the context
when optimization across combined expressions is wanted. They do not supply
every Symbolica numeric domain or automatic error-controlled precision. Invalid tags/arities produce descriptive
errors or panics; numeric callbacks have no `Result` return channel. Stack and
license requirements below apply to either route. See
[MASTER_SYMBOL_AUDIT.md](MASTER_SYMBOL_AUDIT.md) for verification status.

To compile the bare Symbols while retaining their cached evaluation hooks, omit
the function map. This example keeps the squared scale as a runtime parameter:

```rust,ignore
use symbolica::prelude::*;

let parameters = [parse!("p"), parse!("m0"), parse!("m1"), parse!("mu_squared")];
let calls = [0, -1, -2].map(|tag| {
    oneloop::B0().call((tag, &parameters[0], &parameters[1],
                      &parameters[2], &parameters[3]))
});
let exact = Atom::evaluator_multiple(&calls, &parameters)
    .direct_translation(true).build()?;
let mut evaluator = exact.jit_compile::<Complex<f64>>(oneloop::jit_settings())?;
let input: Vec<_> = (0..1024).flat_map(|_| [-1., 1., 1., 4.])
    .map(|value| Complex::new(value, 0.)).collect();
let mut output = vec![Complex::new(0., 0.); 1024 * 3];
evaluator.batch_evaluate(&input, &mut output, 1024);
```

The outer evaluator is compiled here; each native family backend is prepared
on its first numerical call. For transparent expressions, use
`OneLoopExpressions::jit_evaluator` or `MappedLaurentSeries::jit_evaluator`.
These methods lower complex square roots and branch conditions to registered
numeric callbacks before JIT compilation. Calling Symbolica's `jit_compile`
directly on an exact evaluator bypasses those compatibility fixes. The exact
interpreter uses upstream Numerica logarithms and magnitudes directly; its
mixed-precision fixes preserve component precision. Inspection retains the
original formulas, and native code generation emits the existing logarithm
primitive.

## Numerical evaluators and batches

Each prepared family evaluator returns all three Laurent coefficients. Its
parameters are runtime inputs, including
`mu_squared` as the **last argument of every row**; the scale is not baked in.

| API | Numerical backend | Numeric types |
| --- | --- | --- |
| `NativeEvaluator<T>` | Generated, ahead-of-time Rust arithmetic | `f64`, `DoubleFloat`, `Float` |
| `ScalarEvaluator` | Portable SymJIT O2, with scalar binary64 batches | `f64` |
| `PrecisionEvaluator` | Explicit Native or exact-expression interpreter | `Float` |

`NativeEvaluator` requires `generated-evaluators` and owns reusable constants
and a workspace. Its arithmetic runs
in generated Rust. Binary64 dilogarithms call Symbolica's
`transcendental::dilog_complex_f64` kernel directly (no callback lookup, no
MPFR); `DoubleFloat` and `Float` use Symbolica's registered polylog callback, so
first use can initialize Symbolica and the requested OneLOop backend. Rational
constants are converted directly at construction.
Clones have independent workspaces; batches run the generated scalar function
for each row, without a SIMD promise or hidden worker threads. See
[the native evaluator](src/native/mod.rs) and
[its generic primitives](src/native/primitives.rs).

```rust,ignore
use oneloop::{NativeEvaluator, ScalarIntegral};
use symbolica::domains::float::{Complex, DoubleFloat};

let mut evaluator = NativeEvaluator::<DoubleFloat>::new(ScalarIntegral::B0)?;
let input = [-1., 1., 1., 4.].map(|x| {
    Complex::new(DoubleFloat::from(x), DoubleFloat::from(0.))
});
let zero = Complex::new(DoubleFloat::from(0.), DoubleFloat::from(0.));
let mut output = [zero; 3];
evaluator.evaluate(&input, &mut output)?;
```

### Tensor coefficients

`NativeTensorEvaluator<T>` evaluates OneLOop's tensor coefficients together, in
the same generic domains and with the same batch contract as `NativeEvaluator`
(row-major input with `mu²` **last and squared**, every row validated before any
output is written, shape/domain errors leave the output unchanged):

| `TensorIntegral` | Input row | Output row (Laurent triples, in order) | OneLOop routine |
| --- | --- | --- | --- |
| `An` | `[m², mu²]` | `A0, A00, A0000` | `olo_an`, rank 4 |
| `B11` | `[p², m0², m1², mu²]` | `B0, B1, B11` | `olo_b11` (without `B00`) |
| `B111` | `[p², m0², m1², mu²]` | `B0, B1, B11, B111` | `olo_bn`, rank 3 |

Bubble denominators are `(q²-m0²)((q+p)²-m1²)`, with `q^mu -> p^mu B1` and
`q^mu q^nu -> g^{mu nu} B00 + p^mu p^nu B11`, in OneLOop's normalization; each
coefficient is `[finite, 1/epsilon, 1/epsilon²]`. `B11` is what a rank-`n`
reduction such as Ninja's requests; `B111` serves rank `n+1`, where Ninja's own
Passarino-Veltman reconstruction divides by `p⁶`.

The arithmetic is a hand-written port of OneLOop's `bub0`-`bub111`, `bnlog`,
`solabc` and `tadpn` (see [the module](src/native/tensor.rs)): exact zero
dispatch and explicit real-axis lips replace OneLOop's `EPSN` tolerances and
`IEPS` regulator, and `log(1-1/x)` is shared by all ranks of a root, with
OneLOop's logarithmic series of `bnlog` for `|log(1-1/x)| < 1/2`. It needs no
generated code, so it is available without `generated-evaluators`. For `f64`
the same algorithm runs as a binary64 specialization (real arithmetic for real
kinematics, Estrin's scheme for the series), at most 0.6x of OneLOop's
`olo_b11` time on every bubble call of Ninja's benchmark streams. The
[tensor tests](tests/native_tensor.rs) compare all coefficients with original
OneLOop (366 rows from [the tensor oracle](tests/support/tensor_oracle.f90)), with
320-bit evaluation, with the exact symbolic `an`/`b0`/`b1`/`b11`/`b111`
expressions at 256 bits, and with Ninja's rank-three formulas.

```rust,ignore
use oneloop::{NativeTensorEvaluator, TensorIntegral};
use symbolica::domains::float::Complex;

let mut bubbles = NativeTensorEvaluator::<f64>::new(TensorIntegral::B11)?;
let input = [15625., 29929., 29929., 10000.].map(|x| Complex::new(x, 0.));
let mut output = [Complex::new(0., 0.); 9]; // B0, B1, B11
bubbles.evaluate(&input, &mut output)?;
```

For `NativeEvaluator`, fixed precisions are 53 bits (`f64`) and 106 bits
(`DoubleFloat`); `Float` defaults to 128 bits or accepts explicit working bits
through `with_binary_precision`. Basic DoubleFloat arithmetic and elementary
functions retain both components. Its Li2 primitive alone uses Symbolica Float
at 160 bits and converts back without narrowing the input to f64.

The API also preserves manual expression construction and provides
`OneLoopExpressions::jit_evaluator` / `MappedLaurentSeries::jit_evaluator` for
compiling combinations to SymJIT O2.

```rust,ignore
use oneloop::{ScalarEvaluator, ScalarIntegral};
use symbolica::domains::float::Complex;

oneloop::initialize()?; // All symbols and enabled numeric backends ready.
// Clone the prepared O2 backend; no additional IR loading or JIT compilation.
let mut evaluator = ScalarEvaluator::prebuilt(ScalarIntegral::B0)?;
let rows = [-1., 1., 1., 1.,  // p², m0², m1², mu²
            -1., 1., 1., 4.]; // same kinematics, another squared scale
let input: Vec<_> = rows.into_iter().map(|x| Complex::new(x, 0.)).collect();
let mut output = vec![Complex::new(0., 0.); 2 * 3];
evaluator.evaluate_batch(&input, &mut output, 2)?;
// output[0..3] belongs to the first row; output[3..6] to the second.
```

`evaluate` handles one point; `evaluate_batch` accepts flat row-major arrays and
checks their exact dimensions, including empty batches. SymJIT currently runs
scalar code for each row; its upstream SIMD path fails zero-input batches.
For default manual evaluation without constructing an evaluator object,
`oneloop::evaluate(family, input, output)` and
`oneloop::evaluate_batch(family, input, output, rows)` use `DEFAULT_BACKEND`,
which is `EvaluationBackend::Expression` by default and `EvaluationBackend::Native`
with `generated-evaluators` enabled. Select a route explicitly with
`evaluate_with_backend(family, input, output, EvaluationBackend::Native)` or
`evaluate_batch_with_backend(family, input, output, rows, backend)`; `Expression`
selects Symbolica's native expression interpreter. These convenience routes
prepare only the requested backend and family.
The historical measurements below used an older patched SIMD backend and do
not establish performance for this revision. The selector does not control
the registered Symbol hooks. `ScalarEvaluator::cached(family)` clones the
prepared SymJIT backend into an independently reusable workspace; an explicitly
compiled custom expression keeps its own evaluator. See the
[backend selectors](src/backend.rs).
`to_bytes` / `from_bytes` serialize numerical evaluator instructions **and nested
function definitions**, not machine code or process pointers. Loading recompiles
every level with the same strict JIT settings. Cache compatibility is tied to
the Symbolica IR schema and compilation policy; no backend machine code is
reused across SymJIT versions. Load only trusted artifacts. Float/arbitrary-precision evaluators use
`PrecisionEvaluator` or the original exact-expression API, not these f64 blobs.
Those evaluators are constructed only for the requested precision.
Keep `SYMJIT_TOML` unset and do not place `symjit.toml` in the working directory
when building/loading portable caches: these overrides can pin SymJIT to a
particular machine architecture, so the portable API rejects them. Register
any custom numerical Symbol callbacks before loading a manually composed cache.

### Build size and optional generated evaluators

The default build enables `prebuilt` and uses the expression interpreter.
The opt-in `generated-evaluators` feature compiles the ahead-of-time Rust
formulas for `f64`, `DoubleFloat`, and `Float`. These generated formulas,
especially the box integral, account for most of OneLOop's machine code.
To enable them for faster evaluation at the cost of a larger library:

```toml
oneloop = { git = "https://github.com/alphal00p/oneloopmaster", branch = "main", features = ["generated-evaluators"] }
```

Without `generated-evaluators`, automatic evaluation and the registered master
callbacks use the expression interpreter. Arbitrary precision stays arbitrary
precision, and `DoubleFloat` stays `DoubleFloat`. Explicit `NativeEvaluator`
construction or `EvaluationBackend::Native` returns a configuration error.
Rust callers can query `EvaluationBackend::Native.is_available()` first.
SymJIT remains explicitly selectable. This trades throughput and first-use
formula preparation for a smaller library; it does not change the integral API.

The Python adapter and development host also omit generated evaluators by
default and forward both features. To opt into generated evaluators, use:

```sh
cd python
maturin build --release --features generated-evaluators
```

Community hosts must ensure no dependency enables `generated-evaluators` to
keep the smaller build; Cargo combines features enabled by every caller.
Enable `community,prebuilt` on the adapter for expression interoperability and
cached SymJIT support. A host can expose its own opt-in feature forwarding
`oneloop/generated-evaluators` and `oneloop-python/generated-evaluators`.
Disabling `prebuilt` separately also removes the embedded SymJIT caches, at the
cost of generating them on first use.

To compare warm numerical throughput in all three numeric domains, run
`cargo run --release --features generated-evaluators --example numeric_performance -- 100 5 native` and repeat
with `expression` as the final argument. The survey validates the independent
Fortran fixtures before timing, and excludes construction and input conversion.

The [2026-10-02 Linux measurements](performance/2026-10-02-size/results.json)
compare against commit `06a72bf`: sharing `DoubleFloat` arithmetic in the large
generated graphs reduces generated machine code from 27.10 to 12.38 MiB.
The standalone adapter wheel shrinks from 17.19 to 13.47 MiB with generated
evaluators enabled, or 9.97 MiB with the current defaults (generated evaluators
disabled and `prebuilt` retained).
These are local unstripped x86-64 release builds, not total community-wheel sizes.
Small-family arithmetic stays inline to avoid adding call overhead to A0/B0/dB0.

The `prebuilt` feature is enabled by default and embeds all five evaluators using
`include_bytes!`. Asset provenance, sizes and hashes are in
[assets/evaluators/README.md](assets/evaluators/README.md). To rebuild explicitly,
call `ScalarEvaluator::rebuild` or run:

```sh
cargo run --release --no-default-features --example rebuild_evaluators -- assets/evaluators
```

Without `prebuilt`, first use builds the requested SymJIT family from native
expressions. An optional final family argument writes only that rebuilt asset.
Applications can explicitly replace a shared SymJIT manual cache with
`rebuild_cached_evaluator(ScalarIntegral::B0)`; this does not replace the direct
native Symbol-hook cache. Native Rust source is generated by
`cargo run --release --features rust-codegen --example generate_native` and compiled with the crate,
not loaded from the portable blobs. Cache regeneration and fresh-load batch
checks must pass for each dependency update. A batch
method does not guarantee a speedup for divergent conditional branches.

The opt-in `rust-codegen` feature also exposes the same generator as
`oneloop::rust_codegen::RustEmitter`. It accepts an exact Symbolica evaluator,
retains shared functions and lazy branches, and returns a Rust module with any
entry output count. `constants_source()` emits its shared exact constant table;
`constants()` exposes the descriptors for inspection. Failed exports leave
previously emitted constant indices unchanged. The `generate_native` example
uses this API for all five scalar families.

Generated child modules use the parent imports
`oneloop::rust_codegen::runtime::{C, Context, NativeFloat, ConstantSpec, primitives}`.
Construct `Context::new(generated_constants::CONSTANTS, bits)` once, then call
each child's `evaluate`. Constants are parsed directly at working precision;
inputs are neither converted nor padded. The caller validates input counts and
precision. `emit` emits scalar arithmetic; wrapping that entry in a row loop remains a
scalar control. Shared helper outputs must be scalar, and unsupported external
functions return an export error.

`RustEmitter::emit_batch` emits the scalar entry above and an additional
instruction-major `evaluate_batch(input, output, rows, ctx, workspace)` entry.
Re-export `BatchWorkspace` and `batch` from the same runtime in its parent module,
and construct `BatchWorkspace::new(prototype)` once. Public buffers are row-major;
retained scratch columns are slot-major. Each instruction runs across its active
rows, including inside shared helpers. Storage follows Symbolica's exported parameter,
temporary and output slot counts, rather than the number of emitted statements.
Empty batches do no arithmetic; single-row calls use the scalar entry without
scratch packing or allocation. Dimension/size errors leave output unchanged.

Conditional graphs retain Symbolica's lazy forward branches: inactive singular
arms do no arithmetic. Shared helpers receive the parent's stride and active
lane indices, and consecutive joins preserve the selected component values.
Active/resume scratch is retained per canonical helper. Entirely branch-free
graphs keep their existing contiguous loop path. Malformed or backward control
flow is rejected before export changes the constant table. Larger batches never
fall back to scalar row loops. The existing scalar master evaluator backend
still uses its established row loops; emitting lazy batch kernels does not
automatically replace that backend.

Scratch initialization uses the explicit prototype, while inputs and constants
retain their own precision (including different Float precision per component).
`Context` still owns exact cached constants. A workspace can be reused across
modules, contexts and batch sizes. Symbolica's `Dualizer` supplies flattened
higher/mixed Taylor components in its original ordering; no Dual components are
dropped and no new differentiation rules are implemented. The runtime uses
existing native primitives, so Float operations may allocate even when the
outer scratch storage is warm. No SIMD or speedup guarantee is implied.

Run `cargo test --features rust-codegen --test rust_batch_codegen -- --nocapture`
for compiled f64/DoubleFloat/Float and flattened-Dual comparisons, an explicit
instruction-major primitive trace, warmed f64 allocation checks, and generated
source-size/optimized-compilation measurements.

SymJIT has a compatible minimum requirement of `2.27.3`; the consuming
application's `Cargo.lock` selects the concrete version. The six previously
reported standalone defects passed their regressions in 2.26.0. The newer
minimum supports configurable call argument storage: portable evaluators size
it from all nested function calls, counting two real slots per complex
argument, and reject requests above the backend's accepted capacity. Source
cache loading recomputes this size before compiling on the receiving host.
Production retains ordinary compilation, packed complex arithmetic, `fastmath=false` and scalar
batches while the broader backend matrix is audited. Complex square-root
callbacks remain necessary for extreme scales. See the
[dependency audit](patches/upstream-status-2026-09-21.md) for reproduced remaining
limitations and the standalone cases.


## Python (PyO3)

The separate [Python adapter](python/README.md) supplies scalar and batched
complex-number evaluation through PyO3. It adds no Python/PyO3 dependency to the
core Rust crate. The recommended [development host](python/host/) builds one
`symbolica.core` containing both Symbolica and OneLOop, so all symbolic APIs use
actual Expression and Replacement objects. It is a local host, not an upstream
community-repository registration. The [HEP integration assessment](HEP_INTEGRATION.md)
records Feynkit type reuse and the host registration interface. Install it in an isolated environment without
a competing Symbolica distribution:

```sh
cd python/host
maturin develop --release
```

For example, on an adequately stacked calling thread:

```python
from symbolica import S
from symbolica.community.hepkit import oneloop as olo

bubble = olo.Evaluator(olo.B0)  # Reuses a direct Rust numeric workspace.
values = bubble.evaluate_batch([
    [-1.0, 1.0, 1.0, 1.0],
    [3.0, 0.7 - 0.03j, 1.4 - 0.08j, 4.0],
])
# Each row returns (finite, simple_pole, double_pole); final input is mu_squared.
```

In the community module, `olo.A0`, `B0`, `dB0`, `C0`, and `D0` are callable
Symbolica Expressions exported from this crate's primitive master accessors.
An untagged call such as `olo.A0(m2, mu2)` is suitable for `get_expression`;
`olo.A0(0, m2, mu2).evaluate({m2: 2, mu2: 1})` evaluates its finite
coefficient through the existing native Rust hook. Tags `-1` and `-2` select
the simple and double poles. Supply the squared scale explicitly as the last
argument when constructing a symbolic call. These exports resolve and cache
on first attribute access without building formulas or numerical backends.

Complete tagged calls automatically normalize to numerical Expressions when
all kinematic arguments are numbers and at least one is floating point. For
example, `olo.A0(0, 2.0, 1)` evaluates the finite coefficient during construction.
Mixed exact and floating inputs use the maximum supplied floating-point
precision, including complex components and the scale, without a binary64
conversion of exact inputs. All-exact calls such as `olo.A0(0, 2, 1)` retain
their symbolic meaning; request `.evaluate({})` to choose numerical evaluation.
Untagged calls, partially symbolic calls, and invalid tags, arities, or numerical
domains remain symbolic. Nonfinite numerical results also leave the call intact.


The optional legacy build, `cd python && maturin develop --release`, exposes
numeric-only `oneloop_native`; uppercase names remain numerical aliases there.
Its family-name selectors remain a compatibility
interface, but it no longer accepts or returns symbolic strings. See the
[API audit](API_BRANCH_AUDIT.md) for this change's checks and remaining limits.

Numerical functions `olo.a0/b0/db0/c0/d0(..., mu_squared=1, prec=16)` also
accept `rebuild=True` and `backend="auto"`. Reusable evaluators and their methods
likewise accept a backend choice: `"native"`, `"symjit"` or `"expression"`.
`prec` counts **decimal digits**. With `"auto"`, ordinary Python numbers at the
default precision use binary64 arithmetic and return Python `complex` results.
Higher precision, large integers, `Decimal` inputs, or `DecimalComplex` inputs
use Float and return `DecimalComplex` coefficients. Automatic calls use the
expression interpreter unless `generated-evaluators` enables direct Rust.
`backend="expression"` selects the transparent expression interpreter;
`"symjit"` selects the embedded O2 backend and rejects arbitrary-precision requests.
Use `backend="symjit", rebuild=True` to recompile that backend. Rebuilding a
Native evaluator only prepares fresh constants and a numeric workspace.

```python
from decimal import Decimal
import oneloop_native as olo

mass = olo.DecimalComplex(Decimal("2.0000000000000000000000000000001"),
                          Decimal("-0.0000000000000000000000000000003"))
finite, pole, double_pole = olo.a0(
    mass, mu_squared=Decimal("4"), prec=1000, backend="native")
print(finite.real, finite.imag)  # Decimal components, not binary64 approximations.

bubble = olo.Evaluator("B0", prec=32, backend="native")
values = bubble.evaluate_batch([
    [Decimal("-1"), Decimal("1"), Decimal("1"), Decimal("1")],
    [Decimal("-2"), Decimal("1"), Decimal("1"), Decimal("4")],
])
```

Supply decimal text to `Decimal`, not an already rounded Python float. Decimal
construction and conversion preserve the supplied digits independently of the
global decimal context; subsequent arithmetic on returned Decimal components
uses Python's context as usual. See [the Python precision guide](python/README.md#arbitrary-precision)
for return types, per-call overrides, and complex inputs. Arbitrary precision
uses `Complex<Float>` through the selected Native or Expression backend, not
SymJIT; the binary64
performance measurements below do not bound high-precision runtime.
The optional [Python performance survey](python/performance_survey.py) defaults
to 1024-row batches and records Python conversion/allocation overhead separately
from initialization. The completed 1024-row survey includes the one-row tail,
all three outputs, and matched Rust/Fortran measurements; see the
[performance report](performance/2026-09-08-native/README.md).

## Stack and precision

Explicit full warmup and generic C0/D0 compilation traverse a large conditional expression graph.
The Rust regression tests run the entire Symbolica workload on a thread with a 128 MiB
stack. For a standalone program, create that thread **before constructing any
Symbolica values**, and run construction/evaluation inside it. A restricted
Symbolica license permits only one active Symbolica thread per user; do not
overlap Symbolica work across threads or processes under that license.
Applications already using Symbolica must arrange an adequate calling-thread
stack rather than silently spawning another Symbolica worker. Both Python cold
startup checks pass on main threads capped at 8 MiB, but this is not a general
stack bound for arbitrary expressions or explicit rebuilding. Sequential thread
handoff is separately tested; the initializing thread need not remain alive for
the verified cached-evaluator use case.
Raw `NativeEvaluator<T>` uses the calling thread and its owned numeric workspace.
Its first dilogarithm resolves Symbolica's public polylog callback and can trigger
global initialization. The same calling-thread stack and license requirements
therefore apply to that first use, as well as Symbol access and Python import.

`PrecisionEvaluator` evaluates the same formulas at a fixed arbitrary working
precision, with explicit Native or Expression selection and no additional direct
dependency. `new` defaults to Expression, or Native when `generated-evaluators`
is enabled; SymJIT is binary64-only.
An explicit direct native request looks like:

```rust,ignore
use oneloop::{EvaluationBackend, PrecisionEvaluator, ScalarIntegral};
use symbolica::domains::float::{Complex, Float};

let mut evaluator = PrecisionEvaluator::with_backend(
    ScalarIntegral::B0, 1000, EvaluationBackend::Native)?;
let bits = evaluator.binary_precision();
let input: Vec<_> = ["-1", "1", "1", "4"].into_iter().map(|text| {
    Complex::new(Float::parse(text, Some(bits)).unwrap(), Float::new(bits))
}).collect();
let mut output = core::array::from_fn::<_, 3, _>(|_| {
    Complex::new(Float::new(bits), Float::new(bits))
});
evaluator.evaluate(&input, &mut output)?;
println!("{} + ({}) i", output[0].re, output[0].im);
```

The complete [arbitrary-precision example](examples/arbitrary_precision.rs)
arranges the calling-thread stack and evaluates two rows; run it with
`cargo run --release --example arbitrary_precision -- 1000`.

`new` / `with_backend` add fixed guard bits to the requested decimal precision;
`with_binary_precision` / `with_binary_precision_and_backend` instead use exactly
the supplied bit count. Raw `NativeEvaluator::<Float>::with_binary_precision`
also uses exactly its supplied bits, without eager symbolic startup.
The API supports reusable scalar and flat row-major batched evaluation
with the same argument order as `ScalarEvaluator`. Exact rational coefficients
are converted directly at working precision. Public master-Symbol hooks support
`Complex<Float>` with native workspaces keyed by the maximum input-component
precision (including the squared scale), and `Complex<DoubleFloat>` at 106 bits.
Map a bare-Symbol evaluator with `map_coeff_with_prec` to use Float hooks, or map
the transparent expressions directly.

Use exact rational or decimal-text inputs when testing precision convergence.
Mapping coefficients and supplying values as `Complex<f64>` cannot recover digits
already rounded out of the input. Increasing precision provides more working
digits, not a guarantee of that many accurate output digits at every kinematic
point. Outputs retain Symbolica/Numerica's computed precision; we do not pad lost
digits back to the requested precision. No custom adaptive precision or
near-degeneracy gate is implemented here. See [PRECISION.md](PRECISION.md) for
the precision contract and validation coverage.

## Verification

Run local checks from this directory:

```sh
cargo test -- --test-threads=1
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

The original acceptance table has 293 points / 879 Laurent coefficients. The
broader retained audit has 342 accepted Fortran points / 1,026 coefficients;
run its fixed-f64 comparison explicitly with:

```sh
ONELOOP_AUDIT_FIXTURES=tests/data/scalar_audit.txt cargo test --test parity \
  audit_extra_fortran_fixtures -- --ignored --test-threads=1 --nocapture
```

The earlier generated-evaluator release suite passed **114 tests**, with six ignored
entries (four opt-in diagnostics and two subprocess helpers invoked by parent
tests). The direct native backend passes the broader table at 256 bits; its
explicit binary64 diagnostic still fails on 14 small-momentum B0/dB0 coefficients:

```sh
cargo test --release --test native native_expanded_binary64_diagnostic -- \
  --ignored --test-threads=1 --nocapture
```

This is a retained failing test, not an accepted tolerance exception. No
comparison tolerance has been relaxed
and no disagreement is converted into a passing assertion. Passing the sampled
tables or default suite is not a full-parity certificate. The acceptance and
expanded tables have 613 unique ordered inputs, not 635 independent points.
See [the exceptional-vacuum audit](VACUUM_SYMMETRY_AUDIT.md) for the exact disputed
permutations, independent analytic targets, and gaps in region coverage. Analytic
symmetry takes precedence over inconsistent reference outputs there. See the dated audit
for the latest completed runs and additional exact-boundary findings.
The restored strict-v3 scalar caches separately pass all three outputs for the
293 acceptance rows and 31 benchmark rows. The final release Python suite passes
all **27 tests**, including native/SymJIT selection, arbitrary precision, Decimal
inputs, inspection and 1024-row batches. Arbitrary-precision, callback, mapped-expression,
portable-cache and batched routes are tested separately, not inferred from one
another. Bare master calls compiled without a FunctionMap pass every acceptance
row in canonical tag order, scalar and five mixed/repeated batch layouts, before
and after portable restoration. Five additional tag orders cover controls from
each family. The exact repeated-root box regression passes at 128/256 bits; the
nearby width-1e-30 diagnostic at 128 bits asserts finiteness only, with its target
comparison performed at 256 bits.

The performance survey compares the same 31 numeric fixtures, all three Laurent
outputs, and 256 warmup calls against the unchanged original Fortran source.
It separates explicit full warmup from warmed evaluation, checks uniform and mixed rows,
and reports batches 1/4/32/256/1024 individually:

```sh
python3 tests/benchmark_fortran.py /path/to/OneLOop --iterations 2000 \
  --warmup 256 --repetitions 7 --output /tmp/oneloop-fortran.json
cargo run --release --example performance_survey -- 2000 7 all jit \
  > /tmp/oneloop-native.tsv
python3 tests/compare_performance.py /tmp/oneloop-fortran.json /tmp/oneloop-native.tsv
```

Use `all native` instead of `all jit` to measure the generated Rust evaluator;
its batch loop is scalar. The [current report](performance/2026-09-08-native/README.md)
measures this backend and native-only hooks separately from manual SymJIT.
Use `all symbols` to compile public master calls with their
transparent native FunctionMap, using the same batch sizes. Use `all hooks-jit`
to JIT-compile bare master calls without that map, invoking the automatically
registered direct native callbacks. `all hooks` keeps that bare outer evaluator
interpreted and is scalar-only (compare with `--batch 1`). These are distinct
measurements; mapped compilation is not the automatic callback path.
Binary64 hooks reuse one backend result for distinct Laurent tags at
bitwise-identical inputs.
A bounded four-entry FIFO accommodates tag-major SIMD lanes, including identical
points, but each tag is consumed only once per group: repeated same-tag calls
cannot be memoized away. Seven counting-backend tests cover every tag order,
successive 1/2/4-lane groups, signed zeros, changed inputs, eviction, explicit
invalidation and failures during evaluation. More heavily interleaved points or
threads can reduce reuse. The outer JIT callbacks still enter the scalar family
backend per point; they do not combine those callbacks into a family SIMD batch.

The optional fifth survey argument selects batches, for example
`4097 7 all jit 1,1024`. This evaluates four full 1024-row batches plus a one-row
tail per repetition. Compare each selected batch separately with `--batch`.
The Python survey in `python/performance_survey.py` defaults to batches of 1024;
its timings include Python input conversion and output-object allocation, unlike
the Rust-only route. Compare runs on the
same machine without competing benchmark/build jobs, with matching call counts
and repetitions. The comparator rejects incomplete or invalid measurements and
returns exit 2 if any selected median exceeds 1.5 times the Fortran time.
The 1.5× target is not met across all workloads; full measurements and limitations
are in the [performance report](performance/2026-09-08-native/README.md).
The comparator's independent protocol tests
run with `python3 -m unittest discover -s tests -p test_compare_performance.py -v`.

The permanent Fortran fixture generator is standard-library-only Python. Compile
the small test driver against an existing, unchanged original OneLOop build:

```sh
gfortran -I /path/to/OneLOop tests/support/oracle.f90 \
  /path/to/OneLOop/libavh_olo.a -o /tmp/oneloop-oracle
python3 tests/generate_fixtures.py /tmp/oneloop-oracle
```

It prints the fixture table to standard output and rejects rejected/incomplete or
nonfinite Fortran results. The original example program does not print the full
A0/B0/dB0 outputs needed by this table, so use this test driver. The original
library is not a production dependency. Set its compiler runtime library path if required
by your local toolchain. The checked-in fixture table makes normal tests independent
of Fortran, Python, and the standalone numerical port.

For a whole-table high-precision exploratory run, also set
`ONELOOP_AUDIT_BITS=256` (or another precision) in the audit command above.
This selects one numeric domain for every point, not an error-triggered fallback.
`branch_limits`, `triangle_limits` and `box_limits` separately retain focused
precision/branch regressions; see the dated audit for which checks pass.

See `AUDIT.md` and `STATUS.md` for verified coverage and known incorrect or
unsupported cases. This is a development implementation, not a claim of global
Fortran equivalence. Cargo publishing remains disabled. See
`COMMUNITY_INTEGRATION.md` for the inspected Spenso-style integration plan.

## Citations

In the shared Python host, `symbolica.get_citations()` returns the package credit
and the two papers below after OneLoopMaster is used, alongside Symbolica’s
own citation. Entries are Symbolica `Citation` objects. Use
`citation.to_bibtex()` to export an entry. Reading citations does not initialize
numerical backends.

## Credits and references

The mathematical and algorithmic basis of this transcription is
[OneLOop, by Andreas van Hameren](https://helac-phegas.web.cern.ch/OneLOop.html).
The original Fortran implementation is the independent reference used by the
scalar parity tests. This repository is a separate, still-incomplete Rust /
Symbolica transcription, not an official replacement or a claim of endorsement
by the original author. The original
[OneLOopBridge wrapper](https://github.com/SecretGmG/OneLOopBridge) remains unchanged.

Please credit the original work when using or building on these expressions:

- A. van Hameren, **OneLOop: for the evaluation of one-loop scalar functions**,
  *Computer Physics Communications* **182** (2011) 2427–2438,
  [arXiv:1007.4716](https://arxiv.org/abs/1007.4716),
  [doi:10.1016/j.cpc.2011.06.011](https://doi.org/10.1016/j.cpc.2011.06.011).
- A. van Hameren, C. G. Papadopoulos and R. Pittau,
  **Automated one-loop calculations: a proof of concept**,
  *JHEP* **09** (2009) 106,
  [arXiv:0903.4665](https://arxiv.org/abs/0903.4665),
  [doi:10.1088/1126-6708/2009/09/106](https://doi.org/10.1088/1126-6708/2009/09/106).

[Symbolica](https://symbolica.io/) supplies the symbolic expressions, function
maps and numerical evaluation infrastructure. It has its own licensing terms.
Permission for the intended distribution/relicensing of this transcription is
being handled separately; it has not been recorded as granted here. The manifest's
license identifier is provisional, not a grant of rights to the original code.
Cargo publication is disabled pending permission and technical readiness.
