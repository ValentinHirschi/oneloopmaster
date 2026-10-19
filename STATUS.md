# Native OneLOop status

This page preserves the September 8 development snapshot. See the
[Symbolica main migration report](SYMBOLICA_3_MIGRATION.md) for the current
dependency setup and validation.

## Contour box near a vanishing Gram determinant (2026-10-19)

Boxes with two nearly parallel external momenta (their invariant next to its
threshold, e.g. two Higgs bosons nearly at rest relative to each other in
gg → HHH) take the contour representation. The NinjaRS gg → HHH box of
`tests/d0_vanishing_gram.rs` lost 3.4e-11 to 3.6e-11 in binary64 (native,
SymJIT) and 7.5e-12 interpreted; it is now within 1.9e-14 to 4.1e-14 on every route.

- **Cause.** Two of the canonical relabelling's three `tfun` have a contour
  direction nearly parallel to their line (`g + beta h` 1.4e-5 of `g`), and
  one pole within 5e-6 to 7e-5 of the root of each line edge; each edge's
  `s3` term rounded `y - z` differently. The offset now comes from Vieta's
  product with the poles' own `a1 = g + beta h` where it cancels by more than
  64 ([PRECISION.md](PRECISION.md)); since the section below, only where `a1`
  itself cancels by more than 64.
- **NaN at an exactly vanishing Gram determinant in binary64.** `K` of one
  contour is `b (g + beta' h)` of the other; where it rounded to zero
  binary64 returned NaN, on every route: the reported box gave NaN in native
  code only, but at 09b2c74 28% to 40% of the 1-ulp perturbations of all six
  invariants of that box gave NaN on every binary64 route (hand-written,
  generated, interpreted and SymJIT). The second contour is used then: none
  of the 3,000 perturbations (1 to 6 ulps, 1e-16 to 1e-9 relative, equal
  `p1², p2², p4²`) is non-finite on any route, all within 1.7e-13.
- **Scan.** `tests/data/near_threshold_boxes.txt` (200 random physical boxes,
  `tests/generate_near_threshold_rows.py`, OneLOop 3.6 quadruple-precision
  references): every binary64 route within 5e-12 on 199 rows (3.2e-12 at
  most), where 34 were beyond it and 22 beyond 1e-9 before (OneLOop's double
  precision: 32 and 16). On 1,000 more such boxes (seeds 2 to 6) 4 to 6 per
  route remain beyond 1e-11 (at most 3.2e-11), against 110 to 117 before; on
  fresh seeds (101 to 105) up to about 8e-11.
- **Remaining.** Line 73 of the scan keeps 1.1e-10 to 3.2e-10, as before
  (OneLOop 5.5e-11): a nearly double contour direction in one `tfun` and
  nearly coincident poles next to an edge's end in the others, whose
  quadratic `s3` terms cancel.

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
  boxes (|Im m²| from 1e-3 to 0.5 of Re m²) 131 are wrong by 0.3 to 7 of the
  value, and the same 131 with the widths divided by 1e6: OneLOop's D0 does
  not tend to its real-mass value there.
- **The expressions take Denner and Dittmaier's sheets** for a mass with a
  width (arXiv:1005.2076, eqs. (3.82)-(3.83), stated there for any real
  momenta and complex masses): the unnormalized roots `u = X sqrt(m²)/|sqrt(m²)|`
  and the ratios `Y23/Y02`, `Y13/Y01`, `Y03/m²` each with its principal
  logarithm (`box_integral::complex_mass_sheets`); real masses keep
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
  41 wrong, the others up to 2.4e-7); with |Im m²| of 1e-11 of Re m²,
  2000 boxes within 3e-9 of the real-mass expressions. A direct
  Feynman-parameter integral (primary sectors, contour deformation, nested
  tanh-sinh) reproduces the reproducer at P = -1e-3 to 2e-11 and at
  P = -1e-6 to 1e-9. `tests/d0_one_mass_width.rs` checks 15 rows
  (`tests/data/one_mass_width.txt`, references from the formula at 40
  digits) in every relabelling. Real-mass one-mass boxes and every other
  mass pattern are bitwise unchanged on every route.
- Two-mass boxes with complex masses: see the next section.

## Two-mass boxes with complex masses (2026-10-19)

- **OneLOop's `boxf2` (and `boxf5`, which relabels opposite masses into it)
  is wrong for complex masses**, in the same way as `boxf1`: its roots `X`
  are normalized by the complex `sqrt(m3²)` and `sqrt(m4²)`, which turns them
  against Denner and Dittmaier's `u = X sqrt(m4²)/sqrt(m3²)` across the
  negative real axis. On that sheet the continued dilogarithm of an argument
  next to one (`X r24/r12` = 1.0005 in the reproducer) takes the wrong side
  of its branch point, and the terms are off by `4 pi²`. In
  `4 1.0 0.0462 -4.297 12.27 10.74 -31.21 -13.06 0 0 0 0 0.2161 -0.00173
  1.556 -0.00161` OneLOop 3.6 (both precisions) and these expressions gave
  0.0782 + 0.0831i for -7.16e-4 + 0.0831i (COLLIER quad, Denner and
  Dittmaier's formula in mpmath, and a direct Feynman-parameter integral to
  2e-8). Random sweeps, against COLLIER 1.2.9 in quadruple precision: per
  mass pattern (00MM, 0MM0, MM00, M00M, 0M0M, M0M0) 8 to 18 of 2000 boxes with
  |Im m²| from 1e-3 to 0.5 of Re m² wrong by 0.6 to 17 of the value (0.4% to
  0.9%; masses 0.1 to 20 and invariants of either sign from 0.03 to 30,
  log-uniform; the two-mass verifier's sampling, masses 0.03 to 10 and
  invariants uniform up to 30 or log-uniform from 1e-6 to 50, gave 63 to 83
  of 3000 per pattern, 2.1% to 2.8%, with OneLOop's double and quadruple
  precision wrong in the same rows: the rate depends on the sample); 2 to
  10 of 500 per pattern with one real mass; 3 to 11 of 1000 with the widths
  divided by 1e6; 4 to 10 of 500 next to a collinear limit and to
  `Y02 Y13 = Y01 Y23` (by up to 9e4, and 2 to 10 more lost 1e-8 to 0.1 in
  binary64); 2 of 200 boxes with a real mass on shell (all relabellings).
  OneLOop is wrong in the same rows.
- **The expressions take Denner and Dittmaier's two-mass sheets** whenever a
  mass has a width (eqs. (3.78)-(3.80), stated for general momenta and
  complex masses): `u`, the ratios `Y23/Y02` and `Y13/Y01` and the roots
  `r031 = r14/rotation`, `r032 = 1/(r14 rotation)` of `m4² + Y03 x + m3² x²`,
  each with its principal logarithm, through the same
  `box_integral::complex_mass_sheets` as the one-mass box (including its
  winding correction for roots on opposite lips). Those terms are one
  registered function of the roots, ratios and complements
  (`__olo_box_two_mass_sheets`, `box_integral::two_mass_sheets`), called
  from every two-mass sector: inlined, every occurrence of a root carried
  its whole expression, and the definitions grew from 97 MB to 209 MB (sector
  5 alone 98 MB), past the bounds of `tests/expression_construction.rs`;
  they are now 97.7 MB (+0.4%). With a real `m4²` the
  logarithms of `Y01/Y12` and `Y02/m4²` are differences of those of the
  `-i0` invariants. A real `m4²` on shell with a massless line (`Y01` or
  `Y02` zero) takes the formula with lines 0 and 3 exchanged, as Denner and
  Dittmaier prescribe (the same roots inverted); a real `m3²` on shell
  (`Y13` or `Y23` zero) drops that dilogarithm pair. Real masses keep
  OneLOop's sheets and their bits.
- **Validation** against COLLIER in quadruple precision (COLI and DD
  branches, fed the exact binary64 inputs; where they agree to 1e-10) and
  the formula at 40 digits: after the change no row of any of these sweeps
  is wrong. Binary64 worst cases on the native, generated, interpreted and
  SymJIT routes: 4e-13 for random and one-real-mass boxes, 1e-12 with
  narrow widths, 6e-9 to 1e-8 next to collinear limits (where the old
  expressions lost up to 6e-7 in rows they got right: the new ones are more
  accurate by more than 10 in 140 of those rows and less accurate by more
  than 10 in 5, all below 3e-11); at 256 bits all rows agree with the
  formula to the 18 printed digits. The 24 relabellings of all 255 formerly
  wrong rows agree within 4e-14 (collinear rows 1e-8).
  `tests/d0_two_mass_width.rs` checks 17 rows (`tests/data/two_mass_width.txt`)
  in every relabelling on every binary64 route, in DoubleFloat and at 128
  bits. One-mass, three-mass and four-mass boxes and all real-mass boxes are
  bitwise unchanged on every route (2400 boxes of all 16 patterns, 2000 real
  one-mass, 4800 real two-mass including collinear ones, 2400 real on-shell
  relabellings, the D0 calls of the recorded NinjaRS master streams, and
  9000 complex one-, three- and four-mass boxes), and the native code of the
  other families is byte-identical. D0 oracle sweeps (`tests/oracle_sweep.py`,
  seeds 1-4, native, generated and expression routes) give the same failure
  lists as before. Retired instructions per binary64 D0 call (perf stat,
  difference of 2001 and 1 repetitions) against 09b2c74: the 22 recorded
  NinjaRS master streams within 0.002%; the other figures depend on the
  sample. Real two-mass boxes +0.4% (23,673 to 23,763; 500 random boxes);
  complex two-mass boxes -2.5% to +4.4%: on this stream's 500 random boxes
  per sample -1.2% to -1.5% (00MM), -0.1% to -0.3% (0M0M) and -2.3% to -2.5%
  (00MM with one real mass), on the two-mass verifier's random samples +4.4%
  (00MM, 28,403 to 29,649), +2.7% (0M0M) and +0.3% (one real mass); the
  on-shell boxes with one real mass, in all relabellings, +14% to +24%
  (19,938 to 22,810 on this stream's 500 rows, 17,852 to 22,156 on the
  verifier's).
- **A reference corrected**: the second box of
  `tests/d0_collinear_limits.rs::collinear_limits_at_every_precision`
  (complex `m2²`, real `m3²`, next to a collinear limit) had OneLOop's
  quadruple-precision value, 0.0338 + 0.1177i, which is wrong: COLLIER's two
  branches and the formula give -0.1774 + 0.1177i (to 1e-19), and so do
  these expressions now. Its reference is the formula's value; the
  tolerance is unchanged.
- **NinjaRS**: D0 values of complex two-mass boxes change in their last bits
  (where they were right). Six cases of the rank <= n+1 bitwise fixture
  (`edge-low-rank-n5-r0` to `-r3`, `edge-many-n7-r6`, `-r7`: their `mixed`
  lanes) contain such boxes (44 in all; right before and after, within
  2e-14 of COLLIER), so their bits change once NinjaRS uses this revision.

## Contour-box regressions of the combined change (2026-10-19)

The D0 oracle sweeps (seeds 1 to 4, every row checked against its 256-bit
value) found three contour boxes that the line offsets and second contour
above made worse than 09b2c74. Worst binary64 error over the 24 relabellings
and the hand-written, generated, interpreted and SymJIT routes
(`tests/d0_contour_regressions.rs`, references COLLIER 1.2.9 in quadruple
precision fed the exact binary64 inputs):

| Box | 09b2c74 | line products (9a898ee) | now |
|---|---|---|---|
| three complex masses (`4 1.0 0.0 596.64 -14.72 1.0 9.0 -0.411 ...`) | 2.9e-12 | 3.9e-10 | 2.9e-12 |
| three real masses (`4 1.0 0.01348 0.01433 50.50 -0.4886 ...`) | 5.8e-13 | 1.7e-11 | 5.7e-13 |
| `p4 = 0`, `p1² = s12 = Re m1² = Re m4²`, `s23` one ulp below, `m2 = 0` | 130 times the value | NaN | 7.7e-12 |
| the same with `s23` 1e-9, 1e-4 (relative) below `p1²` | 1.7, 8.1e-6 | 4.3, 2.5e-5 | 2.1e-12, 1.3e-12 |
| relabelled, `p3²` 1e-11, 1e-7 below `s12` | 71, 1.5e-2 | 170, 5.9e-2 | 1.2e-12, 1.2e-12 |

- **The line products only for a line nearly parallel to the contour
  direction.** A pole also lies next to a line edge's root where the conic
  nearly vanishes at the line's intersection with that edge; there `a1` is of
  the size of its terms, the product's polynomial `Q(u) - ratio B(u)`
  cancels at least as much as `y - z` (by 4e5 against 1e4 in the three-mass
  boxes), and the line's term cancels against the conic's term on the same
  edge, whose root lies next to the same pole: both must keep `y - z` from
  the same rounded pole. The products are now taken only where `a1 = g +
  beta h` cancels by more than 64 (`4096 |a1|² < |g|² + |beta h|²`), and not
  on an edge whose intersection with the line lies on the conic to within
  the square root of the rounding (`2^52 |Q - ratio B|²` below the sum of the
  conic coefficients' norms; with the product there the one-ulp box lost
  1e-5 at 128 bits). The first two boxes are back to their 09b2c74 errors
  (the second bitwise on the hand-written route).
- **A pole at a Feynman-parameter corner.** In the third box the line passes
  within an ulp of the corner (1, 0) of `zero_external`'s triangle, where the
  massless `m2` makes the conic vanish: a pole lies within 4e-31 of it, at an
  `x` next to 1 that binary64 resolves only to eps (the line product, exactly
  zero there, gave `log 0`). `zero_external` now takes the reflected labelling
  `(p1, p4, p3, p2, s23, s12; m1, m0, m3, m2)` where `m2` is massless and
  the line passes nearer (1, 0) than (0, 0): the corner is then at the
  origin, where the pole keeps its relative precision. In some relabellings
  the pole then lies next to x = 0 of an edge whose line root is at the other
  end, and `r1`'s dilogarithms take the complements `1 - x` of their
  arguments from the poles (`y/z`, `y/(y - z)`, `(1 - y)/(1 - z)`,
  `(y - 1)/(y - z)`, through `sheet_exact::complement_difference`): from the
  rounded quotient a pole 7e-16 from x = 0 lost 8% of `log(1 - x)` on a
  winding sheet (1.2 times the value in 8 relabellings, at 09b2c74 too). The
  hand-written kernel forms a complement only where its dilogarithm takes
  `log(1 - x)` on a winding sheet or has an `x` that rounds to one.
- **An `a1` at the rounding level.** Where the selected contour's `a1` rounds
  to zero, the products take `2^-52 g` in its place: the offsets are
  proportional to it, and the logarithm of the common factor cancels between
  the line edges; with exact zero offsets each edge took `log 0`. Of 400
  boxes with the Gram determinant at the rounding level (seeds 201, 202 of
  `tests/generate_near_threshold_rows.py`), 112 were beyond 1e-6 on some
  route at 09b2c74 and 5 at 9a898ee (four of them O(1) on one or two
  routes); now 3, at most 5.1e-3 (below).
- **`K` guard.** A first contour whose `K` vanishes where the second's does
  not is no longer selected when it encloses a pole either (the switch above
  applied only when it encloses none); no such box was observed.
- **Checks** (against 9a898ee and 09b2c74):
  - D0 oracle sweeps, seeds 1 to 4, native, `--generated` and `--expression`
    routes: all 12 failure lists identical to both; OneLOop-unstable rows
    counted as binary64-unstable drop from 364/364/376 (seed 3) to 348 and
    from 404 to 392 (seed 4).
  - Every one of the 37,642 captured sweep rows on the four binary64 routes
    against its 256-bit value: 112 to 156 rows per route improve by more
    than 10 times (most of them three-mass boxes); the only rows worse than at
    9a898ee are the 8 relabellings of one seed-1 box (below), wrong at both
    revisions.
  - The near-threshold scans: `tests/data/near_threshold_boxes.txt` as
    before (line 73 below); seeds 2 to 6: 4 to 7 rows per route beyond 1e-11,
    at most 4.0e-11 (9a898ee: 4 to 6, 3.2e-11; such rows move by 1e-11 under
    1-ulp perturbations on either revision); seeds 101 to 105: 3 to 5, at
    most 4.4e-11 (9a898ee: 7.8e-11).
  - The 309 distinct D0 calls of NinjaRS's master streams: worst 1.9e-12
    hand-written, 1.7e-12 generated, 1.8e-12 interpreted, 5.2e-12 SymJIT
    (9a898ee: 3.6e-12, 3.7e-12, 8.8e-12, 1.7e-11). Their binary64 bits
    change on 117 (hand-written), 202 (generated), 208 (interpreted) and 214
    (SymJIT) of them against 9a898ee.
  - Retired instructions per D0 call on the 22 NinjaRS `*.masters.json`
    streams against 09b2c74: geometric mean -0.12%, at most +0.29%
    (syn-n4-r5), at least -0.63% (9a898ee: +0.21%, +0.47%, -0.06%).
  - Definitions: 100,424,043 bytes (gate 100,663,296).
- **The sweep rows.** At 9a898ee seeds 1 to 4 of `tests/oracle_sweep.py`
  gave the same failure lists as 09b2c74 on every route, but comparing every
  evaluated row (all 24 relabellings) with its 256-bit value, four boxes had
  lost accuracy against 09b2c74: the three contour boxes of the table above
  (rows 1 to 3 of `tests/data/contour_regressions.txt`) and a two-mass box,
  fixed in the next section. At 9a898ee:
  - `4 1.0 0.0 596.6380017856613 -14.719542747253795 1.0 9.0
    -0.41134991169772017 0 0 0.17229848844808035 -0.0042830249130855775 9.0
    -0.48397362776576125 1.0 -0.00011332140135807317` (three complex masses,
    contour box). It lost up to 2.3e-10 hand-written, 1.1e-10 generated,
    3.9e-10 interpreted and 2.4e-10 SymJIT, where 09b2c74 kept 2.9e-12.
    OneLOop is unstable here (its two precisions differ by 1e-10).
  - `4 1.0 0.013477197248422343 0.014331576181943157 50.495439982011966
    -0.4886017458582275 0.014331576181943155 4.0 0 0 1.1039110268637602 0
    0.014331576181943157 0 4.0 0` (three real masses, contour box). It lost
    up to 9.1e-12 hand-written, 1.7e-11 generated and interpreted, and
    5.4e-12 SymJIT, where 09b2c74 kept 5.8e-13.
  - `4 1.0 0.01545852083001525 2.25 -0.32983766127195113 0.0
    0.01545852083001525 0.015458520830015247 0.01545852083001525
    -0.00013458303781508371 0 0 2.25 0 0.01545852083001525
    -0.00013458303781508371` (contour box, several on-shell coincidences). It
    was NaN on every binary64 route, where 09b2c74 was finite but wrong by 44
    to 128 times the value. OneLOop's double precision is off by 55 times.

## Two-mass boxes with both masses at or next to their shells (2026-10-19)

- **The previous section's sheets lost binary64 accuracy where both massive
  lines sit at or next to their shells**: a real mass with an invariant
  (with a massless line) equal to m² or within about 1e-6 of it, and a
  complex mass whose real part is at or next to its shell, or two complex
  such masses. Denner and Dittmaier's mass pairs take the arguments
  `u r031 = X r14` and `u r032 = X/r14`, and their quadratic at a root
  `-1/r03k` of `m4² + Y03 x + m3² x²` is `(Y01 - Y13/r03k)(Y02 - Y23/r03k)`:
  small there, so one argument is next to one. `1 - X r14`, formed from the
  rounded argument, kept only `|1 - X r14|` relative digits, and the
  dilogarithms on wound sheets (which carry `log(1 - x)`) and their divided
  differences lost the rest. The combined candidate's sweep regression
  `4 1.0 -4.400608205605747 3.4425584296359024 0.3420655554270402
  77.52023396188905 0.16951572574727045 0.06566799954955663 0 0
  0.06566799954955661 -0.0002531972626742262 0.34206555542704026 0 0 0`
  (`p3² = m3²` to 1.6e-16, `s23 = Re m2²` to 2e-16) lost 3.3e-11 native,
  6.8e-12 interpreted and 2.8e-11 SymJIT (09b2c74: 1.2e-15); the two-mass
  verifier's row with a real mass next to its shell in two invariants (1e-6
  and 4e-8) 4.8e-9; boxes exactly on shell 4.9e-11. On 4000 random boxes of
  this class (masses 0.1 to 10, adjacent or opposite, one real or none,
  |Im m²|/Re m² from 1e-7 to 0.3, one or two invariants per mass at
  Re m² (1 + d) with d = 0 or 1e-16 to 1e-4) 9a898ee lost more than 1e-8
  in 56 (up to 0.16 of the value), and in 60 of the 3600 relabellings of 150
  more; 09b2c74's sheets (OneLOop's, right on these rows) lost at most
  2e-10.
- **The mass pairs' complements now come from the channels'.** At a root,
  `(1 - X r14)(1 - X/r14) = X² - k14 X + 1` equals
  `(1 - X r34/r13)(1 - X r24/r12) r12 r13/r23` (the quadratic is
  `(r34 X - r13)(r24 X - r12) - r23 (X² - k14 X + 1)`). The caller forms the
  channel complements from the kinematics, as before, and passes
  `r12 r13/r23` (`Y01 Y02/(Y12 m4²)`) as a twelfth argument of
  `__olo_box_two_mass_sheets` (bound to an existing symbol of box_complex's
  definitions, so no symbol is created). A complement below 1/16 is that
  product over the other complement where the other is above 1/4; otherwise
  both are formed directly, as before. Next to `r14 = 1` (the threshold
  `Y03² = 4 m3² m4²`: equal masses with a vanishing invariant between them)
  both complements are small and `r14` is known only to the square root of
  its rounding; the direct forms share that error and it cancels between the
  pairs, which the product does not (two sweep boxes there lost 1e-9 with the
  product; the direct forms keep 5e-16 and 1.6e-12, as before). The mass
  pairs then take `complement_difference` with these complements, as the
  channels do. The
  exchanged formula (a real `m4²` exactly on shell) takes the complements of
  the inverse arguments as `-w/y` with `y = X r34/r13` itself rather than
  `-w/(1 - w)`, which lost `y`'s digits for a small `y`. The one-mass box is
  unchanged.
- **Denner and Dittmaier's two formulas are still chosen on exact zeros.**
  With the complements fixed, choosing the line-exchanged formula by
  conditioning (`|Y01 Y02| m3²` against `|Y13 Y23| m4²`, as the two-mass
  verifier suggested, or without the masses) left the near-shell rows within
  4e-14 either way, but moved 38 rows of the other scans by more than ten
  times (collinear rows up to 7.5e-11 from 2.1e-12; the worst collinear row
  went from 9.6e-9 to 2.8e-9). It is not used.
- **Validation** against COLLIER 1.2.9 in quadruple precision (COLI and DD
  branches with the exact binary64 inputs, taken where they agree to 1e-10)
  on 128,349 rows: the two-mass stream's and the two-mass verifier's scans
  (every complex two-mass pattern, one real mass, narrow widths, shells,
  near shells, collinear limits, equal masses, the 24 relabellings of
  on-shell boxes, and the 3-, 4- and one-mass patterns) and the 7,624
  near-shell rows above. No row is less accurate by more than ten times
  (and 1e-13) on any route (native hand-written and generated, interpreted,
  SymJIT, native and interpreted 256 bits); 1216 rows are more accurate
  (native). No row is off by more than 1e-8 (9a898ee: 116), and the
  near-shell rows are within 1.3e-13 (9a898ee: 0.16). The remaining rows
  beyond 1e-10 are the collinear ones of the open limitations below.
  `tests/d0_two_mass_width.rs` now also checks the 9 boxes of
  `tests/data/two_mass_shell.txt` (the sweep regression, the verifier's row,
  an exactly on-shell box and the worst random boxes by pattern; references
  from the formula in mpmath at 90 digits, which COLLIER's COLI branch
  reproduces to 1e-17) in all 24 relabellings: binary64 routes within 1e-12
  (measured 6.8e-16), DoubleFloat within 1e-15, and native and interpreted
  `Float` at 128 and 256 bits within their tracked precision (they keep at
  least 80 and 208 bits there; see the open limitations). 9a898ee fails both
  tests (876 failures). Every real-mass and one-mass row, the three- and
  four-mass rows and the NinjaRS master rows are bitwise unchanged on every
  route, including 256 bits: 0-, 1- and 2-mass real rows against 09b2c74,
  complex one-mass rows against bd619ca, all of them against 9a898ee
  (37,642 captured sweep rows and 29,366 curated ones, their complex
  two-mass rows aside). The D0 oracle sweeps
  (`tests/oracle_sweep.py`, seeds 1-4, native, generated and expression
  routes) give the summaries and failure lists of 9a898ee; against its
  captured rows, no row regresses and the 12 relabellings of the sweep
  regression above are back to 4.5e-16. Retired instructions per binary64
  D0 call (perf stat, difference of 2001 and 1 repetitions): the 22 recorded
  NinjaRS master streams within 0.0006% of 9a898ee (against 09b2c74 still
  +0.21% geometric mean, at most +0.471%); complex two-mass boxes +4.5% to
  +4.7% against 9a898ee on the samples of this section (1000 random 00MM and
  0M0M boxes each, 1000 each of 00MR and 0R0M with one real mass, the 4,000
  random near-shell boxes and the 3,600 relabelled on-shell rows above;
  29,535 to 30,858 for 00MM), real two-mass boxes (1000 random) unchanged.
  The
  definitions take 100,349,287 bytes (9a898ee: 100,402,551; the gate of
  `tests/expression_construction.rs` is 100,663,296).
- **NinjaRS**: of the 44 complex two-mass sub-boxes of the six V0 edge cases
  (`edge-low-rank-n5-r0` to `-r3`, `edge-many-n7-r6`, `-r7`), 10 change in
  their last bits on the native route (6 generated, 4 interpreted); they are
  within 7e-15 of COLLIER (1.9e-14 before). The re-bless that 9a898ee needs
  covers them.

## Open D0 limitations after these changes (2026-10-19)

Found while validating the sections above together. Unless a bullet
says otherwise, each one is the same at 09b2c74. References are COLLIER 1.2.9
in quadruple precision (COLI and DD branches, fed the exact binary64 inputs)
and the 256-bit expressions, which agree on every row quoted here unless a
bullet says otherwise.

- **Relabelled one-mass boxes next to a collinear limit and to
  `r13 r24 = r12 r34` at once.** Take 1500 complex one-mass rows (727 in
  the generated labelling, 773 relabelled) with `p1²` from 1e-12 to 1e-4,
  `s12` within 1e-8 to 1e-2 of `p2²`, and `s23` within the same range of
  `p4²`. When the small invariant
  lands in the `p2²` or `s12` slot of the one-mass formula, binary64 loses up
  to 3.4e-4 (native, generated and interpreted; SymJIT 1.5e-4). The labelling
  that keeps it in `p1²` stays within 7.6e-11 (SymJIT 9.0e-11). 09b2c74 loses
  the same 3.4e-4 on that box. With a real mass one row is wrong outright:
  `4 1.0 1.3735597590592743 1.3735605975942236 -3.685517126287593
  -3.685509990953621 3.784390017645558e-12 0.6002365033125651 0 0 0 0 0 0
  0.2512196489044931 0.0` gives -2.107 + 79.14i in binary64 on the native
  (hand-written and generated) and interpreted routes, for -5.6547 - 0.5810i.
  SymJIT is within 1.4e-6.
- **A real negative `m²` is not the `-i0` limit.** The domain check accepts
  it, and it takes the real-mass expressions. For `4 1.0 1.5939214530073218
  -2.5916131266951306 0.02997860964984738 -15.99681513243264
  0.06920134388722839 -0.16447700103649665 0 0 -0.08665390287240995 0 0 0 0
  0`, every route gives 1.9204 + 1.9252i. With Im m² from -1e-8 down to
  -1e-20, the complex-mass expressions give 0.9341 + 1.9252i, and COLLIER
  agrees at Im m² = -1.4e-12 (both branches) and with its DD branch at the
  real mass; its COLI branch gives 1.4273 + 0.5241i there. The
  real-mass value is bitwise that of 09b2c74.
- **SymJIT for |Im m²| below about 1e-150 of Re m².** Take 400 one-mass boxes
  with |Im m²|/Re m² from 1e-40 to 1e-300. The native and interpreted routes
  stay within 7.2e-15 of the real-mass value. SymJIT departs from it by more
  than 1e-6 in 9 boxes at 1e-150, in 149 at 1e-154 and in 330 from 1e-170
  down, by up to 37 times the value. 09b2c74's SymJIT does the same (250 from
  1e-170 down). Such widths are not physical.
- **Opposite-mass two-mass boxes next to a collinear limit** (complex 0M0M,
  M0M0) lose up to 1e-8 in binary64. On 500 such boxes per pattern: 0M0M
  6.2e-9 (native, generated) and 9.6e-9 (interpreted, SymJIT); M0M0 2.4e-9
  to 3.3e-9. The adjacent patterns keep 1.7e-10. The previous expressions
  lost up to 6e-7 on the rows they got right. The loss comes from the
  collinear machinery, not from the sheets.
- **Contour box next to a vanishing Gram determinant.** Line 73 of
  `tests/data/near_threshold_boxes.txt` keeps 1.3e-11 (interpreted) to
  3.2e-10 (hand-written; OneLOop's double precision 5.5e-11), as at 09b2c74.
  One `tfun` has a nearly double contour direction; in the others, nearly
  coincident poles sit next to an edge's end, and their quadratic `s3` terms
  (about 5e3 each) cancel to about 3e2. Of 1,000 more such boxes (seeds 2 to
  6 of `tests/generate_near_threshold_rows.py`), 4 to 7 per route remain
  beyond 1e-11, at most 4.0e-11 (09b2c74: 110 to 117 per route), and of
  1,000 on seeds 101 to 105 3 to 5, at most 4.4e-11. Such rows move by
  1e-11 under 1-ulp perturbations of their invariants. One cause is a
  cancellation of about 1000 between the three `tfun`; another, a pole next
  to an integration endpoint, now takes the complements `1 - x` of the `r1`
  dilogarithms' arguments where their sheets wind (section above).
- **Contour box with the Gram determinant at the rounding level.** Where the
  selected contour's `a1 = g + beta h` rounds to zero and the masses are
  equal and complex, the poles are real in binary64 and lose the imaginary
  parts of order `a1 Im m²` that choose their lips: `4 1.0 2.7819050504399936
  2.3889525049365044 2.388952504936503 105.83204530964998 51.91801797559696
  9.555810019746128` with `m² = 2.2406833551152996 - 0.0005207803841710876i`
  four times keeps 5.1e-3 on the hand-written and interpreted routes in every
  relabelling (the generated code and SymJIT 8e-14; 0.34 at 09b2c74).
  Two more of those 400 boxes keep 9.7e-4 and 4.0e-3 on one route. The box is
  the known-failure fixture of `tests/d0_contour_regressions.rs`.
- **A conic singular at a massless corner.** `4 1.0 0.1113435943976259
  108.23851883958666 0.11134359439762591 -0.012515862621692797 9.0 0.0 0 0
  0.11134359439762591 0 9.0 0 0.11134359439762591 0` (`p1²` one ulp below
  `m2² = m4² = p3²`, massless `m1`): in the contour box's labellings the
  conic's value and gradient vanish at the massless corner to within
  rounding. Every binary64 route is wrong in every relabelling, as at
  09b2c74 (then NaN in 8, by 1.4e-2 in 4 and by 0.3 in 12; now NaN in 16 and
  1.4e-2 in 8, so 8 relabellings of the seed-1 sweep went from 0.3 to NaN).
  Known-failure fixture of `tests/d0_contour_regressions.rs`.
- **128 bits next to a corner.** At 128 bits the boxes of
  `tests/d0_contour_regressions.rs` with a pole within `t²` of a corner keep
  `2^-128/t²` of its offset: the relabellings agree only to 1e-16 to 1e-29,
  and the one-ulp box is within 1.2e-11 (native) to 3.8e-11 (interpreted) in
  16 relabellings (20 at 09b2c74). Further from on shell, Float's tracked
  precision of the poles' small components costs the labelling with the
  corner at the origin up to 2.5e-7 at 128 bits (2.4e-7 at 09b2c74; `p3²`
  3e-16 below `s12`), which the reflection now also takes for 8 of the
  relabellings that kept 1e-15 at 128 bits at 09b2c74 (and gave O(1) errors
  in binary64). 256 bits keep 1e-15 (and agree to 1e-30).
- **Native `Float` at the two-mass threshold `r14 = 1`.** In `4 1.0
  -7.158359142869336 -0.006328911565103391 0 0 0 0.01281606414601131 0.25
  -4.478335170728946e-05 0 0 0.25 -4.478335170728946e-05 0 0` (equal complex
  masses on opposite lines, a vanishing invariant between them, two more
  vanishing invariants) the native 256-bit route gives 3.072e4 + 19i,
  9.728e3 - 2.307e7i or 1.9456e4 - 4.614e7i, numbers of a few bits, for
  8.6969 + 5.8296i in 12 of the 24 relabellings, as at 9a898ee (09b2c74:
  2.8672e4 + 19i and the like). Binary64 on every route and the interpreted
  256-bit expressions are right (5e-16 and 7e-17 against COLLIER's COLI
  branch, whose DD branch stops with error -10, and the formula in mpmath).
- **`Float` keeps the precision it tracks.** Symbolica's `Float` rounds the
  difference of nearly equal operands to the bits it considers significant
  (as in `tests/permutation_regressions.rs`). Next to a shell an invariant
  and a mass differ in their last binary64 bits: in
  `tests/data/two_mass_shell.txt` the native and interpreted results keep 80
  of 128 bits and 208 of 256 (errors 1.5e-23 and 1e-61), within that tracked
  precision.

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
