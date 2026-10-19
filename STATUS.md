# Native OneLOop status

This page preserves the September 8 development snapshot. See the
[Symbolica main migration report](SYMBOLICA_3_MIGRATION.md) for the current
dependency setup and validation.

## Float products and factors of ±1 (2026-10-19)

The `Float` routes (native code and interpreted expressions, as NinjaRS's
arbitrary-precision masters use them) were wrong on rows where binary64 was
right. In the D0 oracle sweeps (`tests/oracle_sweep.py`, seeds 1 to 4, all 24
relabellings: 37,642 rows, 25,645 of them with OneLOop's two precisions within
1e-9), the following numbers of rows were more than 1e-8 off OneLOop's
quadruple precision while binary64 was within 1e-8:

| Route | 128 bits | 256 bits | 512 bits | 1024 bits |
|---|---|---|---|---|
| native, before | 82 | 62 | 56 | 59 |
| native, now | 12 | 0 | 0 | 5 |
| interpreted, before | 96 | 44 | 68 | 46 |
| interpreted, now | 12 | 0 | 4 | 4 |

The errors ranged from 5e-2 to 4e6 of the value, for example 3.07e4 + 19i
for 8.70 + 5.83i at 256 bits. At 256 bits 09b2c74 and the combined branch
fail on the same rows.

In the C0 sweeps (A0, B0, dB0 and C0, seeds 1 to 4: 17,555 rows, 17,198
stable) the counts were:

| Route | 128 bits | 256 bits | 512 bits | 1024 bits |
|---|---|---|---|---|
| native, before | 62 | 46 | 33 | 52 |
| native, now | 23 | 7 | 6 | 12 |
| interpreted, before | 62 | 33 | 33 | 39 |
| interpreted, now | 20 | 12 | 8 | 5 |

A0, B0 and dB0 had none before or after.

- **Cause.** Symbolica's numerica gives a `Float` its own precision, and
  values carry different ones: a square root has one bit more than its
  argument, and additions track cancellations. Two of its operations then
  broke identities that binary64 keeps exactly.
  - A product by reference first rounded the more precise factor to the
    other's precision. So `a b` and `b a` could differ in the last bit, and
    `z conj(z)` kept an imaginary residue.
  - Evaluators form `a - b` as `a + (-1) b`, with the `-1` at the working
    precision. The product therefore rounded a `b` that carried more bits.
    For about half of all square roots, `x - |x|` (the sign test of
    `sheet_exact::sign_nonnegative`) and `x - conj(x)` (the imaginary part)
    were then one-bit residues instead of zero.

  Exact sign and zero tests branched on these residues. For example, a tie
  of `|b + d|` and `|b - d|` in `ordered_negative_quadratic_roots` came out
  `-2.8e-76` at 256 bits and swapped the roots of the two-mass box above. In
  other rows they flipped the lip of a logarithm or the sheet of a
  dilogarithm, or let a few-bit value dominate the result.
- **Fix, in Symbolica.** Products and quotients are rounded once, to the
  smaller precision, and a factor of exactly ±1 keeps the other factor
  exact. This is symbolica-dev/symbolica#97. The pin is now `9cdc2949` (branch
  `claude/oneloopmaster-integration-float-products` of ValentinHirschi/symbolica:
  `e679e8c1` plus that commit). Binary64 is bitwise unchanged on every route
  (hand-written, generated, interpreted, SymJIT): on the 37,642 D0 and
  17,555 A0 to C0 sweep rows. DoubleFloat evaluates its dilogarithm in
  `Float` at 160 bits. It is bitwise unchanged on those rows and on the
  7,645 NinjaRS calls below (62,842 rows in all).
- **Regression test.** `tests/float_exact_products.rs` (`tests/data/float_exact_products.txt`)
  covers 7 D0 and 4 C0 rows from the sweeps in all their relabellings.
  - References: Denner and Dittmaier's two-mass box, or the C0
    Feynman-parameter integral, in mpmath at 40 digits; otherwise COLLIER
    1.2.9 in quadruple precision, fed the exact decimal expansions of the
    inputs.
  - Routes: binary64 (hand-written and generated native, interpreted,
    SymJIT) and DoubleFloat, plus native and interpreted `Float` at 128,
    192, 256, 320 and 512 bits.
  - It fails with the previous pin in 605 checks on every row, and passes
    with the new one.
- **NinjaRS.** Take all 3- and 4-propagator subsets of every lane of the
  rank n+2 suite (cases, controls, identity, scan) and of the final review
  (cases, cases8, cases-eucl, cases-scan, dotcases). That gives 7,645
  distinct C0/D0 calls.
  - With the previous pin, all 8 `Float` lanes (native and interpreted, 128
    to 1024 bits) agree with each other to 1e-8 on every call. None of these
    calls is affected.
  - The fix moves their values by at most 2.9e-32 at 128 bits and 4.7e-71 at
    256 bits relative.
  - NinjaRS pins Symbolica separately (`alphal00p/symbolica` `e679e8c1`).
    Its own `Float` masters change only after its pin moves to a commit with
    #97.

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
  64 ([PRECISION.md](PRECISION.md)).
- **NaN at an exactly vanishing Gram determinant in binary64.** `K` of one
  contour is `b (g + beta' h)` of the other; where it rounded to zero the
  native binary64 code (hand-written and generated) returned NaN. The second
  contour is used then.
- **Scan.** `tests/data/near_threshold_boxes.txt` (200 random physical boxes,
  `tests/generate_near_threshold_rows.py`, OneLOop 3.6 quadruple-precision
  references): every binary64 route within 5e-12 on 199 rows (3.2e-12 at
  most), where 34 were beyond it and 22 beyond 1e-9 before (OneLOop's double
  precision: 32 and 16). On 1,000 more such boxes (seeds 2 to 6) 4 to 6 per
  route remain beyond 1e-11 (at most 3.2e-11), against 110 to 117 before.
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
  |Im m²| from 1e-3 to 0.5 of Re m² wrong by 0.6 to 17 of the value; 2 to
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
  NinjaRS master streams within 0.002%; real two-mass boxes +0.4% (23,673 to
  23,763); complex two-mass boxes -1.2% (00MM) and -0.1% (0M0M), with one
  real mass -2.3%; the on-shell boxes with one real mass, in all
  relabellings, +14% (19,938 to 22,810).
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

## Open D0 limitations after these changes (2026-10-19)

Found while validating the three sections above together. Unless a bullet
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
  `tests/data/near_threshold_boxes.txt` keeps 1.1e-10 to 3.2e-10 (OneLOop's
  double precision 5.5e-11), as before. One `tfun` has a nearly double
  contour direction; in the others, nearly coincident poles sit next to an
  edge's end, and their quadratic `s3` terms (about 5e3 each) cancel to
  about 3e2. Of 1,000 more such boxes (seeds 2 to 6 of
  `tests/generate_near_threshold_rows.py`), 4 per route remain beyond 1e-11
  (6 on SymJIT), at most 3.2e-11 (09b2c74: 110 to 117 per route). One cause
  is a cancellation of about 1000 between the three `tfun`. Another is a
  pole next to an integration endpoint: there the `r1` dilogarithm arguments'
  `1 - x` would need complements, as in `sheet_exact::complement_difference`.
- **Regressions against 09b2c74 in the D0 oracle sweeps.** Seeds 1 to 4 of
  `tests/oracle_sweep.py` give the same failure lists as 09b2c74 on every
  route. Comparing every evaluated row (all 24 relabellings) with the 256-bit
  value, four boxes lose accuracy:
  - `4 1.0 0.0 596.6380017856613 -14.719542747253795 1.0 9.0
    -0.41134991169772017 0 0 0.17229848844808035 -0.0042830249130855775 9.0
    -0.48397362776576125 1.0 -0.00011332140135807317` (three complex masses,
    contour box). It loses up to 2.3e-10 hand-written, 1.1e-10 generated,
    3.9e-10 interpreted and 2.4e-10 SymJIT, where 09b2c74 kept 2.9e-12.
    OneLOop is unstable here (its two precisions differ by 1e-10).
  - `4 1.0 0.013477197248422343 0.014331576181943157 50.495439982011966
    -0.4886017458582275 0.014331576181943155 4.0 0 0 1.1039110268637602 0
    0.014331576181943157 0 4.0 0` (three real masses, contour box). It loses
    up to 9.1e-12 hand-written, 1.7e-11 generated and interpreted, and
    5.4e-12 SymJIT, where 09b2c74 kept 5.8e-13.
  - `4 1.0 0.01545852083001525 2.25 -0.32983766127195113 0.0
    0.01545852083001525 0.015458520830015247 0.01545852083001525
    -0.00013458303781508371 0 0 2.25 0 0.01545852083001525
    -0.00013458303781508371` (contour box, several on-shell coincidences). It
    is NaN on every binary64 route, where 09b2c74 was finite but wrong by 44
    to 128 times the value. OneLOop's double precision is off by 55 times.
  - `4 1.0 -4.400608205605747 3.4425584296359024 0.3420655554270402
    77.52023396188905 0.16951572574727045 0.06566799954955663 0 0
    0.06566799954955661 -0.0002531972626742262 0.34206555542704026 0 0 0`
    (two adjacent masses, one complex, both next to on shell). It loses
    3.3e-11 on the native routes, 2.8e-11 on SymJIT and 6.8e-12 interpreted,
    where 09b2c74 kept 1.2e-15.

  The first three come from the contour box's line offsets and second
  contour. The last comes from the two-mass sheets.
- **`Float` rows that the numerica fix leaves wrong** (Float products and
  factors of ±1, above). These are D0 sweep rows where OneLOop is stable and
  binary64 is right. Most other sign and zero tests on cancelled values
  still branch on rounding residues. The difference is that these residues
  come from expressions that are equal only analytically.
  - `4 1.0 0.0 12.25 12.25 -0.04070524207988602 -0.6390609176495762 0 0 0
    2.25 -0.0003109069463992361 4.0 0 2.25 -0.0003109069463992361`: in 12
    of its 24 relabellings, both `Float` routes at 128 bits give 8.0 + 8.5i
    (4 tracked bits) for 8.238 + 8.283i. From 192 bits up the routes are
    right. Binary64 loses 7.8e-12 here.
  - `4 1.0 0.01230693334712481 0 0 0 0 2.5453032460105445 0 0 0.06371291035718132
    0 0.06371291035718132 0 4.0 0`, in 4 relabellings: off by 1.0 relative
    (interpreted, 512 bits) and by 5.6e4 and 1.1e5 (native and interpreted,
    1024 bits).
  - `4 1.0 -1.3135178794812885 4.616075471676815 1.939351244333556 0.0
    -0.37504160373662826 0.04613395493643435 0.011533488734108588 0
    0.011533488734108588 0 2.25 0 0.011533488734108588 0`, one relabelling:
    native 1024 bits off by 6.4.
  - Unrelated to precision: 96 OneLOop-unstable D0 sweep rows have high-precision
    routes that disagree with each other. These are degenerate (on-shell, IR or
    threshold) boxes where binary64 and OneLOop fail too, and some grow with the
    precision. For example, `4 1.0 0.0 5.1835752237236346 2.314559146961827
    0.054988940641060316 3.083060962167901 -7.793554644603385 0 0
    3.083060962167901 0 3.083060962167901 0 0.054988940641060316 0` gives
    -1.6e9 at 128 bits and 4.0e13 at 4096 bits.

## Open C0 limitations (2026-10-19)

The C0 oracle sweeps cover A0, B0, dB0 and C0, seeds 1 to 4, 17,555 rows.
- **Equal masses with on-shell invariants and one vanishing invariant.** Of
  the eight route and precision pairs (native and interpreted, 128 to 1024
  bits), three triangles are wrong on all eight and a fourth on four, by
  0.07 to 1.5 relative. Binary64 is right to 1.5e-15, and so is OneLOop.
  The four triangles are:
  - `3 1.0 4.0 4.0 0.0 4.0 0 4.0 0 4.0 0`: -0.1406 + 0.0039i at 256 bits
    (4 tracked bits) for -0.15115;
  - `3 1.0 0.0 3.2461522396961473 3.2461522396961473` with three masses
    3.2461522396961473;
  - `3 1.0 2.25 2.25 0.0 0.19617369498967158 0 2.25 0 0.19617369498967158 0`;
  - `3 1.0 28.40228792397598 28.40228792397598 0.0` with three masses
    28.40228792397598 - 0.10470706900097601i.

  The roots of the quadratic are analytically real but formed in complex
  arithmetic. For example, the product of `0.25 + 0.433i` and `0.5 - 0.866i`
  keeps an imaginary residue of 2e-78 at 256 bits, and sign tests branch
  on it.
- **Two-mass rows at single precisions.** Five more triangles are wrong on
  one to three of the eight pairs, by 1e-7 to 0.19 relative.
- **Interpreted `Float`, C0 at vanishing invariants.** At 112, 128 and 224
  bits, `3 1.0 0 0 0 9.0 -0.1667104301703281 21.592408834323358 0 0 0` and its
  rotations keep 14 to 19 tracked bits and are wrong by 5.3e-3 to 1.2e-2.
  The native route is right. The ignored reproducer is
  `open_interpreted_float_defect` in `tests/float_exact_products.rs`.
- **Every route, binary64 included, with the masses reflected.** The same
  triangle with the massless line between the complex and the real mass is
  wrong even though, with vanishing invariants, the value cannot depend on
  the order of the masses. For `3 1.0 0 0 0 9.0 -0.1667104301703281 0 0
  21.592408834323358 0`:

  | Evaluator | Value |
  |---|---|
  | hand-written and generated native | -2.94e-5 + 1.49e-3i |
  | interpreted | 0.0191 - 2.2e-6i |
  | SymJIT | -0.0239 + 2.5e-4i |
  | DoubleFloat | -0.234 - 0.512i |
  | OneLOop quadruple precision | -0.2859 - 0.1114i |
  | correct (mpmath Feynman-parameter integral, 40 digits) | -0.069489 - 5.509e-4i |

  OneLOop's two precisions disagree on this row, so the sweep counts it
  as unstable.

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
