# Native OneLOop status

This page preserves the September 8 development snapshot. See the
[Symbolica main migration report](SYMBOLICA_3_MIGRATION.md) for the current
dependency setup and validation.

## The final candidate after the fresh sweep (2026-10-21)

`claude/olm-final-3`: the candidate of the next section (59b96e5) with three
contour-box fixes for the silent errors of the second verifier's fresh
sweep, regression tests for them, the C0 SymJIT asset rebuilt, and one
regeneration (`generated_d0.rs` 2,592,725 bytes, the D0 asset 1,558,784
bytes, the C0 asset 85,547 bytes; the definitions take 95,753,728 bytes
against the gate's 100,663,296). fs is that sweep (the verifier's own
generator, seed 20261021, 36,840 rows: 12,000 random boxes over all 81 mass
patterns with widths 1e-14 to 0.7 of the mass and near-shell, threshold,
collinear and corner legs; 2,000 physical 2 -> 2 boxes from real momenta;
3,000 boxes with `|Gram|/max|inv|³` from 1e-4 to an exact root; 1,000
planar boxes with a transverse kick; 420 special and 240 soft boxes in all
24 relabellings; 3,000 C0), rl the 24 relabellings of its 30 finding boxes;
its references are COLLIER 1.2.9 in quadruple precision (COLI where its DD
branch agrees to 1e-13, else where the 256-bit values of both builds agree
to 1e-30 and COLI with them to 1e-13). The other sets, "worse", "silent"
and the routes are as in the next section; new NaN counts every row,
referenced or not (SW: 24 of its 32 per route have no reference).

- **The other chart's K a rounding residue** (`box_complex::tfun_chart`,
  `tfun_checked`). Where the selected contour is flat and the other's chart
  disagrees with it, the other chart was taken as regular. Next to an exact
  vanishing Gram determinant its K is itself a rounding residue: fs 16695
  (0CRR) took it in one `tfun` and was wrong by 5.3e-6 on every binary64
  route in all 24 relabellings; its selected chart is right to 1.4e-9 and
  moves by 1.3e-9 with its line, the other is right to 8e-7 and moves by
  2e-6 with its K at the rounding level. Where the other chart's K is below
  2^-52 of its terms, the selected chart is now kept where K at its
  rounding level moves the other more than 2^10 times as much as the
  selected moves with its line, and must then keep 20 bits under that move
  (NaN otherwise). fs 16695 is back at 09b2c74's 1.0e-7 (hand-written) and
  2.5e-8 (the other routes), and fs 15250 (NaN on every route at 59b96e5)
  at 09b2c74's 6.3e-9 to 6.8e-9.
- **Unscaled line products probed** (`box_complex::products_kept`,
  `tfun_edge_units`). A line product more than 2^10 times `a1` from its
  root no longer drops every product of its chart; it is probed as a
  product that disagrees with its pole (at least 128 units), so the chart
  takes `y - z` where the products move it more under the pole's shift.
  Reflected, fs 15926 (R0CR, `|Gram|` 1e-13) was wrong by 4.9e-6 to 1.4e-5
  and keeps 2.4e-14; fs 15326 goes from 8.0e-6 to 8.6e-13 (hand-written;
  interpreted 2.6e-5 to 2.8e-6), fs 14025 from 6.2e-5 to 8.6e-10 on every
  route; v1 8793 and 18833, for which the old rule was made, keep 59b96e5's
  values (09b2c74's).
- **Zero-external orientation** (`box_complex::zero_external_suspicious`,
  `zero_external_state`, `tfun_flat_state`, `vertex::reflection_unresolved`).
  The reflection is skipped where the reflected chart selects a flat contour
  without a regular alternative (in `Float`, which takes the selected
  contour unchecked, any flat contour) and the original chart selects no
  flat contour. fs 15189 (CRR0, `|Gram|` 1e-15, mu² 0.103) was wrong by
  2.9e-2 (hand-written) or NaN reflected, and keeps 9.7e-11 on every route
  in the 12 relabellings that reach the reflection test, as 09b2c74; fs
  15244 is back at 09b2c74's 1.7e-13 to 2.6e-13 (59b96e5 2.1e-12 to
  4.3e-12) and no relabelling of it is worse than 09b2c74; fs 15926 keeps 09b2c74's 6.4e-10 in the 12 relabellings that
  are no longer reflected. At 128 bits fs 15486 (C0CC, exact vanishing Gram
  determinant) kept 20 tracked bits reflected (wrong by 1.2e-6 to 1.9e-5)
  and keeps 64 (3.3e-16, as 09b2c74) in the 12 relabellings that reflect,
  and v1 7952 (F-R1 below) goes from 2.1e-11 at 40 tracked bits to 3.3e-18
  at 73.
- **Symbol order.** The states are formed in place, not as registered
  functions: two new function symbols renumbered every later symbol, and
  the interpreted and SymJIT routes then ordered, and rounded, sums
  elsewhere differently (fs 7795 interpreted 9e-10 to 4.2e-8, SymJIT on four
  relabellings of SW 22062 from 1.5e-16 to 4.6e-10, v1 2601 SymJIT 3.6e-2
  to 1.5).
- **C0 SymJIT asset.** Rebuilt for the first time since PR #9 changed C0's
  real sheet signs (85,547 bytes, `7a7bd442...`; it was 102,902,
  `f8c7109a...`, built before PR #9). It keeps the binary64 SymJIT values
  of the 10,668 C0 rows of the two verifiers' sweeps and the NinjaRS master
  calls (MS 534, NR 4,134, v1 3,000, fs 3,000) and of the 9,516 C0 rows of
  the oracle sweeps (seeds 1 to 4), bit for bit but for the sign of one NaN
  (fs 34317).
- **Tests**: `tests/data/contour_fresh_sweep.txt` with the five boxes above
  (fs 16695, 15926 and 15250 in all 24 relabellings on every binary64 route;
  fs 15189 in the 12 relabellings that reach the reflection test; fs 15486 in
  binary64 in all 24 and at 128 bits, native and interpreted, in the 12 that
  keep their orientation); each of the three tests fails at 59b96e5.

### Results (2026-10-21)

Rows worse / better / silent / new NaN against 09b2c74, as in the next
section:

| set | hw | gen | expr | jit | dd |
|---|---|---|---|---|---|
| fs | 15/4164/0/8 | 15/4163/1/8 | 12/4144/1/6 | 12/4167/0/8 | 1/3364/0/2 |
| rl | 267/72/0/152 | 295/68/24/164 | 219/64/24/124 | 219/60/0/168 | 11/48/0/24 |
| v1 | 5/566/1/5 | 3/586/2/3 | 5/572/1/4 | 10/597/2/6 | 0/105/0/0 |
| CT | 9/673/0/9 | 11/678/0/9 | 11/677/0/9 | 13/676/0/9 | 1/42/0/1 |
| MS | 0/1/0/0 | 0/1/0/0 | 0/3/0/0 | 0/1/0/0 | 0/0/0/0 |
| NR | 0/0/0/0 | 0/0/0/0 | 0/2/0/0 | 0/1/0/0 | 0/0/0/0 |
| g1 | 0/133/0/0 | 1/142/0/0 | 1/150/0/0 | 0/148/0/0 | 0/61/0/0 |
| TM | 2/2133/0/0 | 2/2133/0/0 | 1/2132/0/0 | 1/2119/0/0 | 0/1844/0/0 |
| SW | 8/1420/0/32 | 8/1436/0/32 | 8/1400/0/32 | 8/1400/0/40 | 0/1300/0/24 |
| TD | 0/1096/0/0 | 0/1120/0/0 | 0/1088/0/0 | 0/1084/0/0 | 0/12/0/0 |

(TM better counts with the second verifier's reference rule; the next
section's 2122/2122/2121/2106 used the first verifier's.)

- **Against 59b96e5**: no row of any set is more than ten times worse and
  above 1e-12, or newly NaN, on any binary64 route, in DoubleFloat or at
  `Float` 128 and 256 bits. Better on fs 7/7/7/7 (hw/gen/expr/jit), rl
  96/96/120/96 and at 128 bits fs 1/1, rl 12/12 and v1 1/1. One silent
  route entry: fs 16941 interpreted, NaN at 59b96e5 and now 5.87e-6, the
  value of 09b2c74 (4.7e-6 to 6.0e-6 on the other routes there), within 10
  kappa_1ulp (1.36e-6; class J), in all 24 relabellings; the other routes
  stay NaN.
- **The fresh sweep's findings**: no silent route entry outside class J or
  J* is left. Its three silent boxes are fixed (above; fs 15189 in the 12
  relabellings that 09b2c74 got right). Of its 14 other rows outside the
  accepted classes, fs 15244 and 15250 are fixed; the other 12 rows (25
  binary64 route entries: 20 finite, 5 NaN) are listed below; there were
  45 route entries on 17 rows at 59b96e5. fs 15744 (J*) improves on SymJIT
  (2.2e-7, at 59b96e5 1.6e-6) and keeps 1.55e-6 on the generated and
  interpreted routes.
- **`Float`, rows worse than 09b2c74** (native/interpreted): fs 3/1 at 128
  bits (fs 14026 and 15340 below, and fs 1255, which is NaN where 09b2c74 was
  wrong) and 1/1 at 256 (fs 1255), none finite and wrong by more than 1e-6
  where 09b2c74 was right or non-finite; rl 59/11 at 128 and 11/11 at 256.
  fs 15486 at 128 bits keeps 1.15e-6 to 1.9e-5 at 20 tracked bits in the 12
  relabellings that put the box into the reflected orientation directly,
  where 09b2c74 was wrong by 3.7e-2 to 1.4e6 at 2 to 11 tracked bits.
- **Routes that still decide differently** (on these sets): in the 12
  relabellings of fs 15189 that put the box into the unreflectable
  orientation directly, the flat contour's line test passes on the
  hand-written route in 8 (2.9e-2) and on the generated route in 4
  (2.8e-2) and fails everywhere else (NaN), where 09b2c74 was wrong by
  3.4e-2 to 0.23 on every route. A 13-bit line test makes those values NaN
  on the hand-written route (and v1 2057, 5112 and 7159, wrong by 1.1e-3 to
  0.15 at 59b96e5), but also makes row 3 of
  `tests/data/contour_unresolved.txt` NaN on the generated and interpreted
  routes (3.7e-14 and 7.3e-14 now). fs 16941: interpreted 5.87e-6 (J), the
  other routes NaN.

### Gates (2026-10-21)

- **Full suite** (`cargo test --release --locked --features
  generated-evaluators --no-fail-fast -- --test-threads=1`, with the
  checkout's local root `Cargo.lock`, which is untracked as at 09b2c74):
  267 passed, 0 failed, 9 ignored in 62 test binaries.
- **Oracle sweeps** (seeds 1 to 4, C0 and D0, native, `--generated` and
  `--expression`, 128, 256 and 512 bits): all 36 failure lists equal
  59b96e5's, and differ from 09b2c74's only by the two boxes documented in
  the next section.
- **Instructions per binary64 D0 call** (perf stat instructions:u, 2001
  minus 1 repetitions, the 22 NinjaRS master streams, against 09b2c74):
  geometric mean -0.303%, at most +0.444% (syn-n6-r7; syn-n6-r5 +0.340%,
  syn-n6-r6 +0.330%), within the 0.5% gate. The new tests run only where the
  existing chart tests already hold.
- **Definitions** 95,753,728 bytes (gate 100,663,296; headroom 4,909,568);
  `tests/full_expansion.rs`'s budget holds.
- **Regeneration**: a second `generate_native`, `rebuild_evaluators D0` and
  `rebuild_evaluators C0` from the committed sources gives identical SHA-256
  for every asset and generated file; `cargo metadata --locked` passes in
  the three build roots.

### Accepted regressions against 09b2c74 (2026-10-21)

Every class of the next section stays accepted. The fresh sweep adds these
rows (inputs as in `rows/fs.rows`; errors now / 09b2c74 / 59b96e5, which
equals f2adce4 unless given; kappa_1ulp as defined below; none is finite
and wrong by more than 1e-6). They come from commits before 59b96e5, mostly
before f2adce4, that the first verifier's sets did not exercise:

- **Reflected zero-external charts** (14dd171): fs 1531, every binary64
  route, 6.42e-11 / 1.1e-13 to 2.2e-13 / 6.42e-11 (kappa 1.1e-13), in 12
  relabellings: `4 1.0 281.8956276731394 -70.80559947833059
  0.024785543819142837 0.0 -5.07823123522089 -1.2811539927651582e-08
  281.8956255159561 -8.003183925535893e-06 0.0 0.0 0.024785543819142837 0.0
  0.08416889009264034 -4.526136808732091e-05`; fs 3566, every route, 8.04e-11
  / 2.0e-12 / 8.04e-11 (kappa 9.8e-14, 2.1e-16 without on-shell-breaking
  moves), in 12: `4 1.0 13.17243136110834 0.868149785640474
  2.0820213756049387e-14 0.0 507.626265255123 2.125745649397817e-12
  13.17243136110834 0.0 0.0 0.0 67.18916795391507 -2.7507402731017164e-12
  0.022546450469187566 -1.1750917351676846e-15`. Unreflected they keep
  1.1e-13 and 2.1e-12. The reflection tests hold there (a selected pole
  within 1/3.4e3 and 1/5.9e2 of its size from the corner's vertex), and the
  reflected chart's own vertex test is no worse than the original's, so the
  vertex quantities do not separate these rows from those the reflection
  fixes; vetoing by them undoes 12 rows the reflection improves 10 to 1e3
  times (fs 15728, 14832, 16997) and keeps these.
- **Two-mass sheets** (a5da789): fs 6653, hw/gen/expr, 5.48e-7 / 1.94e-8 /
  5.48e-7 (SymJIT 8.9e-8; kappa 4.6e-8, 11.9 kappa), in all 24:
  `4 1.0 0.0 0.009553960128414992 0.11344809090507872 0.11344809090507897
  0.7741832343245684 0.1888466722914802 0.0 0.0 0.009553960128414992
  -5.059779680935713e-11 0.0 0.0 0.11344809090507872 0.0`; fs 11533,
  interpreted, 1.71e-10 / 1.35e-11 / 1.71e-10 (kappa 1.4e-11, 12.4 kappa;
  the hand-written and generated 7.7e-11 are J), in 12: `4 1.0
  0.15506279215964633 0.2568972092772443 -2.909664211316309e-11
  0.022168254755127226 0.4902406748862591 0.29449095474135034 0.0 0.0
  0.15506279215964633 -1.0087679869580854e-06 0.0 0.0 0.022168254755127226
  0.0`.
- **Vieta line-edge offsets** (dbb2217): fs 10746, hand-written, 5.69e-12 /
  4.8e-13 / 5.69e-12 (f2adce4 3.8e-12; kappa 1.7e-16), in 8: `4 1.0
  297.619542010204 0.352413261575437 1.2575023301933924e-07 0.0
  0.014552824695359588 0.05516351936371848 0.015954861766150564 0.0
  0.352413261575437 0.0 0.04930018064244129 -6.689849933305886e-11
  0.00818647885741143 0.0`; fs 13316 (a physical 2 -> 2 box), hw/gen/expr,
  1.9e-12 / 4.1e-14 to 4.4e-14 / 1.9e-12 (kappa 2.1e-16), in 16: `4 1.0
  62.49314912512144 0.8401490188682588 0.0907637725398564 0.0
  86.57492489348401 -2.498445106837112 0.028691482250267415 0.0
  62.49314912512144 -2.2093188209890053 0.030251225874558578
  -8.026858086222152e-07 0.030251225874558578 -8.026858086222152e-07`; fs
  16627, hand-written, 1.60e-12 / 1.2e-13 / 1.60e-12 (f2adce4 1.41e-12;
  kappa 6.9e-14), in 24: `4 1.0 16.97289510395186 -18.130071000049583
  0.009894566586397417 3.4579339967648344 19.31861875469964
  5.023825017222337 3.457933996764513 -8.283223276663037e-12
  16.973017472593632 -8.452762982897872e-08 6.4299924328639815
  -5.073073133072303e-05 4.015267637589443 0.0`.
- **Route rounding** (regenerations): fs 4668, SymJIT, 3.92e-12 / 3.9e-13 /
  3.92e-12 (f2adce4 4.1e-12; kappa 1.2e-16), in 4: `4 1.0 0.0
  0.007445054013565638 0.0 0.017185642795225723 -0.0720362936036969
  -0.05716382232405437 50.80414335012342 0.0 0.007445054013565638
  -1.9907093696601063e-05 13.305176728220536 -1.8398052009486262e-05
  4.598584079263386 -5.6783041077624115e-09`; fs 14473, generated, 1.77e-12 /
  1.0e-13 / 1.77e-12 (f2adce4 1.72e-12; kappa 1.3e-16), in 8: `4 1.0
  -21.059190935626223 0.0 3.212826502040012e-10 -0.005388485655161519
  0.023710975762189373 -2.3244568484745787e-07 4.673081630052052
  -2.4789513702080443e-07 9.97010965981591 -1.4111835911486574e-05
  0.6967910957518927 0.0 0.0 0.0`; fs 14158, SymJIT, 1.85e-12 / 7.9e-15 /
  1.85e-12 (f2adce4 1.1e-14, so from 59b96e5's regeneration; kappa
  1.0e-16), in 24 on some route: `4 1.0 0.0 -30.54329623783176
  -0.004775170598803929 1.2811282643436634 -1.7388541862237215
  -27.47779900584406 0.13891511859525946 0.0 0.13450611158294978
  -3.1883451577633146e-06 2.0281683293083486 -0.9363637293480552
  2.3140068531658553 -0.0012961339967214824`. The R3 rows of the next
  section are the same class.
- **NaN where 09b2c74 was within 1e-6** (the chart tests of cdb4318, 97efa0f on claude/olm-final-2;
  f2adce4 had 09b2c74's values): fs 14026, generated and interpreted, NaN /
  5.2e-14 (the hand-written route keeps 4.9e-14; SymJIT NaN where 09b2c74
  was wrong by 0.98, N-bbad), in all 24: `4 1.0 -5.085267836351783e-13
  53.9782906826277 -172.99009918320144 -1.0049754268051605e-12
  -0.0855449903579419 -3.454085894812672e-12 0.0 0.0 36.19673827634061 0.0
  0.02284707832426455 0.0 0.8773920678951894 -3.389275623047775e-08`; fs
  14323, hw/expr/jit, NaN / 4.0e-7, 9.5e-7, 9.5e-7 (kappa 3.75e-7; the
  generated route's 1.35e-6 at 09b2c74 is N-bbad), in all 24: `4 1.0
  0.16992060055218752 -0.04652908500956608 7.465665692184823
  0.16992060055218752 0.012400486648567587 129.10366892402206
  0.16992060055218752 0.0 0.16992060055218752 0.0 7.465665692184823
  -4.028993948238839e-10 0.0 0.0`. In fs 14323 the other chart's K keeps
  4 digits (not a residue): taking the selected chart there as well (any
  K below 2^-26 of its terms) gives 09b2c74's values but made SymJIT wrong
  by 0.13 on v1 11845. A NaN is detected and rescued by NinjaRS.
- **`Float` at 128 bits** (F-R1 kind; 256 bits right): fs 14026, native,
  2.45e-8 at 33 tracked bits / 6.6e-11 at 42 / the same (interpreted 5.8e-9
  against 8.8e-10, not ten times), in 24; fs 15340, native, 6.0e-11 at 32
  tracked bits / 5.2e-12 / the same (kappa 1.2e-16), in 24: `4
  464.47260840716086 -1.3219919554631232e-14 -0.0033320115769994514
  -0.005202785845388104 -265.92474425838407 -10.333806188306571
  -255.59405979100657 0.017130019894277274 0.0 0.0 0.0 0.027353075628164823
  -7.463400754425057e-08 5.146807288512516 -0.03522151085733567`.
- **Within 10 kappa_1ulp (J and J*)**: fs 289 hw/gen/expr 6.76e-10 /
  6.3e-11 (kappa 5.0e-10); fs 11533 hand-written and generated 7.7e-11 /
  6.1e-12 (1.4e-11); fs 15780 SymJIT 4.86e-11 / 3.2e-12 (8.0e-11); fs 16941
  interpreted 5.87e-6 / 5.87e-6 / NaN (1.36e-6; silent against 59b96e5); fs
  15744 generated and interpreted 1.55e-6 / 1.2e-7 (J*, 5.2e-7, 1.0e-8
  without on-shell-breaking moves); fs 14323 generated 1.35e-6 at 09b2c74
  and NaN now (N-bbad).
- **NaN where 09b2c74 was wrong by more than 1e-6 (N-bbad)**: fs 1255 (all
  routes, 09b2c74 9.0e2 to 6.0e5), 14026 SymJIT (0.98), 14323 generated
  (1.35e-6), 14923 hw and SymJIT (1.75, 1.2e-2), 15753 (all four routes,
  0.55 to 2.7), 15754 (hw/gen/SymJIT, 11 to 13), 16134 (hw/gen/expr, 9.8e3
  to 7.5e4), 16941 hw/gen/SymJIT (4.7e-6 to 6.0e-6), and fs
  15189 in the 12 relabellings where 09b2c74 was wrong.

On the relabellings of the 30 finding boxes (rl), the same boxes: fs 1531
and 3566 in 12, 4668 in 4, 6653 in 24, 10746 in 8, 11533 in 12 (J in 24),
13316 in 16, 14158 in 24, 14473 in 8, 16627 in 24 (finite, at most the
errors above), 14026 and 14323 NaN in 24.

**Soft-box limitations found by the fresh sweep** (the same at 09b2c74 or
NaN where it was wrong): fs 1798 (a soft box whose two on-shell partners
are joined by a zero leg, `p2 = 0` between `m2 = 39.37` and `m3 = 28.37`,
with a complex third mass: `4 1.0 56.77097630775773 0.0 28.37080267809542
0.00379870900697696 0.012896451236664963 39.36912403067149
3.425293625843942 -3.3786226566159404e-08 39.36912403067149 0.0
28.37080267809542 0.0 0.0 0.0`) is NaN in all 24 relabellings on every
route, DoubleFloat and `Float` 512 included (09b2c74: NaN in 12, wrong by
about 150 in 12; COLLIER's COLI gives 1.8500976744e-3 in all 24, its DD
branch fails): `labelling::soft_order` refuses the exchange because the
target labelling would have a zero channel. fs 7225 (`(0, m, m', M)` with
`p1 = m2`, `s12 = p2 = m3`, `p4 = 0`: `4 1.0 0.016797359183989644
0.09530536762169647 133.98200026057233 0.0 0.09530536762169647
48.26591678889862 0.0 0.0 0.016797359183989644 0.0 0.09530536762169647 0.0
236.76399409806788 -0.0001185860571787181`) is wrong by 20 to 7e3 on every
route of both builds, 256 bits included, in all 24 relabellings. fs 68, 753,
9027 and 14233 (soft rows) are wrong in binary64 by 1.3e-4 to 1.0 at
09b2c74 and now, identically; DoubleFloat is right.

## The final candidate (2026-10-20)

`claude/olm-final-2`, on 09b2c74: the combined candidate of the section
below (f2adce4) and the following changes, with one regeneration
(`generated_d0.rs` 2,474,110 bytes, the D0 asset 1,459,142 bytes; the
definitions take 93,730,882 bytes against the gate's 100,663,296, f2adce4
100,370,779). References are COLLIER 1.2.9 in quadruple precision (COLI and
DD branches, exact decimal inputs), "worse" means more than ten times worse
and above 1e-12, and the routes are hand-written, generated, interpreted,
SymJIT and DoubleFloat (hw/gen/expr/jit/dd).

- **PR #9 and PR #10** of alphal00p/oneloopmaster are folded in with their
  tests: exactly real contour squared norms and real sheet signs at mixed
  `Float` precision (`tests/box_limits.rs`, `tests/mixed_precision_sheet.rs`,
  `box_complex::tests::composite_norms_keep_real_comparisons`), and the
  reclassification of exact infrared box limits after on-shell chart
  exchanges (`tests/d0_one_mass_chart.rs`). PR #9's norm keeps its body as
  an inlined registered function, so that each use holds one copy of its
  operand; written out, it moved 15 ill-conditioned rows of the verifier's
  sweep on the generated and interpreted routes by up to 1.1.
- **Contour charts with the Gram determinant at the rounding level** (the
  regression of f2adce4's zero-external reflection, below): a `K` that rounds
  to zero puts its far pole 2^52 beyond the near one; where the selected
  contour is flat (`2^56 |a1|² < |g|² + |beta h|²`) the other contour's chart
  is evaluated too and taken unless the two agree to 2^-26; the reflection is
  skipped where a reflected conic has a double contour direction; and in
  binary64, where no chart is regular, a chart that does not keep 10 bits
  under a one-unit move of its line (or 13 bits with its `K` at the rounding
  level) returns NaN, which NinjaRS detects and rescues
  ([PRECISION.md](PRECISION.md)). The thresholds are fixed powers of two, so
  that hw, gen, expr, jit and dd take the same chart on binary64 inputs;
  `Float` takes the selected chart. Both statements hold on the sets of this
  section; the second verifier's fresh sweep found boxes where the routes
  decide differently (fs 15189: 2.9e-2 hand-written, NaN on the other
  routes) and boxes whose chart is finite and wrong on every route (fs 16695
  and 15926), fixed in the section above.
- **One- and two-mass sheets** (Denner and Dittmaier's, with a width): the
  one-mass mass-pair complement from the roots' offsets, a small two-mass
  complement from the product in binary64, and in both the root nearer the
  branch cut taken across it, in registered functions
  (`tests/d0_one_mass_width.rs`, `tests/d0_two_mass_width.rs` with
  `tests/data/two_mass_conditioning.txt`).
- **Contour line edges**: no line product at an exactly double real pole
  (four relabellings of a box of the streams' sweep returned 7.6e18 times the
  value); in binary64 a line product that disagrees with its pole's rounded
  offset is probed, and the products are kept only where the poles' offsets
  scale with `a1` and not for a pole within 8 offsets of an end of its edge
  (`tests/data/contour_pole_errors.txt`,
  `tests/data/contour_line_products.txt`; since 2026-10-21 a product that
  `a1` does not account for is probed instead of dropping the chart's
  products, above).
- **Zero-external reflection** only for a line nearer the corner by more
  than 1/128 and a pole that one of the original chart's `tfun` selects with
  an offset from the corner's vertex that the rounding resolves poorly
  (`zero_external_reflects_only_*` in `tests/d0_contour_regressions.rs`).
- **Soft boxes in every labelling.** The D0 expressions, like OneLOop's d0,
  find soft limits on the four external invariants only. Of the 24
  labellings of a soft-divergent box with a massless line whose leg partners
  are on their mass shells, the 16 that put such an on-shell pair on a
  channel (`s12` or `s23`) took the finite formulas: no pole in binary64, a
  wrong value or NaN in DoubleFloat and `Float` (pyAmpliCol's
  e+ e- -> t t~ h box with lines (m, m, 0, m), m = 173). The numerical D0
  entry points now evaluate such a point in the labelling of the same box,
  of the two that exchange the channel pair, whose channels are nonzero and
  hold no such pair (`labelling::soft_order`, with the expressions' own
  exact zero tests); every other point keeps its labelling. Covered:
  `NativeEvaluator` in every carrier (binary64 after its hand-written
  sectors, whose three-mass sector leaves these points to the generated
  code), `ScalarEvaluator` (SymJIT), the Expression backend,
  `PrecisionEvaluator`'s interpreted `Float` and the public masters' numeric
  callbacks. Symbolic uses do not pass an entry point:
  `OneLoopExpressions::d0`, `crate::d0`, `d0_mapped` and evaluators built
  from them, the transparent master definitions inlined through the function
  map, and inspection. **Callers that evaluate the D0 expressions themselves
  must relabel the point with `oneloop::regular_d0_labelling` first**
  (NinjaRS's arbitrary-precision `SymbolicaMasters` does). Six such boxes in
  all 24 labellings on every route (`tests/d0_infrared_relabelling.rs`,
  OneLOop 3.6 quadruple-precision references; COLLIER agrees in all 24) now
  pass everywhere, where 09b2c74 failed 16 (or 8) labellings; on the
  verifier's sweep 32 rows are better on every binary64 route and at `Float`
  128 and 256 bits, on the streams' sweep 344, and none is worse anywhere.
- **Test model** of `d0_two_mass_width`'s 256-bit shell check: it asserts the
  input-cancellation floor (at least `bits - 56` tracked and accurate bits)
  instead of accuracy within 8 bits of the tracked precision, which
  numerica's estimate does not give: against a 1024-bit reference (native
  and interpreted agree to at least 970 bits), 09b2c74 falls 8.23, 9.09,
  9.63 and 7.75 bits short at 128, 192, 256 and 320 bits on row 6 of
  `tests/data/two_mass_shell.txt` (native; 253 tracked and 243.37 accurate
  bits at 256), and at most 5.85 bits on row 0 (the rows' 90-digit
  references cannot measure more than about 299 accurate bits, so the
  16.1 bits at 320 bits quoted earlier were an artifact of the reference).
  The minimum accuracy over the 9 rows, 24 relabellings and both routes is
  unchanged (202.53 bits at 256 against 202.33 at 09b2c74); single
  relabellings lose up to 40.3 bits against 09b2c74 (row 0; 356 of the 1,728
  entries at 128 to 320 bits more than 10 bits), all above the floor and the
  same since f2adce4 ("`Float` keeps the precision it tracks" below).
  `tests/float_exact_products.rs` makes the same 8-bit claim; on its rows the
  claim holds wherever it binds (at most 4.2 bits short), and it is
  unchanged.

### Results

Rows worse / better / silent / new NaN against 09b2c74 (silent: finite and
wrong by more than 1e-6 where 09b2c74 was non-finite or right). v1 is the
verifier's fresh sweep (32,600 rows, seed 777013: every mass pattern, widths
1e-14 to 0.7, corners and vanishing Gram determinants, 400 corner boxes in
all 24 relabellings, 3,000 C0), CT the contour stream's 3,387 scan rows, MS
the 843 distinct NinjaRS master calls, NR the 7,644 calls of NinjaRS's rank
n+2 suite with their mu², g1 12,000 generic rows (seed 424242), TM the
128,462 two-mass rows, SW the 37,642 rows of the oracle sweeps (seeds 1 to
4, all 24 relabellings), TD the tests' rows in their 24 relabellings:

| set | hw | gen | expr | jit | dd |
|---|---|---|---|---|---|
| v1 | 5/566/1/5 | 3/586/2/3 | 5/573/1/4 | 11/597/2/6 | 0/105/0/0 |
| CT | 9/673/0/9 | 11/678/0/9 | 11/677/0/9 | 13/676/0/9 | 1/42/0/1 |
| MS | 0/1/0/0 | 0/1/0/0 | 0/3/0/0 | 0/1/0/0 | 0/0/0/0 |
| NR | 0/0/0/0 | 0/0/0/0 | 0/2/0/0 | 0/1/0/0 | 0/0/0/0 |
| g1 | 0/133/0/0 | 1/142/0/0 | 1/150/0/0 | 0/148/0/0 | 0/61/0/0 |
| TM | 2/2122/0/0 | 2/2122/0/0 | 1/2121/0/0 | 1/2106/0/0 | 0/1844/0/0 |
| SW | 8/1420/0/8 | 8/1436/0/8 | 8/1400/0/8 | 8/1400/0/8 | 0/1300/0/0 |
| TD | 0/1096/0/0 | 0/1120/0/0 | 0/1088/0/0 | 0/1084/0/0 | 0/12/0/0 |

(The TM better counts follow the first verifier's reference rule; with the
second verifier's, COLI where its DD branch agrees to 1e-13 and else where
the 256-bit values agree, they are 2133/2133/2132/2119/1844. New NaN counts
only referenced rows here.) Every one of these 114 worse or silent route
entries is an accepted exception listed below. Against f2adce4: v1 27/33/31/41/0 rows worse, of
which 15/26/26/24 are new NaN (the chart tests return to NaN boxes that
f2adce4 got right but that were NaN or wrong at 09b2c74; the others are
f2adce4's gains given back to 09b2c74's level); CT 0/1/1/0/0; SW 4 on the
interpreted route (the double real pole, now NaN as at 09b2c74); TD 48 per
route (`tests/data/contour_unresolved.txt`, NaN by design); MS, NR, g1 and
TM none. Finite and wrong by more than 1e-6 where f2adce4 was right: v1 395
(generated, SymJIT; J* below) and v1 14736 on the generated route (0.65, the
value of 09b2c74; its hand-written route is NaN, R8 below).

- **Vanishing Gram determinants** (the verifier's finding on f2adce4). Box A
  (`4 1.0 7.83411141296995 0.0 0.05539849041093991 1.2025720308924596e-10
  0.44413758479240445 0.9771700504003817 0.18581901508118223 0.0
  0.05539865304260027 -6.430669304543917e-07 0.05539865304260027
  -6.430669304543917e-07 0.0 0.0`, -25.40966982844103 + 0.884429096461921i):
  f2adce4 O(30) wrong on every binary64 route in all 24 relabellings, now
  within 3.3e-14/3.6e-14/5.1e-14/6.4e-14 (hw/gen/expr/jit) in all 24
  (09b2c74: 12 within 1.39e-13, 12 non-finite). Box B (`4 1.0 0.0
  3.598382043529174 24.453495802712766 24.453495802712766
  -1.255476097339391e-09 502241.0517849417 0.6808521119960873
  -0.02266926845550123 0.02844863278465547 0.0 0.0 0.0 24.453495802712766
  0.0`): f2adce4 4.7e-3 in all 24, now 2.2e-9 (SymJIT 4.9e-9) in all 24
  (09b2c74: 2.2e-9 in 12, non-finite in 12). Both are rows 1 and 2 of
  `tests/data/contour_gram_rounding.txt`. Of the 26 boxes of v1 that 09b2c74
  returned as NaN and f2adce4 as finite values wrong by up to 2.8e2, 11 are
  NaN again, 14 right and one (v1 5151, 4.1e-6) within 10 kappa_1ulp
  (below); 7 boxes that 09b2c74 already got wrong by more than 1e-6 are NaN,
  none that it got right; 23 boxes that f2adce4 got right (NaN or wrong at
  09b2c74) are NaN again. On the 1,649 referenced v1 rows with
  `|Gram|/s³ < 1e-15` (hw/gen/expr/jit): beyond 1e-6 at 09b2c74
  55/56/57/57, f2adce4 61/59/62/64, now 40/40/39/41, every one of them also
  beyond 1e-6 at 09b2c74; non-finite 141/142/138/141, 11/11/11/9, now
  16/18/19/17. The largest errors there are 4.3 (hw), 1.6 (gen) and 9.0
  (expr, jit), on rows 09b2c74 gets as wrong (5.5e2 at 09b2c74 and
  f2adce4).
- **The tests' rows, worst of the 24 relabellings** (hw/gen/expr/jit): the
  vanishing-Gram box of `tests/d0_vanishing_gram.rs` 2.6e-14/1.9e-14/
  3.4e-14/3.0e-14 (09b2c74 3.6e-11); the seven boxes of
  `tests/data/contour_regressions.txt` at most 7.7e-12 (09b2c74 up to 1.3e2);
  the nine of `tests/data/contour_gram_rounding.txt` at most 1.6e-11 but row
  2 (box B); `tests/data/two_mass_conditioning.txt` at most 5.7e-14; 199 of
  the 200 boxes of `tests/data/near_threshold_boxes.txt` within 3.92e-12 and
  line 73 at 3.2e-10/5.6e-11/1.35e-11/2.4e-10.
- **NinjaRS master calls.** Of the 309 distinct D0 calls (MS), 99 (hw), 191
  (gen), 191 (expr), 219 (jit) and 200 (dd) change their bits against
  09b2c74 (f2adce4: 107, 206, 217, 215, 221); the worst errors against
  COLLIER are 1.43e-12, 1.42e-12, 1.29e-12, 5.17e-12 and 4.6e-15 (f2adce4
  1.85e-12, 1.72e-12, 1.83e-12, 5.17e-12; 09b2c74 3.60e-11, 3.59e-11,
  8.89e-12, 3.41e-11). The 534 C0 calls keep their binary64 and DoubleFloat
  bits (worst 8.3e-14, 9.1e-14, 6.9e-14, 6.2e-14). On NR, D0 (3,510 calls)
  worst 7.53e-11, 7.48e-11, 7.51e-11, 2.58e-11 (f2adce4 the same, SymJIT
  2.43e-11; 09b2c74 SymJIT 1.66e-11) and C0 (4,134) unchanged. No call of MS
  or NR is worse than f2adce4 or 09b2c74 on any route. At `Float` 256 bits
  (native/interpreted; errors against the 512-bit value): MS D0 moves on
  308/303 of 309 calls against 09b2c74, by at most 4.5e-71/1.3e-70, and on
  157/160 against f2adce4, by at most 5.3e-72/8.4e-72; the C0 calls equal
  f2adce4's; NR D0 moves on 3,331/3,361 of 3,510 against 09b2c74 (at most
  2.0e-66/1.7e-67, where 09b2c74 was that far off) and 2,552/2,603 against
  f2adce4 (at most 5.3e-69/1.1e-69). No MS or NR call is more than ten times
  worse than f2adce4 at 256 bits (above 1e-70).
- **`Float`, rows worse than 09b2c74 against COLLIER** (native/interpreted):
  v1 3/27 at 128 bits and 0/0 at 256; CT 9/9 at 128 (eight relabellings of an
  on-shell scan point within 10 kappa_1ulp through an on-shell-breaking move,
  and CT 337) and 1/1 at 256 (CT 337); MS none; g1 none at 256 bits; TM
  none at 128 and 256; SW 4/4 at 128 (four relabellings of a box with a mass
  one ulp from an invariant, within 10 kappa_1ulp) and none at 256, where
  1,184 (native) and 1,180 (interpreted) values that 09b2c74 returned
  non-finite are now right (the soft-box relabelling). None is worse than
  f2adce4, and none is finite and wrong by more than 1e-6 where 09b2c74 was
  right or non-finite. The rows are listed below.
- **Still wrong, as at 09b2c74** (not regressions; the largest binary64
  errors of the sets): v1 6671 (`0CCC` collinear: hand-written 4.3,
  interpreted and SymJIT 9.0, generated right), the C0 v1 30180 (SymJIT 40,
  the other routes NaN), v1 14736 (generated 0.65, SymJIT 0.44) and 9649
  (generated 1.6); the streams' sweep box `4 1.0 0.08938220565009944
  -24.09653762111702 -0.8946129797968647 0.0 0.1941275888116677
  0.08938220565009944 0.6383507753193867 0 0.25 0 0 0 0.6383507753193867 0`
  (SW 22062 and its relabellings, against COLLIER where its COLI and DD
  branches agree, in 16 of 24: every binary64 route wrong by 2.8 to 4.4 times
  the value in 4, SymJIT by 2.5e2 in 4 more where the other routes are NaN;
  DoubleFloat within 1.2e-11), SW 24778 (generated 7.2e6, SymJIT 2.0e7;
  09b2c74 1.2e8 and 2.0e7) and SW 23555 (`4 1.0 0.04370275218626788
  8.017760183204489 -3.3132567693387864 0.04370275218626788
  0.018500872985981687 0 0.04370275218626788 -0.0039588657598020865 0 0 0 0 0
  0`: every route, DoubleFloat and `Float` included, returns a finite part of
  exactly 0, where COLLIER's two branches agree on 496 + 370i).

### Gates

- **Full suite** (`cargo test --release --locked --features
  generated-evaluators --no-fail-fast -- --test-threads=1`, with the
  checkout's local root `Cargo.lock`, untracked as at 09b2c74): 264 passed, 0
  failed, 9 ignored in 62 test binaries (f2adce4: 242 passed and one
  failure, the shell check whose model is corrected above).
- **Oracle sweeps** (`tests/oracle_sweep.py`, seeds 1 to 4, C0 and D0, native,
  `--generated` and `--expression`, the high-precision line at 128, 256 and
  512 bits): the 36 failure lists are 09b2c74's, but for two documented boxes.
  On seed 2, in all nine configurations, 8 relabellings of the complex
  one-mass box `4 1.0 -0.7698139610464303 -1.0396864914415214
  -29.998775868631494 81.22865927661125 8.756308210934074 0.25 0 0 0.25
  -0.018035862285877242 0 0 0 0` fail as `analytic`: OneLOop, the sweep's
  oracle, is wrong there (one-mass boxes with a complex mass, above) and
  COLLIER agrees with the candidate to 6e-16. On seed 3 at 128 bits, 16
  relabellings of `4 1.0 0.25 0.25 1.0 4.653123451185273 2.2370332401120216
  2.25 2.2370332401120216 -0.0013515238817967694 0.25 -0.0010321989957060157 0
  0 1.0 0` move from `analytic` to `binary64` on the three routes: their
  binary64 values keep 8.4e-9 against OneLOop, as at 09b2c74, and their
  128-bit values, off by 1.6e-3 at 09b2c74, are now right. Every other list is
  identical.
- **Instructions per binary64 D0 call** (perf stat instructions:u,
  difference of 2001 and 1 repetitions) on the 22 NinjaRS master streams
  against 09b2c74: geometric mean -0.305%, at most +0.443% (syn-n6-r7;
  syn-n6-r5 +0.336%, syn-n6-r6 +0.329%, the other 19 streams -0.007% to
  -0.770%), within the 0.5% gate on every stream.
- **Definitions** 93,730,882 bytes (gate 100,663,296; headroom 6,932,414);
  the complete D0 expansion fits `tests/full_expansion.rs`'s 10M budget.
- **Regeneration**: a second `generate_native` and `rebuild_evaluators D0`
  from the committed sources gives identical SHA-256 for every asset and
  generated file; `cargo metadata --locked` passes in the three build roots.
  The C0 SymJIT asset is not rebuilt (as in PR #9's own regeneration): from
  the current expressions `rebuild_evaluators C0` gives 85,547 bytes
  (`7a7bd442...`) instead of the committed 102,902 (`f8c7109a...`); the
  binary64 C0 values of the NinjaRS master calls (MS, NR) are bitwise those of
  09b2c74 on every route, SymJIT included. (Rebuilt on 2026-10-21, above.) `src/sheet_exact.rs` still reads
  the debugging switches `ZZ_GEN_LOG` and `ZZ_GEN_DD` (from 3cbe08d).

### Accepted regressions against 09b2c74 (decision of 2026-10-20)

The user accepted the following route entries as documented exceptions. No
other binary64 route entry is more than ten times worse than 09b2c74 and
above 1e-12, and no other finite value is wrong by more than 1e-6 where
09b2c74 was right or non-finite, on the sets above (the first verifier's
sets, the streams' sweep and the tests' rows). The second verifier's fresh
sweep found 45 more binary64 route entries on 17 of its 36,840 rows, 3 of
them silent; after the fixes of the section above 25 entries on 12 rows
remain, none silent, listed there. kappa_1ulp is the
largest relative change of D0 when one of the ten real binary64 inputs (the
six invariants and the real and imaginary parts of the four masses; not
mu²) moves by one ulp in either direction, evaluated at 256 bits where
09b2c74 and the candidate agree: a row within 10 kappa_1ulp is justified
(class J), as no algorithm that rounds its inputs once can promise better;
J* rows are within 10 kappa_1ulp only through a move that breaks an exact
input equality (an on-shell or equal-mass relation).

- **Route rounding of the contour chart's arguments (R3)**, the same at
  f2adce4 (errors now / 09b2c74 / f2adce4):
  - v1 19098, SymJIT, 2.89e-12 / 1.38e-13 / 2.89e-12 (kappa 1.6e-15):
    `4 1.0 0.028456754273287505 28.53500714956581 11.340618430715825
    -0.2108511787556433 96.16211930117166 -59.84315524140586
    19.930972340236043 -6.579379864209245e-05 19.930972340236043
    -6.579379864209245e-05 28.53500714956581 -0.42841797262874637
    3.8975973963186417 -8.175209021663547e-12`.
  - v1 19320, interpreted, 1.18e-11 / 1.06e-12 / 1.18e-11 (kappa 4.4e-15):
    `4 1.0 -0.011597221745223632 25.285481364080724 6.585614318476484
    0.030754490595793042 242.05052338832093 -210.8487715706494
    25.285481364080724 -3.4401499278722525e-11 25.285481364080724
    -3.4401499278722525e-11 0.0 0.0 0.030754490595793042 0.0`. The line
    products also lose on the hand-written route here (5.2e-13 against
    4.4e-14 with `y - z`), but no product rule separates this box from the
    vanishing-Gram boxes that need them.
  - CT 2168, SymJIT, 4.54e-12 / 2.55e-13 / 5.32e-12 (kappa 3.8e-16):
    `4 0.01 0.0027763524551391383 0.0013364056189293424 0.05120181816720371
    0.0008105366415512052 0.03220975168078065 0.006587112684859578
    0.027065858857493085 0.0 0.01994971052861607 0.0 0.0035726174395936216
    0.0 0.030735674662001472 0.0`.
  Each loss is reproduced by the route's own rounding of the chart's
  arguments (`ordered` on the hand-written route's chart inputs); a fix
  needs those arguments evaluated identically on every route. The
  candidate's worst binary64 route over 09b2c74's worst route is 4.50
  (v1 19098), 11.05 (v1 19320) and 0.99 (CT 2168); on the R6 rows below
  3.25, 1.00, 1.00 and 1.52 (CT 38, 43, 48, 55) and 1.25 (g1 9239).
- **The CT 31-70 scan around a vanishing Gram determinant (R6)**, the same at
  f2adce4: `4 1.0 6.329201971579061 102.98802777639872 3.3205818912242475
  3.320581891224248 13.28232756492699 51.33798400005867 3.343448813357665 0.0
  3.343448813357665 0.0 3.343448813357665 0.0 0.0 0.0` (CT 38) and its
  neighbours in the last digits: CT 38 generated 5.56e-12 / 6.15e-14 /
  5.48e-12 and interpreted 5.83e-12 / 2.71e-13 / 5.81e-12; CT 43 generated
  1.26e-11 / 1.15e-12 / 1.26e-11; CT 48 SymJIT 2.58e-11 / 3.02e-13 /
  2.58e-11; CT 55 SymJIT 1.97e-11 / 1.03e-12 / 1.96e-11 (kappa 5.549e-15).
  On the 40 rows of this scan the geometric mean of the error per binary64
  route is 1.1e-11 to 1.3e-11 at 09b2c74 and 9.7e-12 to 1.3e-11 at f2adce4
  and now (also without the `r1` complements), single rows 6e-14 to 6.5e-11:
  09b2c74 had one lucky route on these rows. The error comes from the
  rounding of the three `tfun`'s arguments in `general`, not from `alpha`.
- **g1 9239 (R6)**, generated 2.58e-11 / 1.26e-12 / 2.58e-11 and interpreted
  2.57e-11 / 1.48e-12 / 2.57e-11 (kappa 1.0e-16): `4 1.0 0.0
  -121.2552217479067 -138.88837706189003 -0.2362666283315987
  0.1454848278533424 -0.04327231551251378 0.0 0.0 0.6217945335910897 0.0
  4.0386278704465886 0.0 0.1414204746299639 -9.099246028823919e-07`. In one
  `r1` call, the complements `1 - x` of 14dd171 break an error cancellation
  between the logarithms and the dilogarithm (2.4e-8 with them, 3.5e-12
  without); taking them only where the hand-written kernel does fixes this
  row but makes v1 6738 (to 9.7e-10), 13553 (to 5.7e-8), NR 5532 (to
  1.2e-11) and 8 g1 boxes worse.
- **NaN beats wrong.** NinjaRS detects a non-finite master and rescues it in
  higher precision, but cannot detect a silently wrong one. So a NaN where
  09b2c74 returned a finite value already wrong by more than 1e-6 is
  accepted:
  - N-bbad, 9 entries, where 09b2c74 was wrong by 9.8e-6 to 0.77 (now NaN):
    v1 2057 SymJIT (09b2c74 9.3e-4), v1 3794 hw and SymJIT (1.2e-4,
    3.4e-4), v1 9649 interpreted and SymJIT (0.77, 9.8e-6), v1 18033 on
    all four routes (5.6e-4, 1.4e-3, 6.7e-4, 6.6e-4); the chart tests cannot
    resolve these charts in binary64.
  - R8, 78 binary64 entries where 09b2c74 was wrong by 0.25 to 3.5e2: v1
    10082 (all four routes, 3.5e2), 14736 (hw, 1.1) and 17188 (all four,
    0.92 to 5.8); CT 337 (all five routes; open limitations below) and eight
    relabellings each in CT (rows 344 to 363) and SW (rows 4062 to 4082) of
    the conic singular at a massless corner (0.25 on every binary64 route;
    open limitations below); and at
    `Float` 128 and 256 bits CT 337 (4 entries; 09b2c74 wrong by 81 to
    3.4e2).
- **Within 10 kappa_1ulp (J)**, now / 09b2c74 / kappa_1ulp: v1 333 SymJIT
  1.79e-12 / 1.01e-13 / 7.4e-13; v1 5151 on all four routes 4.14e-6 /
  non-finite / 8.4e-3 (the only finite value wrong by more than 1e-6 where
  09b2c74 was non-finite); v1 16427 SymJIT 2.14e-11 / 1.28e-13 / 7.9e-11;
  v1 18577 SymJIT 1.14e-11 / 2.23e-13 / 1.7e-10; TM 65774 hw and gen
  2.66e-12 / 1.58e-13 / 2.4e-12; TM 103232 hw, gen and expr 5.85e-12 /
  2.92e-13 / 7.7e-12; TM 103983 SymJIT 1.42e-12 / 1.06e-13 / 1.9e-12; and at
  128 bits v1 16793 native 8.75e-11 / 2.31e-12 / 5.1e-7 (4 tracked bits) and SW 36841, 36844,
  36855 and 36858 (four relabellings of `4 1.0 0.0 11.923773121437657
  10.57262939762882 0.04061179735123128 11.923773121437657
  0.04061179735123127 0 0 0 0 11.923773121437659 0 0.04061179735123128
  -6.341526100574037e-06`, a mass one ulp above two invariants) native
  3.29e-12 and interpreted 1.17e-12 / 7.8e-17 / 18.6 (0.51 without
  on-shell-breaking moves; 10 and 9 tracked bits), the same at f2adce4; and
  at 128 bits the box `4 1.0 -8.493998272567758 40.673484401155186
  0.045772794745771575 0.6178812442530255 1.0 1.3123023694541682 1.0 0.0
  1.0 0.0 0.0 0.0 0.04577279474577158 0.0` of the seed-6 oracle sweep (`p3²`
  one ulp below `m4²`, `m3 = 0`) in 8 relabellings, native, generated and
  interpreted, 6.7e-8 and 1.8e-7 / 1.3e-14 and 7.9e-15 against OneLOop in
  quadruple precision (at least 22.4 accurate bits over the 24
  relabellings, 46.1 at 09b2c74), kappa_1ulp 0.94 (the move that puts `p3²`
  onto `m4²`; 7.8e-2 without it), the same at f2adce4, 192 and 256 bits
  unchanged. At `Float` precision the binary64 inputs are exact, so binary64
  conditioning alone does not bound these rows' errors: they are justified
  by that conditioning (the inputs carry binary64's rounding) together with
  numerica's tracked precision, which reports the loss (the tracked bits
  given).
- **Within 10 kappa_1ulp only through an on-shell-breaking move (J*)**, now /
  09b2c74 / kappa_1ulp (without such moves): v1 395 generated 5.72e-6 /
  6.04e-7 / 1.1e-6 (7.8e-8) and SymJIT 3.65e-6 / 3.18e-7 (silent: finite
  and wrong by more than 1e-6 where 09b2c74 was within 1e-6); CT 778
  interpreted 1.11e-12 / 9.58e-14 and SymJIT 1.13e-12 / 6.40e-14 / 1.8e-11
  (6.4e-15); at 128 bits the eight relabellings CT 76, 78, 79, 80, 85, 87, 95
  and 96 of an on-shell scan point, native and interpreted, 3.43e-12 /
  3.4e-15 / 1.8e-11 (6.4e-15), 43 tracked bits (about 38 accurate).
- **`Float` at 128 bits next to a vanishing Gram determinant (F-R1)**, from
  14dd171's expression changes (now / 09b2c74 / f2adce4, tracked bits now):
  v1 7952 native 2.13e-11 / 3.3e-18 / 1.03e-11 and interpreted 1.85e-11 /
  3.3e-18 / 1.85e-11 (40 bits; 09b2c74 76; fixed on 2026-10-21, 3.3e-18 at
  73 bits); v1 12595 interpreted 1.89e-10 /
  1.40e-11 / 2.06e-10 (25 to 39 bits, as at 09b2c74); v1 15453 native
  5.24e-12 / 4.3e-16 / 5.18e-12 and interpreted 2.21e-11 / 4.3e-16 /
  2.19e-11 (38 to 41 bits; 09b2c74 55 to 57). 256 bits are right.
- **The 0RRC collinear box at 128 bits from the Symbolica pin (F1)**: `4 1.0
  50.261921198182165 0.15614357216746133 82.77728497126748
  1.6038370763317959e-13 0.2547782935219607 -12.782604724181319 0.0 0.0
  16.407592918401335 0.0 82.77728497126748 0.0 15.201391682789128
  -0.00022669601519667755` in all 24 relabellings, interpreted, 8.34e-9 /
  2.75e-10 / 8.00e-9 at 30 tracked bits (the old pin e679e8c1 gives
  2.9e-10); the native route is within ten times of 09b2c74 and 256 bits
  are right.

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
  most; 3.92e-12 in the final candidate, SymJIT), where 34 were beyond it and
  22 beyond 1e-9 before (OneLOop's double
  precision: 32 and 16). On 1,000 more such boxes (seeds 2 to 6) 4 to 6 per
  route remain beyond 1e-11 (at most 3.2e-11), against 110 to 117 before; on
  fresh seeds (101 to 105) up to about 8e-11.
- **Remaining.** Line 73 of the scan keeps 1.3e-11 (interpreted) to
  3.2e-10 (hand-written), as before (OneLOop 5.5e-11; the final candidate
  1.35e-11 to 3.2e-10): a nearly double contour direction in one `tfun` and
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
  routes); now 3, at most 5.1e-3 (below; the final candidate's chart tests
  keep that box within 4.8e-13).
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
  least 80 and 208 bits there; the 256-bit check now asserts at least
  `bits - 56` accurate bits, see the open limitations). 9a898ee fails both
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
(9a898ee) fail on the same rows.

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
  exact. This is symbolica-dev/symbolica#97. The pin is now `9cdc2949`
  (`e679e8c1` plus that commit), taken from branch `ninjars-pin` of
  alphal00p/symbolica, where NinjaRS takes its Symbolica from too (also
  branch `claude/oneloopmaster-integration-float-products` of
  ValentinHirschi/symbolica). Binary64 is bitwise unchanged on every route
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
- **The two-mass threshold `r14 = 1`** is the example above. In `4 1.0
  -7.158359142869336 -0.006328911565103391 0 0 0 0.01281606414601131 0.25
  -4.478335170728946e-05 0 0 0.25 -4.478335170728946e-05 0 0` (equal complex
  masses on opposite lines, a vanishing invariant between them, two more
  vanishing invariants) the native 256-bit route gave 3.072e4 + 19i,
  9.728e3 - 2.307e7i or 1.9456e4 - 4.614e7i, numbers of a few bits, for
  8.6969 + 5.8296i in 12 of the 24 relabellings, with the two-mass sheets
  and before them (09b2c74: 2.8672e4 + 19i and the like), while binary64 on
  every route and the interpreted 256-bit expressions were right (5e-16 and
  7e-17 against COLLIER's COLI branch, whose DD branch stops with error -10,
  and the formula in mpmath). With the pin the native route is right too:
  in all 24 relabellings, at 128, 256 and 512 bits, the native and
  interpreted values agree with COLLIER's COLI branch to the binary64
  rounding of the output.
- **NinjaRS.** Take all 3- and 4-propagator subsets of every lane of the
  rank n+2 suite (cases, controls, identity, scan) and of the final review
  (cases, cases8, cases-eucl, cases-scan, dotcases). That gives 7,645
  distinct C0/D0 calls.
  - With the previous pin, all 8 `Float` lanes (native and interpreted, 128
    to 1024 bits) agree with each other to 1e-8 on every call. None of these
    calls is affected.
  - The fix moves their values by at most 2.9e-32 at 128 bits and 4.7e-71 at
    256 bits relative.
  - NinjaRS pins Symbolica separately (`alphal00p/symbolica`, `e679e8c1`
    when this was written). Its own `Float` masters change only after its pin
    moves to a commit with #97, such as `9cdc2949`.

## The combined candidate (2026-10-20)

The sections above together, with one regeneration: `generated_d0.rs` is
2,794,680 bytes and the D0 asset 2,044,124 bytes (each 9a898ee's size plus
both branches' growth; the other families are unchanged), and the definitions
take 100,370,779 bytes (gate 100,663,296). Against 09b2c74 and against the
heads of the three branches (contour-box regressions, two-mass complements,
Symbolica pin):

- **No interaction.** On the D0 and A0 to C0 sweep rows (seeds 1 to 4), the
  NinjaRS calls (843 distinct master calls and 7,644 calls of the rank n+2
  suite and its review), the contour branch's 3,387 scan rows and the
  two-mass branch's 128,462 rows, every binary64 route (hand-written,
  generated, interpreted, SymJIT) and DoubleFloat give the bits of the
  contour branch wherever the two-mass branch leaves a row unchanged and
  those of the two-mass branch wherever the contour branch does; no row is
  changed by both, and the pin changes no binary64 or DoubleFloat bit.
- **Real 0-, 1- and 2-mass boxes** (18,270 sweep rows and 362 NinjaRS
  calls) are bitwise those of 09b2c74 on every binary64 route and in
  DoubleFloat. Their 256-bit values move only with the pin: 784 (native) and
  816 (interpreted) by at most 1e-30 relative, 16 towards OneLOop's
  quadruple precision, and 24 OneLOop-singular boxes stay non-finite.
- **`Float` routes.** The counts of the pin's section hold on this build:
  of the 25,645 stable D0 sweep rows, 12, 0, 0 and 5 (native) and 12, 0, 4
  and 4 (interpreted) at 128, 256, 512 and 1024 bits; of the 17,198 stable
  A0 to C0 rows, 23, 7, 6 and 12, and 20, 12, 8 and 5.
- **Oracle sweeps** (`tests/oracle_sweep.py`, D0 and A0 to C0, seeds 1 to 4,
  native, generated and expression routes, the high-precision line at 128,
  256 and 512 bits): all 72 failure lists are those of 09b2c74 and of the
  three branches. Seed 2's largest stable binary64 error is 1.09e-9, as at
  09b2c74 (9a898ee: 9.7e-10 native, 8.3e-10 generated and interpreted), on
  a row of that list. Row by row, against OneLOop's quadruple precision or
  the 256-bit value, 112 to 132 D0 rows per binary64 route are more
  accurate than at 09b2c74 by more than ten times, and 8 less: the
  relabellings of the conic singular at a massless corner (below). This
  holds on the oracle-sweep rows only. On the verifier's fresh sweep (32,600
  rows, seed 777013, COLLIER quad references) this candidate had 27 to 33
  rows per binary64 route more than ten times worse than 09b2c74 and above
  1e-12 (47 distinct rows, 14 of them beyond 1e-6), against 519 to 549
  better, and 26 boxes at vanishing Gram determinants where 09b2c74's NaN
  became a finite value wrong by up to 2.8e2; the final candidate (above)
  repairs them.
- **Retired instructions per binary64 D0 call** (perf stat, difference of
  2001 and 1 repetitions) on the 22 NinjaRS master streams against 09b2c74:
  geometric mean -0.12%, at most +0.29% (syn-n4-r5), at least -0.63%. On
  the samples of the two-mass sections and of the two-mass verifier (200 to
  4,000 boxes each): complex two-mass boxes +2.5% to +9.8%, boxes with a
  real mass at or next to its shell +16% to +31%, complex one-mass boxes
  -14% to -15%, real two-mass boxes +0.4%.
- **NinjaRS master calls.** Of the 309 distinct D0 calls, 107 (hand-written),
  206 (generated), 217 (interpreted) and 215 (SymJIT) change their binary64
  bits against 09b2c74, 221 their DoubleFloat bits (by at most 1.4e-31) and
  almost all their 256-bit values (309 native and 303 interpreted, by at most
  3.9e-71 and 1.0e-70, in the verifier's full-precision comparison);
  the worst errors against OneLOop's quadruple precision (COLLIER's
  quadruple precision agrees with it to 4.6e-15) are 1.9e-12, 1.7e-12,
  1.8e-12 and 5.2e-12 (09b2c74: 3.6e-11, 3.6e-11, 8.9e-12, 3.4e-11). The 534 distinct C0 calls
  keep their binary64 and DoubleFloat bits; 510 (native) and 496
  (interpreted) of their 256-bit values move with the pin, by at most 4.5e-72
  and 6.1e-73.
- **Tests.** `cargo test --release --locked --features generated-evaluators
  --no-fail-fast -- --test-threads=1`: 242 passed, 9 ignored and one
  failure, the 256-bit check of the near-shell two-mass boxes (`Float`
  keeps the precision it tracks, below), whose model the final candidate
  corrects.

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
  the known-failure fixture of `tests/d0_contour_regressions.rs`. **Fixed in
  the final candidate:** the box is row 3 of
  `tests/data/contour_gram_rounding.txt`, within 4.8e-13 (hand-written),
  8.7e-14 (generated), 4.1e-13 (interpreted) and 4.3e-14 (SymJIT) in all 24
  relabellings.
- **A conic singular at a massless corner.** `4 1.0 0.1113435943976259
  108.23851883958666 0.11134359439762591 -0.012515862621692797 9.0 0.0 0 0
  0.11134359439762591 0 9.0 0 0.11134359439762591 0` (`p1²` one ulp below
  `m2² = m4² = p3²`, massless `m1`): in the contour box's labellings the
  conic's value and gradient vanish at the massless corner to within
  rounding. Every binary64 route is wrong in every relabelling, as at
  09b2c74 (then NaN in 8, by 1.4e-2 in 4 and by 0.3 in 12; now NaN in 16 and
  1.4e-2 in 8, so 8 relabellings of the seed-1 sweep went from 0.3 to NaN;
  the same in the final candidate). Known-failure fixture of
  `tests/d0_contour_regressions.rs`. A second box of this kind, CT row 337 of
  the contour stream's scans, `4 1.0 0.01545852083001525 2.25
  -0.32983766127195113 0.0 0.01545852083001525 0.01545852083001525
  0.01545852083001525 -0.00013458303781508371 0 0 2.25 0 0.01545852083001525
  -0.00013458303781508371` (several on-shell coincidences; -26.2653 +
  15893.75i), is NaN
  in 16 of its 24 relabellings, its own included, on every route including
  DoubleFloat and `Float` at 128, 256 and 512 bits, and right in the other 8
  (binary64 within 8.1e-15, DoubleFloat and `Float` exact), the same since
  f2adce4. 09b2c74 was NaN in 8, right in 8 and wrong in 8 (by 3.4 to 46 in
  binary64, SymJIT 3.4 in 4 of them, 31 in DoubleFloat and up to 4e2 at 256
  bits).
- **128 bits next to a corner.** At 128 bits the boxes of
  `tests/d0_contour_regressions.rs` with a pole within `t²` of a corner keep
  `2^-128/t²` of its offset: the relabellings agree only to 1e-16 to 1e-29,
  and the one-ulp box is within 1.2e-11 (native) to 3.8e-11 (interpreted) in
  16 relabellings (20 at 09b2c74). Further from on shell, Float's tracked
  precision of the poles' small components costs the labelling with the
  corner at the origin up to 2.5e-7 at 128 bits (2.4e-7 at 09b2c74; `p3²`
  3e-16 below `s12`), which the reflection now also takes for 8 of the
  relabellings that kept 1e-15 at 128 bits at 09b2c74 (and gave O(1) errors
  in binary64). 256 bits keep 1e-15 (and agree to 1e-30). In the final
  candidate, boxes next to a vanishing Gram determinant lose more at 128 bits
  than at 09b2c74 (14dd171's expressions): v1 7952, 12595 and 15453 of the
  verifier's sweep keep 5.2e-12 to 1.9e-10 at 25 to 41 tracked bits (09b2c74
  3.3e-18 to 1.4e-11); and the 0RRC collinear box `4 1.0 50.261921198182165
  0.15614357216746133 82.77728497126748 1.6038370763317959e-13
  0.2547782935219607 -12.782604724181319 0 0 16.407592918401335 0
  82.77728497126748 0 15.201391682789128 -0.00022669601519667755` keeps
  8.3e-9 on the interpreted route at 128 bits in all 24 relabellings, at 30
  tracked bits, where 09b2c74 kept 2.7e-10 (the Symbolica pin: the old pin
  gives 2.9e-10). 256 bits are right on all of them.
- **`Float` keeps the precision it tracks.** Symbolica's `Float` rounds the
  difference of nearly equal operands to the bits it considers significant
  (as in `tests/permutation_regressions.rs`). Next to a shell an invariant
  and a mass differ in their last binary64 bits: in
  `tests/data/two_mass_shell.txt` the native and interpreted results keep 80
  of 128 bits and 208 of 256 (errors 1.5e-23 and 1e-61), within that tracked
  precision. With the Symbolica pin (products rounded once) six
  relabellings of the first row (the sweep regression) track 216 of 256 bits
  but are accurate to 207.7 (3.0e-63), 8.3 bits below what they track, where
  `complex_two_mass_boxes_next_to_their_shells_match_the_references_at_128_and_256_bits`
  allowed 8: that test failed in 12 checks (both routes; the final candidate
  asserts at least `bits - 56` accurate bits instead, the input-cancellation
  floor that it already used for the tracked precision: the tracked
  precision is not a bound, and 09b2c74 itself falls 8.23 to 9.63 bits short
  of it at 128 to 320 bits on row 6, against a 1024-bit reference; the
  16.1 bits at 320 bits on row 0 quoted earlier came from its 90-digit
  reference). With the previous pin
  the same relabellings tracked 216 bits and kept 208.3. The other
  relabellings track 208 or 209 bits and keep 204 to 205. At 192 bits the
  same six fall 8.4 bits short of what they track (previous pin: at most
  6.5); at 128 and 320 bits every relabelling stays within 5.2 and 7.1 bits
  on both pins. numerica's estimate does not follow the amplification of
  rounding errors after these cancellations. The minimum accuracy did not
  regress: at least 202.5 bits at 256 (202.3 at 09b2c74) and 267.0 at 320
  (266.9); single relabellings lose up to 40.3 bits against 09b2c74 (row 0),
  above the floor and the same since f2adce4.
  `tests/float_exact_products.rs` keeps its 8-bit claim: it binds on two of
  its rows (128 and 192 bits) and holds there with at most 4.2 bits between
  tracked and accurate bits; elsewhere its tolerances govern.
- **`Float` rows that the numerica fix leaves wrong** (Float products and
  factors of ±1, above). These are D0 sweep rows where OneLOop is stable and
  binary64 is right. Most other sign and zero tests on cancelled values
  still branch on rounding residues. The difference is that these residues
  come from expressions that are equal only analytically.
  - `4 1.0 0.0 12.25 12.25 -0.04070524207988602 -0.6390609176495762 0 0 0
    2.25 -0.0003109069463992361 4.0 0 2.25 -0.0003109069463992361`: in the 12
    of its 24 relabellings where OneLOop's two precisions agree, both
    `Float` routes give 8.0 + 8.5i at 128 bits (4 tracked bits) and
    8.0 + 8.25i at 192 bits for 8.238 + 8.283i; in the other 12 (OneLOop's
    precisions differ by 6.6e-9) they are off by 6 to 19 times the value
    at both precisions. From 256 bits up all 24 are right. Binary64 loses
    7.8e-12 in the first 12 and 3.5e-9 in the others, as at 09b2c74.
  - `4 1.0 0.01230693334712481 0 0 0 0 2.5453032460105445 0 0 0.06371291035718132
    0 0.06371291035718132 0 4.0 0`, in 4 relabellings: off by 7.0e3
    (interpreted, 192 bits), by 1.0 relative (interpreted, 512 bits) and by
    5.6e4 and 1.1e5 (native and interpreted, 1024 bits).
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
