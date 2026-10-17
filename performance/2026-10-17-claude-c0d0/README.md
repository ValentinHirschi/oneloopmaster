# Hand-written binary64 C0/D0 kernels — 2026-10-17

Gate: oneloopmaster's native (direct Rust, `NativeEvaluator<f64>`) masters may
take at most 1.3× the time of OneLOop Fortran for the same calls, per diagram
and per family, on the master-call streams of the benchmark diagrams, without
losing accuracy. Before this change D0 (1.4–2.2×) and several C0 workloads
(up to 2.1×) failed it.

**Result.** All 7 physical diagrams now run at 0.24–0.57× of OneLOop 3.6 in
total, D0 at 0.53–0.64× and C0 at 0.26–0.69×. The synthetic streams are at
0.44–1.02× in total; the only remaining per-family misses are the two-mass
triangles of `syn-n4-*` (C0 1.14–1.42×), which still use the generated code
(see *Remaining misses*).

## What changed

`src/native/handwritten/` adds binary64 kernels for the dominant sectors,
dispatched from `NativeEvaluator::evaluate` before the generated code:

| Kernel | Expression it transcribes | OneLOop analogue | Used for |
| --- | --- | --- | --- |
| `triangle.rs` | `triangle_finite_three_masses` (`select`, `first`, `second`, `ordered`) | `trif3` | C0, three massive propagators, nonzero momenta, Källén ≥ 0 |
| `box_contour.rs` | all of `box_complex.rs` (`box`, `general`, `ordered`, `zero_external`, `opposite*`, `tfun`, `tfun_selected`, `t1`, `t13`, `s3`, `s3_double_real`, `r0`, `r1`, `planar`, `inside`, `distance`, `root_1/2`, `sort`) | `boxc` | D0 with four, or three (finite sector) massive propagators where `finite_box_representation` selects the contour box |
| `boxes.rs` | `box_series` channel exchange, `select_four_series` sectors 7/11/13/14/15, `box_finite_four_masses` | `d0cc` dispatch, `boxf4` | D0 dispatch; four masses outside the contour kinematics |
| `sheet.rs` | `sheet_exact.rs` (products, quotients, `log`, `log_over_one_minus`, continued `Li2(1-x)`, divided difference) | `qmplx`, `li2c2` | all of the above |

Each kernel is a branch-for-branch transcription: the same exact zero tests,
the same `x - |x|` sign tests (also for ±∞/NaN), infinitesimal directions,
root labels and sheets. What makes them 2–3.5× faster than the generated code
is structural: one sector dispatch per call instead of one per Laurent
coefficient, inactive branches never evaluated, real predicates tested on real
values instead of
`x + (-1)·|x|` in complex arithmetic, no duplicated 755-line bodies or
84-argument helpers, and values computed once (e.g. the `beta_data` of both
contour candidates is reused for the selected one). Two deliberate numerical
differences, identical in exact arithmetic: the three-mass triangle takes a
quadratic root that cancels in `(b ± d)/(2a)` from the product of roots
(OneLOop's `solabc`), so that a root rounded to an exact zero cannot reach a
logarithm (the generated code returns NaN at such points, e.g.
`C0(0.2827, 0, 3; 0.2827, 0.3891, 0.1441)`); and arithmetic is not reassociated
the way Symbolica normalizes the expressions, so results differ at the rounding
level.

`DoubleFloat` and `Float` keep the generated code: their tracked precision
follows the expression's grouping, and an exact zero that the expression
obtains only through symbolic simplification would otherwise lose all tracked
bits. A hidden `NativeEvaluator::set_handwritten(bool)` selects either code for
validation. Regenerating the expressions does not update these kernels; an
expression change in these sectors (such as the pending fixes on
`claude/expression-fixes`) has to be transcribed too.

> **Update (expression fixes):** the expressions now take the cancelling
> triangle root from the product of roots as well, so the generated code no
> longer differs there, and the kernels were re-transcribed with the fixed
> expressions (`trif3` roots, chart choice and root lips, `boxf4` root lips,
> the contour box's double-root split and its exactly real residue
> logarithms, exact r-function roots at zero momenta and thresholds).

## Per diagram and per family (Ninja's master-call streams)

µs per phase-space point (one lane's calls, lanes cycled) on exactly the calls
original Ninja makes (`validation/benchmark/ninja/<case>.masters.json` of
NinjaRS), native vs OneLOop 3.6 in the same process, interleaved block by
block. Native and its tensor families as in NinjaRS's `masters_native`
harness (`B` = `B11`/`B111` tensors for Ninja's rank-2/3 bubbles, compared with
`olo_b11` and Ninja's rank-3 construction). Ranges over two rounds; before =
`b39241f` (native C0/D0 code identical to `origin/main` `a7fc9c3`), after =
`2c14e1a` (measured from a snapshot with an identical `src/`).

| case | OneLOop µs/point | native µs/point before → after | all: before → after | D0: before → after | C0: before → after | bubbles | A0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| hgg-top-triangle-0 | 6.55–6.87 | 3.18 → 1.67–1.68 | 0.48–0.49 → **0.24–0.25** | – | 0.58–0.59 → **0.26–0.28** | 0.14 | 0.62 |
| gg-hh-top-triangle-0 | 1.93–1.98 | 1.84–1.85 → 0.79–0.86 | 0.95–0.96 → **0.40–0.44** | – | 1.63–1.64 → **0.51–0.59** | 0.24 | 0.62 |
| gg-hhh-top-triangle-0 | 1.30–1.34 | 1.73 → 0.70–0.77 | 1.33 → **0.52–0.57** | – | 2.06–2.07 → **0.58–0.69** | 0.39 | 0.62 |
| four-photon-top-box-0 | 10.14–10.21 | 15.55–15.62 → 5.35–5.56 | 1.53–1.54 → **0.53–0.54** | 2.17–2.18 → **0.61–0.62** | 1.39–1.40 → **0.46–0.51** | 0.31 | 0.62 |
| gg-hh-top-box-2 | 27.22–27.39 | 32.53–32.76 → 13.15–13.23 | 1.19–1.20 → **0.48–0.49** | 1.44–1.46 → **0.53–0.54** | 1.18–1.19 → **0.49–0.53** | 0.17 | 0.62 |
| gg-hhh-top-box-8 | 22.16–22.23 | 31.79–31.96 → 12.34–12.48 | 1.43–1.44 → **0.56** | 1.75–1.77 → **0.62–0.63** | 1.40–1.41 → **0.54–0.58** | 0.20 | 0.62 |
| gg-hhh-top-pentagon-26 | 91.64–92.38 | 127.31–127.58 → 51.16–51.36 | 1.38 → **0.56** | 1.67–1.68 → **0.63–0.64** | 1.02–1.04 → **0.47–0.48** | 0.17 | 0.60–0.62 |
| syn-gghhh-gHgHH-n5-r4 | 61.22–61.66 | 55.66–55.77 → 26.81–26.89 | 0.90–0.91 → **0.44** | 1.12 → **0.51** | 0.89–0.90 → **0.43** | 0.19 | 0.62 |
| syn-gghhh-gHgHH-n5-r5 | 59.44–59.80 | 55.81–55.87 → 26.66–26.71 | 0.93–0.94 → **0.45** | 1.17 → **0.52–0.53** | 0.91–0.92 → **0.44** | 0.19 | 0.62 |
| syn-gghhh-ggHHH-n5-r4 | 110.12–110.78 | 130.91–131.19 → 56.44–56.77 | 1.18–1.19 → **0.51–0.52** | 1.49–1.51 → **0.60** | 0.82–0.83 → **0.42–0.43** | 0.16 | 0.62 |
| syn-gghhh-ggHHH-n5-r5 | 109.31–109.93 | 131.77–132.06 → 56.45–56.88 | 1.20 → **0.51–0.52** | 1.51–1.54 → **0.60–0.61** | 0.82–0.83 → **0.42–0.43** | 0.16 | 0.62 |
| syn-gghhh-nosmat-n5-r4 | 110.29–110.95 | 130.24–131.01 → 56.31–56.66 | 1.17–1.18 → **0.51** | 1.50–1.51 → **0.61** | 0.79–0.80 → **0.40–0.41** | 0.16 | 0.62 |
| syn-n2-r1 | 0.20–0.22 | 0.20 → 0.20 | 0.98 → **0.91–0.92** | – | – | 1.06–1.08 | 0.61–0.62 |
| syn-n2-r2 | 0.20–0.21 | 0.21 → 0.21 | 1.07–1.08 → **1.00–1.02** | – | – | 1.21–1.22 | 0.62–0.63 |
| syn-n2-r3 | 0.32–0.35 | 0.22 → 0.21 | 0.67 → **0.62–0.64** | – | – | 0.63–0.64 | 0.62 |
| syn-n3-r2 | 4.26–4.28 | 3.69 → 2.30–2.33 | 0.86–0.87 → **0.54–0.55** | – | 1.10 → **0.63** | 0.30 | 0.65 |
| syn-n3-r3 | 3.87–3.88 | 3.72 → 2.30–2.33 | 0.96 → **0.59–0.60** | – | 1.27–1.28 → **0.72–0.73** | 0.30 | 0.64–0.65 |
| syn-n3-r4 | 5.02–5.04 | 3.69 → 2.33–2.35 | 0.73 → **0.46–0.47** | – | 0.96–0.97 → **0.55–0.56** | 0.24 | 0.64–0.65 |
| syn-n4-r3 | 16.81–16.94 | 19.65–19.76 → 11.90–11.99 | 1.16–1.18 → **0.71** | 1.27–1.28 → **0.46** | 1.36–1.37 → **1.14–1.15** | 0.24 | 0.62 |
| syn-n4-r4 | 15.86–15.97 | 24.94–25.18 → 13.55–13.66 | 1.57–1.58 → **0.85–0.86** | 1.69–1.71 → **0.59** | 1.65 → **1.41–1.42** | 0.46 | 0.62 |
| syn-n4-r5 | 12.20–12.38 | 18.59–18.79 → 10.66–10.73 | 1.52–1.54 → **0.87** | 1.92–1.94 → **0.62** | 1.69–1.70 → **1.41–1.42** | 0.28 | 0.62 |
| syn-n5-r4 | 65.49–65.65 | 84.81–85.01 → 35.00 | 1.29–1.30 → **0.53** | 1.55–1.56 → **0.59–0.60** | 1.09–1.10 → **0.49–0.50** | 0.22 | 0.62 |
| syn-n5-r5 | 71.49–71.67 | 83.87–83.91 → 35.63–35.81 | 1.17 → **0.50** | 1.39–1.40 → **0.54–0.55** | 1.01–1.03 → **0.48** | 0.20 | 0.62 |
| syn-n5-r6 | 85.67–85.95 | 103.99–104.01 → 44.79–45.07 | 1.21 → **0.52–0.53** | 1.55–1.56 → **0.61–0.62** | 0.93–0.94 → **0.47** | 0.17–0.18 | 0.62 |
| syn-n6-r5 | 481.28–487.26 | 586.35–587.64 → 302.20–304.85 | 1.20–1.21 → **0.63** | 1.26–1.27 → **0.60–0.61** | 1.06–1.07 → **0.80–0.81** | 0.26 | 0.64–0.65 |
| syn-n6-r6 | 420.38–422.38 | 504.42–506.33 → 258.24–260.57 | 1.20 → **0.61–0.62** | 1.25–1.26 → **0.57** | 1.01 → **0.76–0.77** | 0.25 | 0.63–0.65 |
| syn-n6-r7 | 534.99–536.14 | 613.53–613.86 → 317.30–318.10 | 1.14–1.15 → **0.59** | 1.23 → **0.58–0.59** | 0.83–0.85 → **0.63–0.65** | 0.21 | 0.64–0.67 |

The OneLOop side reproduces NinjaRS's `ninja_bench` masters-only modes
(`validation/benchmark/ninja/SUMMARY.md`, CPU 29) within 1–3 % (e.g. four-photon
10.13, pentagon 91.38 µs).

## Repository row and fixture protocol

The same comparison with this repository's drivers (`examples/rows_survey.rs`,
`examples/performance_survey.rs` against `tests/support/rows_bench_oracle.f90` /
`bench_oracle.f90` built from OneLOop 3.6 with `-O3 -fno-fast-math`, as in
[the review re-measurement](../2026-10-17-claude/README.md#review-re-measurement-independent-builds-three-more-diagrams)),
nine blocks of at least 25 ms (rows) or 20 ms (fixtures), CPU 60, before =
`a7fc9c3`, after = `2c14e1a`. Native / Fortran, batch 1 (batch 1,024 in
`raw/`):

| Workload | Fortran | before | after |
| --- | ---: | ---: | ---: |
| gg → HHH pentagon (graph26 rows), all | 3,427 ns/call (102.8 µs/point) | 1.250 | **0.578** |
| — D0 / C0 per call | 12,475 / 3,247 ns | 1.548 / 0.886 | **0.672 / 0.480** |
| γγ → γγ top box (aaaa rows), all | 567 ns/call | 1.722 | **0.759** |
| — D0 / C0 per call | 3,227 / 910 ns | 2.145 / 1.859 | **0.848 / 0.801** |
| gg → HH top box (gghh rows), all | 2,106 ns/call | 1.083 | **0.527** |
| — D0 / C0 per call | 15,775 / 2,877 ns | 1.415 / 0.924 | **0.629 / 0.492** |
| H → gg top triangle (hgg rows), all | 469 ns/call | 0.608 | **0.343** |
| — C0 per call | 1,506 ns | 1.070 | **0.512** |
| scalar fixtures, mixed D0 | 22,119 ns | 1.196 | **0.571** |
| — `D0_real_timelike` | 1,434 ns | 2.133 | **0.874** |
| — `D0_real_euclidean` / `D0_complex` | 5,956 / 41,789 ns | 0.688 / 1.214 | **0.362 / 0.576** |
| scalar fixtures, mixed C0 | 2,960 ns | 0.842 | **0.772** |
| — `C0_real_timelike` | 2,476 ns | 1.111 | **0.556** |

Native values agree with the Fortran values of every row to 2.2e-14 (`VALUE
parity` lines in `raw/*_comparison_*`). `C0_real_euclidean` and `C0_complex`
(Källén < 0, `trif3HV`) and the other families are unchanged.

## Accuracy

`tests/native_handwritten.rs` (part of the suite):

- **Ninja's streams.** Every distinct C0/D0 call of the 27 recorded streams
  (`tests/data/ninja_master_rows.txt`, 495 C0 and 288 D0 rows, generated by
  `tests/generate_ninja_master_rows.py`), each Laurent coefficient. In the
  hand-written sectors (376 C0, 286 D0 rows), scaled to the largest
  coefficient, median / worst error:

  | | vs OneLOop 3.6 | vs the expressions at 256/1024 bits (C0 232, D0 222 rows that keep ≥ 60 tracked bits) |
  | --- | --- | --- |
  | C0 hand-written | 1.7e-15 / 8.9e-13 | 9.1e-16 / 8.9e-13 |
  | C0 generated | 2.2e-15 / 6.1e-12 | 1.0e-15 / 5.7e-12 |
  | D0 hand-written | 4.1e-15 / 6.8e-10 | 5.2e-15 / 6.8e-10 |
  | D0 generated | 4.3e-15 / 6.8e-10 | 5.4e-15 / 6.8e-10 |

  The 6.8e-10 is one `syn-gghhh-ggHHH-n5-r5` box where both binary64
  evaluations of the expression lose digits (OneLOop agrees with the
  high-precision expression to 4e-12 there). The `olbench check` mode gives
  the same per-case maxima against OneLOop before and after (e.g. pentagon D0
  1.7e-14 → 1.9e-14, C0 8.5e-15 → 1.1e-14), and the native results are
  invariant under the 8 dihedral box and 6 triangle relabellings to ≤ 8e-13
  (`syn-n6-*`; ≤ 1e-13 elsewhere) except the known two-mass rows below.
- **Randomized kinematics.** 3,000 triangles and 600 boxes with "nice" values
  (exact zeros, equal masses, thresholds, complex widths, one massless line),
  binary64 against the generated expressions at 256 (or 512) bits: no gross
  deviation (> 1e-7 and more than 16× the generated code's error) and no
  systematic loss (hand-written more than 16× less accurate than the generated
  binary64 code on 0 points, more accurate on 4); the generic triangle kernel
  also agrees with the 256-bit expressions at 256 bits (worst 1e-38 where both
  keep ≥ 120 tracked bits).
- **Contour box definitions.** `box_contour.rs`'s unit test evaluates each
  ported definition against its registered `__olo_boxc_*` expression
  (Symbolica's binary64 evaluator) on 1,300–3,000 random, often degenerate,
  arguments each, deciding disagreements with the expression at 256 bits:
  no failure for `root_1/2`, `r0`, `r1`, `s3`, `planar`, `t1`, `t13`; for
  `box` 1 of 2,999 points (a three-mass box whose binary64 sign test is
  decided by roundoff — the generated binary64 code gives the same value as the
  port there), while 3 points of `box` and 2 of `t1` are cases where the port
  agrees with the 256-bit expression and Symbolica's binary64 evaluation does
  not.

Unchanged known expression defects (same values before and after): the
two-mass triangles with one massless line (`syn-n4-r3/r4`, `syn-n6-r6`, scaled
deviation 1.5 from OneLOop), the four-mass real box with `s12·s23 < 0`
(`tests/native_d0_four_mass_defect.rs`, reproduced by the DNS kernel since it
transcribes `box_finite_four_masses`).

## Test suite

`cargo test --release --features generated-evaluators -- --test-threads=1`
(`RUST_MIN_STACK=134217728`): 185 passed, 0 failed, 8 ignored, over every
test target (including the native, tensor, parity and hook tests, the new
`native_handwritten` target and the `box_contour` unit test).

## Remaining misses

- `syn-n4-r3/r4/r5` C0: 1.14–1.42× (diagram totals 0.71–0.87×). These are
  two-mass triangles (one massless internal line), still generated code: their
  expression (`triangle_finite_two_masses`) is being repaired on
  `claude/expression-fixes` (NaN and wrong imaginary parts, new root lips), and
  transcribing the current version would freeze the defect. A first
  transcription was 3× faster but, like the current expression, depends on
  exact zeros produced by roundoff (e.g. `a = r34/r24 - r23` at zero momenta),
  so it waits for the repaired expression.
- `syn-n2-r2` bubbles: 1.21× (the tensor `B11` family, unchanged; total 1.00–1.02×).

## Protocol and provenance

- Host: shared AMD EPYC 9754, logical CPU 60 (SMT sibling 316 and the rest of
  its L3 complex idle when measured), no builds running. Rust 1.98, GNU Fortran
  15.3.0.
- Streams: `harness/` (`olbench`), built with NinjaRS's release profile (thin
  LTO, 1 CGU) against `git archive` snapshots of the measured commits and the
  Symbolica snapshot `e679e8c1`. OneLOop 3.6 is `libavh_olo.a` of NinjaRS's
  `o3-fortranrules` Ninja build (sha256 `502929c4…`, gfortran `-O3 -fPIC`),
  called through `bind(c)` wrappers identical to Ninja's `avholo_wrapper.F90`
  (`olo_d0`/`olo_c0`/`olo_b11`/`olo_a0`/`olo_b0` with complex masses and
  `rmu = sqrt(mu²)`), `olo_onshell(1e-9)` as Ninja's `AvHOneLoop::init`, warnings
  off. Each group (all calls in order, then each family alone) calibrates
  blocks of at least 25 ms, runs one warm-up block per side and nine rounds of
  one native and one OneLOop block; medians. Raw: `raw/time_*`.
- Rows and fixtures: `harness/run_rows.sh`, raw TSV/JSON and comparison tables
  in `raw/`. The harness scripts reference the local scratch paths used here
  (`build.sh` builds `olbench` against a oneloopmaster directory; `table.py`
  renders the stream table from the `raw/time_*.json` files).
