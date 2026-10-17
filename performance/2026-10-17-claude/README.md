# Binary64 native masters and tensor coefficients — 2026-10-17

Matched warm binary64 measurements of the direct Rust (`native`) backend against
original OneLOop Fortran, after the changes on `claude/native-masters`:

1. **Dilogarithm.** `NativeFloat for f64` calls Symbolica's new
   `transcendental::dilog_complex_f64` (a Bernoulli-series kernel with Symbolica's
   lips and signed zeros, [PR89](https://github.com/symbolica-dev/symbolica/pull/89))
   directly. The pinned integration `e679e8c1` is the previous pin `2c06f245`
   plus that kernel. Previously every binary64 Li2 went through the registered
   `polylog(2,z)` callback, whose binary64 route runs MPFR with a 64-bit guard.
2. **Reciprocal and norm.** Off the axes and inside `[2^-500, 2^500]`, binary64
   complex reciprocals use one division and `hypot` uses `sqrt(x²+y²)`; axes,
   non-finite values and extreme scales keep the scaled algorithms.
3. **Tensor coefficients.** `NativeTensorEvaluator` families `An`
   (`A0, A00, A0000`), `B11` (`B0, B1, B11`) and `B111` (`B0, B1, B11, B111`),
   a port of OneLOop's `tadpn`/`bub*` routines (see the
   [README](../../README.md#tensor-coefficients)).

## Headline: the gg → HHH top-loop pentagon (graph26)

Per phase-space point, Ninja requests 5 D0, 10 C0, 10 rank-two bubbles
(`olo_b11` ↔ `B11`) and 5 A0 for the loop shifts
`{0, -p1, -p1-p2, -p1-p2+p3, -p1-p2+p3+p4}`, all internal masses `m² = 173²`,
`mu² = 10⁴`, at the three points of `gg-hhh-changing-pool.json`
([rows](graph26_rows.txt), from `tests/generate_bench_rows.py graph26`). Fortran
uses the real-mass interfaces, as Ninja's `*_rm` wrappers do. Times are medians
of nine timed blocks of at least 25 ms over the cycled 90-row stream (`ALL`);
per point = 30 calls.

| Fortran build | Fortran µs/point | native µs/point (batch 1 / 1024) | native / Fortran |
| --- | ---: | ---: | ---: |
| OneLOop 2024 snapshot, `-O2 -fno-fast-math` | 104.9 | 128.0 / 126.8 | **1.221 / 1.209** |
| OneLOop 2024 snapshot, `-O3 -fno-fast-math` | 104.9 | 128.0 / 126.8 | **1.221 / 1.209** |
| OneLOop 3.6, `-O3 -fno-fast-math` | 103.4 | 128.0 / 126.8 | **1.238 / 1.226** |

The graph26 workload is within the 1.3× gate against all three Fortran builds.
Per family on the same rows (batch 1, ns per call, Fortran 2024 `-O3`):

| Family | Fortran | native | ratio | share of native point |
| --- | ---: | ---: | ---: | ---: |
| D0 (5/point) | 12,879 | 19,250 | 1.495 | 75% |
| C0 (10/point) | 3,256 | 2,861 | 0.879 | 22% |
| B11 vs `olo_b11` (10/point) | 734 | 131 | 0.178 | 1% |
| A0 (5/point) | 55.5 | 43.9 | 0.791 | <1% |

The boxes are the remaining gap: removing the dilogarithm work entirely would
save only ~1.6 µs of the 19.3 µs D0, and all logarithms ~2.6 µs (stub
experiments on these rows); the rest is the generated D0 arithmetic and
branching. The rank-two bubbles save more than the C0 surplus costs.

**Correctness of the reference.** Every graph26 row agrees with OneLOop 3.6 to
1.4e-14 (relative, worst output). The 2024 OneLOop snapshot used by this
repository's oracle recipe (`/common/dev/oneloop/OneLOopBridge/oneloop/avh_olo.f90`,
dated 2024-11-23) returns a different real part for one box,
`sqrt_s_600_D0_1234` = `D0(0, 15625, 15625, -104375, -104375, 135625; 173²×4)`:
`-2.0606e-9 + 1.8575e-10 i`, also in its quad precision. OneLOop 3.4 (bundled
with Ninja), OneLOop 3.6, COLLIER (COLI and DD) and this crate at 53, 128, 256 and
512 bits give `2.8638e-10 + 1.8575e-10 i`; `tests/native_graph26.rs` pins it.

## Scalar families (standard benchmark fixtures)

`performance_survey ... native 1,1024` against `tests/benchmark_fortran.py` with
`tests/support/bench_oracle.f90`, nine repetitions, calls per family chosen for
blocks of at least 20 ms (A0 524,288; B0/dB0 262,144; C0 16,384; D0 2,048).
Mixed-input (`HETERO`) medians, ns per call, native / Fortran:

| Family | 2024 `-O2`, b1 | 2024 `-O3`, b1 | 3.6 `-O3`, b1 | 2024 `-O3`, b1024 | 3.6 `-O3`, b1024 |
| --- | ---: | ---: | ---: | ---: | ---: |
| A0 | 0.824 | 0.825 | 0.849 | 0.705 | 0.725 |
| B0 | 0.821 | 0.825 | 0.859 | 0.738 | 0.769 |
| dB0 | 0.883 | 0.902 | 0.914 | 0.869 | 0.880 |
| C0 | 0.847 | 0.858 | 0.871 | 0.828 | 0.841 |
| D0 | 1.110 | 1.111 | 1.184 | 1.122 | 1.195 |

Absolute (batch 1, 2024 `-O3`): A0 55.4 vs 45.7, B0 132.2 vs 109.0,
dB0 166.6 vs 150.2, C0 2,990 vs 2,565, D0 23,574 vs 26,201 ns. All mixed
workloads are within 1.3×. Single-row misses above 1.3× (all three Fortran
builds): `D0_real_timelike` 2.15–2.18× (1.40 vs 3.03 µs), `A0_scaleless`
1.44–1.65× (10–11 vs 16–17 ns, validation and output overhead), and at batch 1
`dB0_one_massless_onshell` 1.35× (43 vs 59 ns). Full tables:
`scalar_comparison_<fortran>_batch<b>.txt`.

### Before/after on the same fixtures

The previous `main` (`c9dd82f`, Symbolica `2c06f245`), same protocol except the
D0 run (64 calls, three repetitions), mixed workloads at batch 1:

| Family | `c9dd82f` | this branch | speed-up |
| --- | ---: | ---: | ---: |
| A0 | 50.8 ns | 45.7 ns | 1.11× |
| B0 | 115.4 ns | 109.0 ns | 1.06× |
| dB0 | 159.1 ns | 150.2 ns | 1.06× |
| C0 | 1,721,503 ns | 2,565 ns | 671× |
| D0 | 8,171,690 ns | 26,201 ns | 312× |

On the previous pin each binary64 dilogarithm cost 33–250 µs (MPFR), so C0 and
D0 were dominated by Li2 (graph26 boxes evaluate about 50 per call). The full
baseline C0/D0 survey was stopped after C0 batch 1 (about 70 minutes per family).

## Tensor families

Rows from `tests/generate_bench_rows.py tensor` ([rows](tensor_rows.txt)):
eight bubble kinematics (zero momentum, one massless propagator, below/above
threshold with top masses, spacelike, complex masses, the small-`p²` series)
and three tadpoles; Fortran calls `olo_b11` (which also computes `B00`),
`olo_bn` at rank 3 and `olo_an` at rank 4. Mixed-row medians, ns per call:

| Family | Fortran 2024 `-O3` | Fortran 3.6 `-O3` | native b1 / b1024 | ratio (2024, b1) |
| --- | ---: | ---: | ---: | ---: |
| `An` vs `olo_an` | 72.5 | 73.7 | 37.0 / 37.2 | 0.510 |
| `B11` vs `olo_b11` | 424.9 | 417.4 | 126.7 / 132.0 | 0.298 |
| `B111` vs `olo_bn(3)` | 526.9 | 529.3 | 142.3 / 148.8 | 0.270 |

Per-row ratios are 0.12–0.63 (`tensor_comparison_*.txt`). Native outputs agree
with the Fortran values of these rows to 9.2e-16.

## Protocol and provenance

- Host: shared AMD EPYC 9754 (load about 40 from other users), every timed
  process pinned to logical CPU 300, no builds running. Rust 1.98.1, release
  profile; GNU Fortran 15.3.0. Fortran object files are unchanged OneLOop
  sources compiled with `-O2 -fno-fast-math` or `-O3 -fno-fast-math`:
  the 2024 snapshot (sha256 `196e3533…`) and OneLOop 3.6
  (`MG5_aMC_v3_7_0/HEPTools/oneloop/OneLOop-3.6`, sha256 `f3cd1a8a…`). Ninja's
  bundled 3.4 archive was used for the value check only (its module aliasing
  does not link with the benchmark drivers).
- Native: `examples/rows_survey.rs` (graph26, tensor rows) and
  `examples/performance_survey.rs` (scalar fixtures), route `native`, `f64`,
  all outputs, batches 1 and 1,024, warmup 256 calls. Fortran:
  `tests/support/rows_bench_oracle.f90` via `tests/benchmark_rows.py fortran
  --real` and `tests/support/bench_oracle.f90` via `tests/benchmark_fortran.py`.
  Comparisons: `tests/benchmark_rows.py compare` and `tests/compare_performance.py
  --maximum-ratio 1.3`. Raw TSV/JSON in [`raw/`](raw/); the two driver scripts
  there reference local scratch paths.
- Each SAME/HETERO/ALL workload calibrates its block size to at least 25 ms
  (rows drivers) before nine timed repetitions; medians are reported.
