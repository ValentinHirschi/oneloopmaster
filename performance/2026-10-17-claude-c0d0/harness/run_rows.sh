#!/usr/bin/env bash
set -euo pipefail
source /tmp/claude-1134/-shared-external-ninja/191aef8e-9579-487f-b79d-13d6de4c94fa/scratchpad/devenv.sh >/dev/null 2>&1
S=${S:?scratch directory with devenv.sh, rows_bench_36o3, bench_oracle_36o3 and snap/olm-a7fc9c3}
P=$S/c0d0perf; O=$P/rows; CPU=${CPU:-60}
NEW=/shared/external/ninja/references/oneloopmaster_claude_perf
OLD=$S/snap/olm-a7fc9c3
R=$NEW/performance/2026-10-17-claude
declare -A ROWS=([graph26]=$R/graph26_rows.txt [aaaa]=$R/review/aaaa_rows.txt [gghh]=$R/review/gghh_rows.txt [hgg]=$R/review/hgg_rows.txt)
cd $NEW
for set in graph26 aaaa gghh hgg; do
  taskset -c $CPU $OLD/target/release/examples/rows_survey ${ROWS[$set]} 1024 9 1,1024 --values > $O/${set}_native_before.tsv
  taskset -c $CPU $NEW/target/release/examples/rows_survey ${ROWS[$set]} 1024 9 1,1024 --values > $O/${set}_native_after.tsv
  python3 tests/benchmark_rows.py fortran --executable $S/rows_bench_36o3 --rows ${ROWS[$set]} --calls 1024 --repetitions 9 --cpu $CPU --real --output $O/${set}_fortran_36o3.json
done
declare -A CALLS=([C0]=16384 [D0]=2048)
for f in C0 D0; do
  taskset -c $CPU $OLD/target/release/examples/performance_survey ${CALLS[$f]} 9 $f native 1,1024 > $O/scalar_before_$f.tsv 2> $O/scalar_before_$f.err
  taskset -c $CPU $NEW/target/release/examples/performance_survey ${CALLS[$f]} 9 $f native 1,1024 > $O/scalar_after_$f.tsv 2> $O/scalar_after_$f.err
  taskset -c $CPU python3 tests/benchmark_fortran.py --executable $S/bench_oracle_36o3 --family $f --iterations ${CALLS[$f]} --repetitions 9 --output $O/scalar_fortran_36o3_$f.json > /dev/null
done
echo ROWS-DONE
