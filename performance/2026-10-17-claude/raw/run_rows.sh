#!/usr/bin/env bash
set -euo pipefail
source /tmp/claude-1134/-shared-external-ninja/191aef8e-9579-487f-b79d-13d6de4c94fa/scratchpad/devenv.sh >/dev/null 2>&1
S=/tmp/claude-1134/-shared-external-ninja/191aef8e-9579-487f-b79d-13d6de4c94fa/scratchpad
O=$S/perf
CPU=${CPU:-300}
cd /shared/external/ninja/references/oneloopmaster_claude_native
T=$PWD/target/release/examples
for set in graph26 tensor; do
  taskset -c $CPU $T/rows_survey $S/${set}_rows.txt 1024 9 1,1024 --values > $O/${set}_native.tsv
  for v in o2 o3 36o3; do
    python3 tests/benchmark_rows.py fortran --executable $S/rows_bench_$v --rows $S/${set}_rows.txt --calls 1024 --repetitions 9 --cpu $CPU --real --output $O/${set}_fortran_$v.json
  done
done
echo ROWS-DONE
