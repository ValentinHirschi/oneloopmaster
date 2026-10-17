#!/usr/bin/env bash
set -euo pipefail
source /tmp/claude-1134/-shared-external-ninja/191aef8e-9579-487f-b79d-13d6de4c94fa/scratchpad/devenv.sh >/dev/null 2>&1
S=/tmp/claude-1134/-shared-external-ninja/191aef8e-9579-487f-b79d-13d6de4c94fa/scratchpad
O=$S/perf
CPU=${CPU:-300}
ROOT=${ROOT:-/shared/external/ninja/references/oneloopmaster_claude_native}
TAG=${TAG:-native}
cd $ROOT
T=$ROOT/target/release/examples
declare -A CALLS=([A0]=524288 [B0]=262144 [dB0]=262144 [C0]=16384 [D0]=2048)
for f in A0 B0 dB0 C0 D0; do
  taskset -c $CPU $T/performance_survey ${CALLS[$f]} 9 $f native 1,1024 > $O/scalar_${TAG}_$f.tsv 2> $O/scalar_${TAG}_$f.err
done
if [ "${FORTRAN:-1}" = 1 ]; then
  for v in o2 o3 36o3; do
    for f in A0 B0 dB0 C0 D0; do
      taskset -c $CPU python3 tests/benchmark_fortran.py --executable $S/bench_oracle_$v --family $f --iterations ${CALLS[$f]} --repetitions 9 --output $O/scalar_fortran_${v}_$f.json > /dev/null
    done
  done
fi
echo SCALAR-DONE
