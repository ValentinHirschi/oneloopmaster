#!/usr/bin/env bash
# build.sh OLM_DIR LABEL  -> builds/LABEL/target/release/olbench
set -e
source /tmp/claude-1134/-shared-external-ninja/191aef8e-9579-487f-b79d-13d6de4c94fa/scratchpad/devenv.sh >/dev/null 2>&1
H=/tmp/claude-1134/-shared-external-ninja/191aef8e-9579-487f-b79d-13d6de4c94fa/scratchpad/olbench
B=$H/builds/$2
mkdir -p $B
sed "s|path = \"OLM\"|path = \"$1\"|" $H/Cargo.toml > $B/Cargo.toml
[ -f $B/Cargo.lock ] || cp $H/Cargo.lock $B/Cargo.lock
ln -sfn $H/src $B/src; ln -sf $H/build.rs $B/build.rs; ln -sf $H/olowrap.f90 $B/olowrap.f90
cd $B && CARGO_TARGET_DIR=$B/target cargo build --release --offline -j ${JOBS:-24} 2>&1 | grep -E "^(warning: unused|error|Finished)|-->" | head -30
ls -la $B/target/release/olbench
