#!/usr/bin/env python3
"""Run and compare the mixed-row throughput benchmarks (stdlib only).

Fortran (tests/support/rows_bench_oracle.f90 built against unchanged OneLOop):
  python3 tests/benchmark_rows.py fortran --executable rows_bench --rows ROWS \
      --calls 20480 --repetitions 9 [--cpu 160] [--real] --output fortran.json
Native (examples/rows_survey.rs, release, generated-evaluators):
  taskset -c 160 target/release/examples/rows_survey ROWS 20480 9 1,1024 --values > native.tsv
Compare medians per workload and the VALUE parity of every row:
  python3 tests/benchmark_rows.py compare fortran.json native.tsv [--batch 1024] [--per-point 30]

Each repetition is one timed block of `calls` evaluations; medians are taken over
repetitions. `--per-point N` additionally reports the ALL workload per N rows.
"""

import argparse
import json
import math
from pathlib import Path
import statistics
import subprocess

KINDS = {0: "ALL", 1: "A0", 2: "B0", 3: "C0", 4: "D0", 5: "An", 6: "B11", 7: "B111"}


def numeric_rows(path):
    names, rows = [], []
    name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith("# name "):
            name = line.split()[2]
            continue
        if not line.strip() or line.startswith("#"):
            continue
        rows.append(line.strip())
        names.append(name or f"row_{len(rows)}")
        name = None
    return names, rows


def run_fortran(args):
    names, rows = numeric_rows(args.rows)
    stdin = "OLO_ROWS_BENCH_V1\n"
    stdin += f"{len(rows)} {args.warmup} {args.calls} {args.repetitions} {1 if args.real else 0}\n"
    stdin += "\n".join(rows) + "\n"
    command = [args.executable]
    if args.cpu is not None:
        command = ["taskset", "-c", str(args.cpu)] + command
    stdout = subprocess.run(command, input=stdin, capture_output=True, text=True, check=True).stdout
    values, samples, clock = {}, [], None
    for line in stdout.splitlines():
        words = line.split()
        if not words:
            continue
        if words[0] == "VALUE":
            values[int(words[1])] = [float(v.replace("D", "E")) for v in words[3:]]
        elif words[0] == "CLOCK":
            clock = int(words[1])
        elif words[0] == "SAMPLE":
            mode, kind, fixture, repetition, calls, distinct = words[1], *map(int, words[2:7])
            seconds = float(words[7].replace("D", "E"))
            name = names[fixture - 1] if mode == "SAME" else "all"
            samples.append({"mode": mode, "kind": kind, "name": name, "repetition": repetition,
                            "calls": calls, "distinct": distinct, "ns_per_call": seconds * 1e9 / calls})
    result = {"executable": args.executable, "rows": str(args.rows), "real_interfaces": args.real,
              "cpu": args.cpu, "clock_ticks_per_second": clock, "warmup": args.warmup,
              "calls": args.calls, "repetitions": args.repetitions, "names": names,
              "values": values, "samples": samples}
    Path(args.output).write_text(json.dumps(result, indent=1))
    print(f"wrote {len(samples)} samples to {args.output}")


def medians(samples):
    grouped = {}
    for s in samples:
        grouped.setdefault((s["mode"], s["kind"], s["name"]), []).append(s["ns_per_call"])
    return {key: statistics.median(v) for key, v in grouped.items()}


def compare(args):
    fortran = json.loads(Path(args.fortran).read_text())
    native_samples, native_values = [], {}
    for line in Path(args.native).read_text().splitlines():
        words = line.split("\t")
        if words[0] == "VALUE":
            native_values[int(words[1])] = [float(v) for v in words[3:]]
        elif words[0] in ("SAME", "HETERO", "ALL") and int(words[3]) == args.batch:
            native_samples.append({"mode": words[0], "kind": int(words[1]), "name": words[2],
                                   "ns_per_call": float(words[6])})
    # Value parity: relative difference per complex output, against the row's
    # largest Fortran output of the same coefficient.
    worst = (0.0, None)
    for index, expected in fortran["values"].items():
        actual = native_values.get(int(index))
        if actual is None:
            continue
        for k in range(0, len(expected), 6):
            scale = max(math.hypot(expected[k + 2 * j], expected[k + 2 * j + 1]) for j in range(3)) or 1.0
            for j in range(3):
                a = complex(actual[k + 2 * j], actual[k + 2 * j + 1])
                e = complex(expected[k + 2 * j], expected[k + 2 * j + 1])
                error = abs(a - e) / scale
                if error > worst[0]:
                    worst = (error, f"{fortran['names'][int(index) - 1]} output {k // 2 + j}")
    if native_values:
        print(f"VALUE parity over {len(native_values)} rows: worst relative difference {worst[0]:.3e} at {worst[1]}")
    f, n = medians(fortran["samples"]), medians(native_samples)
    print(f"batch {args.batch}: median ns per call (native / Fortran)")
    print(f"{'mode':7} {'family':5} {'workload':28} {'Fortran':>11} {'native':>11} {'ratio':>7}")
    for key in sorted(f, key=lambda k: (k[0] != "ALL", k[0], k[1], k[2])):
        if key not in n:
            continue
        mode, kind, name = key
        ratio = n[key] / f[key]
        print(f"{mode:7} {KINDS[kind]:5} {name:28} {f[key]:11.1f} {n[key]:11.1f} {ratio:7.3f}")
        if mode == "ALL" and args.per_point:
            print(f"{'':7} {'':5} {'per point':28} {f[key] * args.per_point:11.1f} "
                  f"{n[key] * args.per_point:11.1f} {ratio:7.3f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("fortran")
    run.add_argument("--executable", required=True)
    run.add_argument("--rows", required=True)
    run.add_argument("--calls", type=int, default=20480)
    run.add_argument("--warmup", type=int, default=256)
    run.add_argument("--repetitions", type=int, default=9)
    run.add_argument("--cpu", type=int)
    run.add_argument("--real", action="store_true", help="use real-mass interfaces for real rows")
    run.add_argument("--output", required=True)
    cmp = sub.add_parser("compare")
    cmp.add_argument("fortran")
    cmp.add_argument("native")
    cmp.add_argument("--batch", type=int, default=1)
    cmp.add_argument("--per-point", type=int)
    args = parser.parse_args()
    run_fortran(args) if args.command == "fortran" else compare(args)


if __name__ == "__main__":
    main()
