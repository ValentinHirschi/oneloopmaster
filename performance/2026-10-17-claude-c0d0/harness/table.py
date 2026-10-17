#!/usr/bin/env python3
"""Combined before/after table over rounds. usage: table2.py base_r1 base_r2 -- after_r1 after_r2"""
import json, sys
args = sys.argv[1:]; k = args.index("--")
load = lambda f: {c["case"]: c["timing_ns_per_point"] for c in json.load(open(f))["cases"]}
B = [load(f) for f in args[:k]]; A = [load(f) for f in args[k+1:]]
order = ["hgg-top-triangle-0", "gg-hh-top-triangle-0", "gg-hhh-top-triangle-0", "four-photon-top-box-0",
         "gg-hh-top-box-2", "gg-hhh-top-box-8", "gg-hhh-top-pentagon-26"]
order += sorted(c for c in A[0] if c not in order)
def rng(vals, fmt="{:.2f}"):
    lo, hi = min(vals), max(vals)
    return fmt.format(lo) if fmt.format(lo) == fmt.format(hi) else fmt.format(lo) + "–" + fmt.format(hi)
def rat(runs, case, g):
    if g not in runs[0][case]: return "–"
    return rng([r[case][g]["ratio"] for r in runs])
print("| case | OneLOop µs/point | native µs/point before → after | all: before → after | D0: before → after | C0: before → after | bubbles | A0 |")
print("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
for case in order:
    ol = rng([r[case]["all"]["oneloop_ns"] / 1e3 for r in A + B])
    nb = rng([r[case]["all"]["native_ns"] / 1e3 for r in B]); na = rng([r[case]["all"]["native_ns"] / 1e3 for r in A])
    def col(g):
        if g not in A[0][case]: return "–"
        return f"{rat(B, case, g)} → **{rat(A, case, g)}**"
    print(f"| {case} | {ol} | {nb} → {na} | {col('all')} | {col('D0')} | {col('C0')} | {rat(A, case, 'B')} | {rat(A, case, 'A0')} |")
