#!/usr/bin/env python3
"""Write mixed master-row benchmark inputs (stdlib only).

  python3 tests/generate_bench_rows.py tensor > tensor_rows.txt
  python3 tests/generate_bench_rows.py graph26 --pool gg-hhh-changing-pool.json > graph26_rows.txt

Rows are `kind mu_squared arguments...` as read by tests/support/rows_bench_oracle.f90
and examples/rows_survey.rs: kinds 1 A0, 2 B0, 3 C0, 4 D0 (OneLOop argument order,
masses as re/im pairs), 5 An, 6 olo_b11 (B0, B1, B11), 7 olo_bn rank 3 (B0..B111).
A `# name` comment precedes each row.

`graph26` is the gg -> HHH top-loop pentagon of NinjaRS: for each phase-space
point, the 5 boxes, 10 triangles, 10 rank-two bubbles and 5 tadpoles that Ninja
requests for loop shifts V = {0, -p1, -p1-p2, -p1-p2+p3, -p1-p2+p3+p4}, with all
internal squared masses 173^2 and mu^2 from the pool.
"""

import argparse
from itertools import combinations
import json
import sys

TOP = 173.0**2


def tensor_rows():
    real = lambda x: (x, 0.0)
    bubbles = [
        ("p0_equal", 0.0, real(TOP), real(TOP)),
        ("p0_unequal", 0.0, real(1.0), real(2.0)),
        ("massless_m0", 3.0, real(0.0), real(1.0)),
        ("below_threshold", 15625.0, real(TOP), real(TOP)),
        ("above_threshold", 2e5, real(TOP), real(TOP)),
        ("spacelike", -3.0, real(1.0), real(2.0)),
        ("complex", 3.0, (0.7, -0.03), (1.4, -0.08)),
        ("small_p_series", 1e-6, real(1.0), real(2.0)),
    ]
    rows = []
    for kind, label in ((6, "B11"), (7, "B111")):
        for name, p, m0, m1 in bubbles:
            rows.append((f"{label}_{name}", [kind, 1.0, p, *m0, *m1]))
    for name, mass in (("top", real(TOP)), ("negative", real(-2.0)), ("complex", (2.0, -0.1))):
        rows.append((f"An_{name}", [5, 1.0, *mass]))
    return rows


def minkowski(v):
    return v[0] * v[0] - v[1] * v[1] - v[2] * v[2] - v[3] * v[3]


def graph26_rows(pool_path):
    pool = json.load(open(pool_path))
    rows = []
    for point in pool["points"]:
        p = point["momenta"]
        mu2 = float(point["scale_squared"])
        neg = lambda v: [-x for x in v]
        add = lambda a, b: [x + y for x, y in zip(a, b)]
        v0 = [0.0] * 4
        v1 = neg(p[0])
        v2 = add(v1, neg(p[1]))
        v3 = add(v2, p[2])
        v4 = add(v3, p[3])
        shifts = [v0, v1, v2, v3, v4]
        s = lambda i, j: minkowski([a - b for a, b in zip(shifts[j], shifts[i])])
        mass = [TOP, 0.0]
        tag = point["id"]
        for i, j, k, l in combinations(range(5), 4):
            args = [s(i, j), s(j, k), s(k, l), s(l, i), s(i, k), s(j, l)]
            rows.append((f"{tag}_D0_{i}{j}{k}{l}", [4, mu2, *args, *mass * 4]))
        for i, j, k in combinations(range(5), 3):
            rows.append((f"{tag}_C0_{i}{j}{k}", [3, mu2, s(i, j), s(j, k), s(k, i), *mass * 3]))
        for i, j in combinations(range(5), 2):
            rows.append((f"{tag}_B11_{i}{j}", [6, mu2, s(i, j), *mass * 2]))
        for i in range(5):
            rows.append((f"{tag}_A0_{i}", [1, mu2, *mass]))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("set", choices=["tensor", "graph26"])
    parser.add_argument("--pool")
    args = parser.parse_args()
    rows = tensor_rows() if args.set == "tensor" else graph26_rows(args.pool)
    out = sys.stdout
    out.write(f"# {args.set} benchmark rows: kind mu_squared arguments...\n")
    for name, values in rows:
        out.write(f"# name {name}\n")
        out.write(" ".join([str(values[0])] + [repr(float(v)) for v in values[1:]]) + "\n")


if __name__ == "__main__":
    main()
