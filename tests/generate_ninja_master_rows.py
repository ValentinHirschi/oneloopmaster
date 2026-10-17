#!/usr/bin/env python3
"""Collect the distinct C0/D0 calls that original Ninja makes through OneLOop.

  python3 tests/generate_ninja_master_rows.py DIR > tests/data/ninja_master_rows.txt

DIR holds `<case>.masters.json` files (NinjaRS validation/benchmark/ninja/):
exact IntegralLibrary call lists recorded from Ninja with OneLOop 3.6 and the
OneLOop results. Output rows (one call each, first-seen order, duplicates
dropped):

  family mu_squared momenta... (re im of each squared mass)... (re im of finite, pole, double pole)

family is 3 (C0: p1 p2 p3 m1 m2 m3) or 4 (D0: p1 p2 p3 p4 s12 s23 m1..m4),
arguments in OneLOop order, as the native evaluators take them.
"""
import glob
import json
import os
import sys


def main():
    directory = sys.argv[1]
    seen = set()
    rows = []
    cases = sorted(glob.glob(os.path.join(directory, "*.masters.json")))
    for path in cases:
        case = os.path.basename(path)[: -len(".masters.json")]
        if case.startswith("madloop-"):
            continue
        data = json.load(open(path))
        for lane in data["lanes"]:
            for call in lane["calls"]:
                kind = call["kind"]
                if kind == "box_cm":
                    family, momenta, masses = 4, call["s"][:6], call["m2"][:4]
                elif kind == "triangle_cm":
                    family, momenta, masses = 3, call["s"][:3], call["m2"][:3]
                else:
                    continue
                values = [float(call["mu_squared"])] + [float(x) for x in momenta]
                values += [float(x) for m in masses for x in m]
                key = (family, tuple(values))
                if key in seen:
                    continue
                seen.add(key)
                result = [float(x) for z in call["result"][:3] for x in z]
                rows.append((case, family, values + result))
    print("# Distinct C0/D0 calls of original Ninja (OneLOop 3.6 results), from tests/generate_ninja_master_rows.py.")
    print("# family mu_squared momenta... mass_re mass_im... finite_re finite_im pole_re pole_im double_re double_im")
    for case, family, values in rows:
        print(f"{family} " + " ".join(repr(v) for v in values) + f"  # {case}")


if __name__ == "__main__":
    main()
