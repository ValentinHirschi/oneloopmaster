#!/usr/bin/env python3
"""Collect the distinct bubble and tadpole calls that original Ninja makes through OneLOop.

  python3 tests/generate_ninja_tensor_rows.py DIR > tests/data/ninja_tensor_rows.txt

DIR holds `<case>.masters.json` files (NinjaRS validation/benchmark/ninja/):
exact IntegralLibrary call lists recorded from Ninja with OneLOop 3.6 and the
OneLOop results. Output rows (one call each, first-seen order, duplicates
dropped):

  kind mu_squared arguments... (re im of finite, pole, double pole per coefficient)

kind 5 is a tadpole (`olo_a0`; arguments m_re m_im; the A0 triple), 6 a
rank-two bubble and 7 a rank-three bubble (arguments p² m0_re m0_im m1_re
m1_im; the B0, B1 and B11 triples of `olo_b11`, in the native order). For
rank three Ninja builds B111 from these, A0 and B0, so only the OneLOop
coefficients are recorded.
"""
import glob
import json
import os
import sys


def main():
    directory = sys.argv[1]
    seen = set()
    rows = []
    for path in sorted(glob.glob(os.path.join(directory, "*.masters.json"))):
        case = os.path.basename(path)[: -len(".masters.json")]
        if case.startswith("madloop-"):
            continue
        data = json.load(open(path))
        for lane in data["lanes"]:
            for call in lane["calls"]:
                kind = call["kind"]
                result = call["result"]
                masses = call["m2"]
                if kind == "tadpole_cm":
                    family, arguments, triples = 5, [masses[0][0], masses[0][1]], [result[0:3]]
                elif kind in ("rank2_bubble_cm", "rank3_bubble_cm"):
                    # ninja_bench records [b11, b1, b0] and [b111, b11, b1, b0].
                    family = 6 if kind == "rank2_bubble_cm" else 7
                    offset = 0 if family == 6 else 3
                    arguments = [call["s"][0], masses[0][0], masses[0][1], masses[1][0], masses[1][1]]
                    triples = [result[offset + 6 : offset + 9], result[offset + 3 : offset + 6], result[offset : offset + 3]]
                else:
                    continue
                values = [float(call["mu_squared"])] + [float(x) for x in arguments]
                key = (family, tuple(values))
                if key in seen:
                    continue
                seen.add(key)
                outputs = [float(x) for triple in triples for z in triple for x in z]
                rows.append((case, family, values + outputs))
    print("# Distinct bubble/tadpole calls of original Ninja (OneLOop 3.6 results), from tests/generate_ninja_tensor_rows.py.")
    print("# kind(5=A0,6=rank-2,7=rank-3 bubble) mu_squared arguments... OneLOop A0 or B0, B1, B11 as re/im per Laurent order")
    for case, family, values in rows:
        print(f"{family} " + " ".join(repr(v) for v in values) + f"  # {case}")


if __name__ == "__main__":
    main()
