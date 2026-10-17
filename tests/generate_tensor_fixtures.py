#!/usr/bin/env python3
"""Regenerate tests/data/tensor_parity.txt from the original OneLOop Fortran.

Build tests/support/tensor_oracle.f90 against an unchanged double-precision
OneLOop object (see that file), then run:

  python3 tests/generate_tensor_fixtures.py --oracle /path/to/tensor_oracle

Rows are `kind mu_squared arguments... expected...`, with kind 5 (An: m_re m_im;
A0, A00, A0000), 6 (B11: p m0_re m0_im m1_re m1_im; B0, B1, B11) or 7 (B111:
B0, B1, B11, B111). Expected values are real/imaginary pairs of the finite,
simple-pole and double-pole coefficients, in coefficient order.
"""

import argparse
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent.parent
OUTPUT = ROOT / "tests/data/tensor_parity.txt"

TOP = 173.0**2
HIGGS = 125.0**2


def bubbles():
    """(name, p², m0², m1²) with complex masses as (re, im) pairs."""
    real = lambda x: (x, 0.0)
    cases = []
    add = lambda name, p, m0, m1: cases.append((name, p, m0, m1))
    # Zero momentum: equal, unequal and swapped masses, one or both massless.
    add("p0_equal", 0.0, real(1.0), real(1.0))
    add("p0_unequal", 0.0, real(1.0), real(2.0))
    add("p0_unequal_swapped", 0.0, real(2.0), real(1.0))
    add("p0_first_massless", 0.0, real(0.0), real(1.0))
    add("p0_second_massless", 0.0, real(1.0), real(0.0))
    add("p0_equal_complex", 0.0, (1.0, -0.1), (1.0, -0.1))
    add("p0_complex", 0.0, (0.7, -0.03), (1.4, -0.08))
    add("p0_complex_swapped", 0.0, (1.4, -0.08), (0.7, -0.03))
    add("p0_opposite_signs", 0.0, real(-1.0), real(2.0))
    add("p0_opposite_signs_swapped", 0.0, real(2.0), real(-1.0))
    add("p0_near_equal", 0.0, real(1.0), real(1.0 + 1e-9))
    add("p0_series_edge_inside", 0.0, real(0.69), real(1.0))
    add("p0_series_edge_outside", 0.0, real(0.68), real(1.0))
    add("p0_top", 0.0, real(TOP), real(TOP))
    # One massless propagator, including on-shell points.
    add("m0_zero_timelike", 3.0, real(0.0), real(1.0))
    add("m1_zero_timelike", 3.0, real(1.0), real(0.0))
    add("m0_zero_spacelike", -3.0, real(0.0), real(1.0))
    add("m0_zero_below", 0.5, real(0.0), real(1.0))
    add("m0_zero_onshell", 1.0, real(0.0), real(1.0))
    add("m1_zero_onshell", 1.0, real(1.0), real(0.0))
    add("m0_zero_near_onshell", 1.0 + 1e-12, real(0.0), real(1.0))
    add("m0_zero_complex", 2.0, real(0.0), (1.0, -0.1))
    add("m0_zero_complex_shell", 1.0, real(0.0), (1.0, -0.1))
    add("m0_zero_small_p", 1e-7, real(0.0), real(1.0))
    add("m0_zero_negative_mass", 2.0, real(0.0), real(-1.0))
    add("massless_timelike", 3.0, real(0.0), real(0.0))
    add("massless_spacelike", -3.0, real(0.0), real(0.0))
    # Two masses: real and complex, both orderings, all threshold regions.
    add("below_pseudo", 0.1, real(1.0), real(2.0))
    add("pseudo_threshold", (2.0**0.5 - 1.0) ** 2, real(1.0), real(2.0))
    add("between_thresholds", 2.0, real(1.0), real(2.0))
    add("normal_threshold", (2.0**0.5 + 1.0) ** 2, real(1.0), real(2.0))
    add("above_threshold", 7.0, real(1.0), real(2.0))
    add("above_threshold_swapped", 7.0, real(2.0), real(1.0))
    add("spacelike", -3.0, real(1.0), real(2.0))
    add("spacelike_swapped", -3.0, real(2.0), real(1.0))
    add("equal_below", 3.0, real(1.0), real(1.0))
    add("equal_threshold", 4.0, real(1.0), real(1.0))
    add("equal_above", 5.0, real(1.0), real(1.0))
    add("complex", 3.0, (0.7, -0.03), (1.4, -0.08))
    add("complex_swapped", 3.0, (1.4, -0.08), (0.7, -0.03))
    add("complex_spacelike", -3.0, (0.7, -0.03), (1.4, -0.08))
    add("complex_above", 7.0, (1.0, -0.1), (2.0, -0.2))
    add("complex_equal", 4.0, (1.0, -0.1), (1.0, -0.1))
    add("negative_mass", 3.0, real(-1.0), real(2.0))
    add("small_p", 1e-6, real(1.0), real(2.0))
    add("small_p_equal", 1e-6, real(1.0), real(1.0))
    add("small_p_negative", -1e-3, real(1.0), real(1.0))
    add("series_p_0p01", 0.01, real(1.0), real(1.0))
    add("series_p_0p2", 0.2, real(1.0), real(1.0))
    add("series_p_0p5", 0.5, real(1.0), real(1.3))
    add("large_p", 1e8, real(1.0), real(1.0))
    add("large_spacelike", -1e8, real(1.0), real(1.0))
    add("hierarchy_light", 1.0, real(1e-6), real(1.0))
    add("hierarchy_heavy", 1.0, real(1.0), real(1e6))
    # Top-loop kinematics of gg -> H(H)(H).
    add("top_higgs", HIGGS, real(TOP), real(TOP))
    add("top_spacelike", -48125.0, real(TOP), real(TOP))
    add("top_s", 450.0**2, real(TOP), real(TOP))
    add("top_tiny_p", 1e-9, real(TOP), real(TOP))
    return cases


def tadpoles():
    cases = []
    for name, mass in [
        ("positive", (2.0, 0.0)),
        ("negative", (-2.0, 0.0)),
        ("complex", (2.0, -0.1)),
        ("zero", (0.0, 0.0)),
        ("top", (TOP, 0.0)),
        ("tiny", (1e-6, 0.0)),
    ]:
        cases.append((name, mass))
    return cases


def rows():
    out = []
    for mu2 in (1.0, 0.01, 1e4):
        for name, mass in tadpoles():
            out.append((f"an_{name}_mu{mu2:g}", [5, mu2, *mass]))
        for name, p, m0, m1 in bubbles():
            for kind in (6, 7):
                label = "b11" if kind == 6 else "b111"
                out.append((f"{label}_{name}_mu{mu2:g}", [kind, mu2, p, *m0, *m1]))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oracle", required=True)
    parser.add_argument("--output", default=OUTPUT)
    args = parser.parse_args()
    cases = rows()
    encode = lambda values: " ".join([str(values[0])] + [f"{v:.17e}" for v in values[1:]])
    text = "\n".join(encode(values) for _, values in cases) + "\n"
    result = subprocess.run([args.oracle], input=text, capture_output=True, text=True, check=True)
    values = [line.split()[2:] for line in result.stdout.splitlines() if line.startswith("VALUE")]
    if len(values) != len(cases):
        raise SystemExit(f"oracle returned {len(values)} rows for {len(cases)} inputs")
    lines = [
        "# Original OneLOop tensor coefficients (tests/support/tensor_oracle.f90).",
        "# kind(5=An,6=B11,7=B111) mu_squared arguments... expected re/im pairs per coefficient and Laurent order",
    ]
    for (name, inputs), expected in zip(cases, values):
        lines.append(f"# name {name}")
        lines.append(" ".join([encode(inputs)] + [v.replace("D", "E") for v in expected]))
    Path(args.output).write_text("\n".join(lines) + "\n")
    print(f"wrote {len(cases)} rows to {args.output}")


if __name__ == "__main__":
    main()
