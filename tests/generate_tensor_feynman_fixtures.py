#!/usr/bin/env python3
"""Regenerate tests/data/tensor_feynman.txt: independent bubble references.

The finite parts of the longitudinal bubble coefficients are integrated
directly in Feynman-parameter space with mpmath (45 significant digits):

  B_{1^n}(p², m0², m1²) = (-1)^(n+1) int_0^1 x^n log((D(x) - i0)/mu²) dx,
  D(x) = p² x² + (m1² - m0² - p²) x + m0²,

for n = 0..3 (B0, B1, B11, B111; OneLOop normalization, denominators
(q²-m0²)((q+p)²-m1²)). This shares nothing with OneLOop's root/bnlog
algorithm. The rows probe small momenta, normal and pseudo thresholds and
their neighbourhoods, massless propagators (on and near shell), complex
masses, large mass ratios, real negative masses and the gg -> HHH top loop.

  python3 tests/generate_tensor_feynman_fixtures.py [--jobs 16]

Requires mpmath. Rows: `# name <label>` then
`p² m0_re m0_im m1_re m1_im mu² B0_re B0_im B1_re B1_im B11_re B11_im B111_re B111_im`.
"""

import argparse
import math
import multiprocessing
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parent.parent
OUTPUT = ROOT / "tests/data/tensor_feynman.txt"
TOP, HIGGS = 173.0**2, 125.0**2


def cases():
    rows = []
    R = lambda x: (float(x), 0.0)

    def add(name, p, m0, m1, mu2=1.0):
        rows.append((name, float(p), m0, m1, float(mu2)))

    small = [1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16, 1e-20, 1e-300]
    for m0, m1, tag in [
        (R(1), R(1), "equal"), (R(1), R(2), "unequal"), (R(2), R(1), "swapped"),
        (R(0), R(1), "m0_zero"), (R(1), R(0), "m1_zero"), ((1, -0.1), (1, -0.1), "complex_equal"),
        ((0.7, -0.03), (1.4, -0.08), "complex"), (R(TOP), R(TOP), "top"),
        (R(1), R(1 + 1e-9), "near_equal"), (R(1e-8), R(1), "hierarchy"),
    ]:
        for e in small:
            for s in (1, -1):
                add(f"small_p_{tag}_{s * e:g}", s * e * (TOP if tag == "top" else 1.0), m0, m1)
    for a, b in [(1.0, 4.0), (4.0, 1.0), (1.0, 2.0), (TOP, TOP), (1.0, 1.0), (1e-6, 1.0), (9.0, 16.0)]:
        for sign, tag in [(1, "normal"), (-1, "pseudo")]:
            threshold = (math.sqrt(a) + sign * math.sqrt(b)) ** 2
            if threshold == 0.0:
                continue
            for d in [0.0, 1e-15, -1e-15, 1e-12, -1e-12, 1e-8, -1e-8, 1e-4, -1e-4, 1e-1, -1e-1]:
                add(f"{tag}_threshold_{a:g}_{b:g}_{d:g}", threshold * (1 + d), R(a), R(b))
    for m in [1.0, 2.0, 49.0, TOP]:
        for m0, m1, tag in [(R(0), R(m), "m0_zero"), (R(m), R(0), "m1_zero")]:
            for f in [1.0, 1 + 1e-15, 1 - 1e-15, 1 + 1e-10, 1 - 1e-10, 2.0, 0.5, -1.0, 1e6, 1 / 7.0]:
                add(f"{tag}_{m:g}_x{f!r}", m * f, m0, m1)
        add(f"p0_m0_zero_{m:g}", 0.0, R(0), R(m))
        add(f"p0_m1_zero_{m:g}", 0.0, R(m), R(0))
    for p in [1.0, -1.0, 3.0, 49.0, -49.0, 98.0, 1e-10, -1e-10, 1e10, -1e10, 1e-300, 1e300]:
        add(f"massless_{p:g}", p, R(0), R(0))
    width = (TOP, -173.0 * 1.5)
    z_boson = (91.1876**2, -91.1876 * 2.4952)
    for m0, m1, tag in [
        (width, width, "top_width"), (z_boson, z_boson, "z_width"), (R(0), width, "zero_top_width"),
        (width, R(0), "top_width_zero"), (z_boson, width, "z_top"), (R(TOP), width, "top_top_width"),
        ((1.0, -1e-12), (1.0, -1e-12), "tiny_width"), ((1.0, -1e-12), (4.0, -1e-12), "tiny_width_14"),
    ]:
        for p in [0.0, HIGGS, -HIGGS, 4 * TOP, (math.sqrt(m0[0]) + math.sqrt(m1[0])) ** 2,
                  (math.sqrt(m0[0]) - math.sqrt(m1[0])) ** 2, 1e-6, 1e8, m1[0]]:
            add(f"complex_{tag}_{p:g}", p, m0, m1)
    for a, b in [(1e-12, 1.0), (1.0, 1e-12), (1e-6, 1e6), (1e6, 1e-6), (1.0, 1e12), (1e12, 1.0)]:
        for p in [1e-12, 1e-6, 1.0, -1.0, 1e6, 1e12, -1e12]:
            add(f"ratio_{a:g}_{b:g}_{p:g}", p, R(a), R(b))
    for a, b in [(-1.0, 2.0), (2.0, -1.0), (-1.0, -1.0), (-1.0, 0.0), (0.0, -1.0), (-2.0, -1.0)]:
        for p in [0.0, 1e-8, 3.0, -3.0, 0.5]:
            add(f"negative_{a:g}_{b:g}_{p:g}", p, R(a), R(b))
    for p in [0.0, HIGGS, -48125.0, 450.0**2, 600.0**2, -104375.0, 135625.0, 83125.0, 1e-9]:
        add(f"top_loop_{p:g}", p, R(TOP), R(TOP), mu2=1e4)
    seen, unique = set(), []
    for row in rows:
        key = tuple(float(v).hex() for v in (row[1], *row[2], *row[3], row[4]))
        if key not in seen:
            seen.add(key)
            unique.append(row)
    return unique


def log_lip(v):
    """log(v - i0): real negative arguments take -i pi.

    A node landing exactly on a (double) root of D carries a negligible
    tanh-sinh weight; give it 0 instead of -inf.
    """
    if v == 0:
        return mp.mpf(0)
    if mp.im(v) == 0 and mp.re(v) < 0:
        return mp.log(-mp.re(v)) - 1j * mp.pi
    return mp.log(v)


def reference(row):
    mp.mp.dps = 45
    _, p, m0, m1, mu2 = row
    p, mu2 = mp.mpf(p), mp.mpf(mu2)
    m0, m1 = mp.mpc(*m0), mp.mpc(*m1)
    if p == 0 and m0 == 0 and m1 == 0:
        return " ".join(["0.0 0.0"] * 4)
    if p != 0:
        # Factor D = p (x - r1)(x - r2) with stably computed roots, so that
        # D keeps its relative precision next to (double) roots.
        b, c = m1 - m0 - p, m0
        root = mp.sqrt(b * b - 4 * p * c)
        q = -(b + root) / 2 if mp.re(mp.conj(b) * root) >= 0 else -(b - root) / 2
        roots = [q / p, c / q] if q != 0 else [mp.mpf(0), mp.mpf(0)]
        if mp.im(m0) == 0 and mp.im(m1) == 0:
            # Real D: drop the rounding noise of complex-conjugate roots,
            # which would otherwise pick a lip at random.
            D = lambda x: mp.re(p * (x - roots[0]) * (x - roots[1]))
        else:
            D = lambda x: p * (x - roots[0]) * (x - roots[1])
        candidates = roots
    else:
        D = lambda x: (m1 - m0) * x + m0
        candidates = [-m0 / (m1 - m0)] if m1 != m0 else []
    points = [mp.mpf(0)]
    points += sorted({mp.re(x) for x in candidates if 0 < mp.re(x) < 1})
    points.append(mp.mpf(1))
    values = [
        (-1) ** (n + 1) * mp.quad(lambda x: x**n * log_lip(D(x) / mu2), points, maxdegree=12)
        for n in range(4)
    ]
    # Format here: mpmath numbers unpickled in the parent round to its precision.
    text = lambda x: mp.nstr(x, 40, min_fixed=1, max_fixed=0)
    return " ".join(f"{text(mp.re(v))} {text(mp.im(v))}" for v in values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--output", default=OUTPUT)
    args = parser.parse_args()
    rows = cases()
    with multiprocessing.Pool(args.jobs) as pool:
        values = pool.map(reference, rows, chunksize=4)
    lines = [
        "# Feynman-parameter bubble references (tests/generate_tensor_feynman_fixtures.py, mpmath).",
        "# p2 m0_re m0_im m1_re m1_im mu2 then re/im of the finite parts of B0 B1 B11 B111",
    ]
    for (name, p, m0, m1, mu2), result in zip(rows, values):
        lines.append(f"# name {name}")
        inputs = " ".join(repr(float(v)) for v in (p, *m0, *m1, mu2))
        lines.append(f"{inputs} {result}")
    Path(args.output).write_text("\n".join(lines) + "\n")
    print(f"wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
