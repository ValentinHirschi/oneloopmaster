#!/usr/bin/env python3
"""Random physical D0 next to a vanishing Gram determinant, with OneLOop references.

Usage:
  python3 tests/generate_near_threshold_rows.py /path/to/sweep_oracle [--seed 1] [--count 200] \
      > tests/data/near_threshold_boxes.txt

Build the oracle from tests/support/sweep_oracle.f90 against the unchanged
original OneLOop library (see tests/oracle_sweep.py). Standard library only.

Each box has real external momenta: a timelike q decays into three particles,
two of which (a, b) move nearly together, their invariant a relative r
(1e-12 to 1e-2, log-uniform) above its threshold (m_a + m_b)². The two momenta
are then nearly parallel, so the Källén function of the pair and the Gram
determinant of the box nearly vanish (gg -> HHH with two Higgs bosons nearly at
rest relative to each other, as in tests/d0_vanishing_gram.rs). The four legs
are shuffled, so the near-threshold invariant is an external mass, s12 or s23.
Internal masses: four equal real masses, four random real masses, three equal
masses and a massless line, four equal complex masses, or random real and
complex masses. Scales 1e-2, 1 and 1e4.

Rows are printed in the scalar fixture order of tests/data/permutation_regressions.txt:
  4 mu² p1² p2² p3² p4² s12 s23 (Re m² Im m²)x4 (Re Im)x3
with OneLOop 3.6's quadruple-precision finite part and poles as references.
Rows where OneLOop's quadruple precision is non-finite or its zero error value
are dropped (none so far).
"""
import argparse
import math
import random
import subprocess


def boost(p, b):
    bx, by, bz = b
    b2 = bx * bx + by * by + bz * bz
    if b2 == 0:
        return list(p)
    g = 1 / math.sqrt(1 - b2)
    bp = bx * p[1] + by * p[2] + bz * p[3]
    g2 = (g - 1) / b2
    return [
        g * (p[0] + bp),
        p[1] + g2 * bp * bx + g * bx * p[0],
        p[2] + g2 * bp * by + g * by * p[0],
        p[3] + g2 * bp * bz + g * bz * p[0],
    ]


def direction(rng):
    while True:
        v = [rng.uniform(-1, 1) for _ in range(3)]
        n = math.sqrt(sum(x * x for x in v))
        if 0.1 < n <= 1:
            return [x / n for x in v]


def square(p):
    return p[0] * p[0] - p[1] * p[1] - p[2] * p[2] - p[3] * p[3]


def add(p, q):
    return [x + y for x, y in zip(p, q)]


def decay(rng, total, ma, mb, mc, k):
    """q (mass `total`) at rest -> a + b + c, with a and b at relative momentum k."""
    ea, eb = math.sqrt(ma * ma + k * k), math.sqrt(mb * mb + k * k)
    mab = ea + eb
    if mab + mc >= total:
        return None
    q = math.sqrt((total**2 - (mab + mc) ** 2) * (total**2 - (mab - mc) ** 2)) / (2 * total)
    n = direction(rng)
    pab = [math.sqrt(mab * mab + q * q)] + [q * x for x in n]
    pc = [math.sqrt(mc * mc + q * q)] + [-q * x for x in n]
    d = direction(rng)
    a = [ea] + [k * x for x in d]
    b = [eb] + [-k * x for x in d]
    beta = [pab[i] / pab[0] for i in (1, 2, 3)]
    return boost(a, beta), boost(b, beta), pc


def boxes(seed, count):
    rng = random.Random(seed)
    while count:
        scale = rng.choice([1.0, 1.0, 1e-2, 1e4])
        mh = math.sqrt(scale) * rng.uniform(0.5, 2)
        ma = mh
        mb = mh * rng.choice([1, 1, rng.uniform(0.3, 2)])
        mc = mh * rng.choice([1, rng.uniform(0.3, 2)])
        total = (ma + mb + mc) * rng.uniform(1.05, 3)
        r = 10 ** rng.uniform(-12, -2)
        point = decay(rng, total, ma, mb, mc, math.sqrt(r * ma * mb))
        if point is None:
            continue
        a, b, c = point
        legs = [[-(a[i] + b[i] + c[i]) for i in range(4)], a, b, c]
        kind = rng.randrange(5)
        top = math.sqrt(scale) * rng.uniform(0.6, 2)
        if kind == 0:
            m = [complex(top * top, 0)] * 4
        elif kind == 1:
            m = [complex((math.sqrt(scale) * rng.uniform(0.3, 2)) ** 2, 0) for _ in range(4)]
        elif kind == 2:
            m = [complex(top * top, 0)] * 4
            m[rng.randrange(4)] = 0j
        elif kind == 3:
            m = [complex(top * top, -top * top * 10 ** rng.uniform(-4, -1))] * 4
        else:
            m = [
                complex(
                    (math.sqrt(scale) * rng.uniform(0.3, 2)) ** 2,
                    -(scale * 10 ** rng.uniform(-4, -1)) * rng.choice([0, 1]),
                )
                for _ in range(4)
            ]
        rng.shuffle(legs)
        p = [square(x) for x in legs]
        p.append(square(add(legs[0], legs[1])))
        p.append(square(add(legs[1], legs[2])))
        # Widths of exactly real masses are +0, as OneLOop's callers pass them.
        fields = ["4", repr(float(scale))] + [repr(float(v)) for v in p]
        fields += [f"{float(z.real)!r} {float(z.imag) + 0.0!r}" for z in m]
        yield " ".join(fields)
        count -= 1


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("oracle")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--count", type=int, default=200)
    args = parser.parse_args()
    rows = list(boxes(args.seed, args.count))
    process = subprocess.run(
        [args.oracle], input="\n".join(rows) + "\n", capture_output=True, text=True, check=True
    )
    lines = [line.split() for line in process.stdout.splitlines()]
    quadruple = [line[1:7] for line in lines if line and line[0] == "qp"]
    assert len(quadruple) == len(rows), (len(quadruple), len(rows))
    print("# D0 next to a vanishing Gram determinant (tests/generate_near_threshold_rows.py,")
    print(f"# seed {args.seed}, {args.count} boxes); references: OneLOop 3.6 in quadruple precision.")
    print("# kind mu_squared p1 p2 p3 p4 s12 s23 (m_re m_im)x4 finite_re finite_im pole_re pole_im double_re double_im")
    for row, reference in zip(rows, quadruple):
        values = [float(v) for v in reference]
        if not all(math.isfinite(v) for v in values) or not any(values):
            continue
        print(row + " " + " ".join(f"{v:.17e}" for v in values))


if __name__ == "__main__":
    main()
