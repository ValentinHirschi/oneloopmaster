#!/usr/bin/env python3
"""Systematic native-versus-OneLOop sweep over masses, invariants and symmetries.

Usage:
  python3 tests/oracle_sweep.py --oracle /path/to/sweep_oracle \
      --native target/release/examples/oracle_sweep [--seed 1] [--scale 1] \
      [--families A0,B0,dB0,C0,D0] [--jobs 16] [--failures failures.txt]

Build the oracle from tests/support/sweep_oracle.f90 against the unchanged
original OneLOop library and the native evaluator with
`cargo build --release --features generated-evaluators --example oracle_sweep`.

Every generated point is evaluated by OneLOop in double and quadruple precision
and by the native backend in binary64 and at 256 bits. All three Laurent
coefficients are compared, relative to the largest reference coefficient.
A row whose two OneLOop precisions disagree is reported as OneLOop-unstable
and judged against the quadruple-precision and 256-bit results instead.
Symmetric images (bubble mass exchange, triangle rotations/reflections and all
24 box propagator relabelings) are generated for every base point. Points where
OneLOop returns non-finite values or its all-zero error result are counted as
oracle-singular and not compared. The failing rows
can be written in the scalar fixture format, with quadruple-precision
references, for promotion into permanent regression tests.
"""
import argparse
import collections
import itertools
import math
import random
import subprocess
import sys
import threading

KINDS = {"A0": 1, "B0": 2, "dB0": -2, "C0": 3, "D0": 4}
NAMES = {value: key for key, value in KINDS.items()}
SHAPES = {1: (0, 1), 2: (1, 2), -2: (1, 2), 3: (3, 3), 4: (6, 4)}

# Agreement required between native binary64 and stable OneLOop.
TOLERANCE = 1e-10
# Agreement required between the two OneLOop precisions to call it stable.
ORACLE_STABLE = 1e-11
# Agreement between 256-bit native and quadruple-precision OneLOop.
PRECISE = 1e-14


def mass_menu(rng, scale):
    """Squared masses: exact squares (exact thresholds), random, complex."""
    squares = [scale * k * k / 4 for k in (1, 2, 3, 4, 6)]
    value = rng.choice(squares) if rng.random() < 0.5 else scale * 10 ** rng.uniform(-2, 1.5)
    if rng.random() < 0.25:
        width = value * 10 ** rng.uniform(-4, -1)
        return complex(value, -width)
    return complex(value, 0)


def sqrt_mass(mass):
    return math.sqrt(abs(mass.real))


# With --generic, avoid exact zeros, on-shell values and (pseudo-)thresholds.
GENERIC = False


def invariant_menu(rng, a, b, scale):
    """External invariant between propagators with squared masses a and b."""
    ra, rb = sqrt_mass(a), sqrt_mass(b)
    threshold = (ra + rb) ** 2
    pseudo = (ra - rb) ** 2
    choice = rng.randrange(5 if GENERIC else 0, 10)
    if choice == 0:
        return 0.0
    if choice == 1:
        return a.real
    if choice == 2:
        return b.real
    if choice == 3:
        return threshold
    if choice == 4:
        return pseudo
    if choice == 5 and threshold > 0:
        return threshold * rng.uniform(0.05, 0.95)
    if choice == 6:
        return (threshold or scale) * 10 ** rng.uniform(0.02, 2)
    if choice == 7:
        return -(threshold or scale) * 10 ** rng.uniform(-2, 2)
    return scale * rng.choice((-1, 1)) * 10 ** rng.uniform(-2, 2)


def with_pattern(rng, count, scale):
    """A random subset of massive propagators, possibly with equal masses."""
    pattern = rng.randrange(1 << count)
    masses = []
    for index in range(count):
        if pattern & (1 << index):
            if masses and any(m != 0 for m in masses) and rng.random() < 0.4:
                masses.append(rng.choice([m for m in masses if m != 0]))
            else:
                masses.append(mass_menu(rng, scale))
        else:
            masses.append(0j)
    return masses


def a0_points(rng, count, scale):
    for _ in range(count):
        yield 1, scale * rng.choice((1, 2, 0.25)), [], [mass_menu(rng, scale)]
    for mass in (0j, complex(scale, 0), complex(-scale, 0), complex(scale, -0.1 * scale)):
        yield 1, scale, [], [mass]


def bubble_points(kind, rng, count, scale):
    grid = [sign * scale * 10 ** (k / 5) for sign in (-1, 1) for k in range(-15, 16)]
    for _ in range(count):
        m0, m1 = with_pattern(rng, 2, scale)
        if rng.random() < 0.4:
            p = rng.choice(grid)
        else:
            p = invariant_menu(rng, m0, m1, scale)
        yield kind, scale, [p], [m0, m1]
        yield kind, scale, [p], [m1, m0]
    # The reported binary64 failures: one massless line on a logarithmic grid.
    for mass in (1e-3, 0.1, 1.0):
        for p in grid:
            for masses in ([complex(mass * scale), 0j], [0j, complex(mass * scale)]):
                yield kind, scale, [p], masses


def triangle_images(p, m):
    """Rotations and reflections preserving the scalar triangle."""
    images = []
    for shift in range(3):
        pr = [p[(i + shift) % 3] for i in range(3)]
        mr = [m[(i + shift) % 3] for i in range(3)]
        images.append((pr, mr))
        images.append(([pr[2], pr[1], pr[0]], [mr[0], mr[2], mr[1]]))
    return images


def triangle_points(rng, count, scale):
    for _ in range(count):
        m = with_pattern(rng, 3, scale)
        # p1 joins m1,m2; p2 joins m2,m3; p3 joins m3,m1.
        p = [
            invariant_menu(rng, m[0], m[1], scale),
            invariant_menu(rng, m[1], m[2], scale),
            invariant_menu(rng, m[2], m[0], scale),
        ]
        if all(v == 0 for v in p) and all(v == 0 for v in m):
            continue
        for image in triangle_images(p, m):
            yield 3, scale, image[0], image[1]
    # The reported top-loop permutations.
    top = complex(173.0 ** 2, 0)
    for p, m in (([0, 0, 490000.0], [top, 0j, top]), ([0, 0, 90000.0], [top, 0j, top])):
        for image in triangle_images(p, m):
            yield 3, 10000.0, image[0], image[1]


def box_images(p, m):
    """All 24 propagator relabelings of the scalar box (D0 is symmetric)."""
    pair = {
        (0, 1): p[0], (1, 2): p[1], (2, 3): p[2], (0, 3): p[3], (0, 2): p[4], (1, 3): p[5],
    }
    invariant = lambda i, j: pair[(min(i, j), max(i, j))]
    images = []
    for order in itertools.permutations(range(4)):
        a, b, c, d = order
        images.append((
            [invariant(a, b), invariant(b, c), invariant(c, d), invariant(d, a),
             invariant(a, c), invariant(b, d)],
            [m[k] for k in order],
        ))
    return images


def box_points(rng, count, scale):
    for _ in range(count):
        m = with_pattern(rng, 4, scale)
        # p1..p4 join (m1,m2),(m2,m3),(m3,m4),(m4,m1); s12 joins m1,m3; s23 joins m2,m4.
        p = [invariant_menu(rng, m[i], m[(i + 1) % 4], scale) for i in range(4)]
        p.append(invariant_menu(rng, m[0], m[2], scale))
        p.append(invariant_menu(rng, m[1], m[3], scale))
        if p[4] == 0 and p[5] == 0:
            p[4] = -scale * 10 ** rng.uniform(-1, 1)
        for image in box_images(p, m):
            yield 4, scale, image[0], image[1]


def row_text(kind, mu2, p, m):
    fields = [str(kind), repr(float(mu2))]
    fields += [repr(float(v)) for v in p]
    for mass in m:
        fields += [repr(float(mass.real)), repr(float(mass.imag))]
    return " ".join(fields)


def run_chunk(command, rows, prefixes):
    process = subprocess.run(
        command, input="\n".join(rows) + "\n", capture_output=True, text=True, check=True
    )
    lines = [line.split() for line in process.stdout.splitlines()]
    lines = [line for line in lines if line and line[0] in prefixes]
    width = len(prefixes)
    assert len(lines) == width * len(rows), (len(lines), len(rows), process.stderr[-2000:])
    result = []
    for index in range(len(rows)):
        entry = {}
        for line in lines[width * index : width * (index + 1)]:
            if line[1] == "error":
                entry[line[0]] = None
            else:
                values = [float(v) for v in line[1:7]]
                entry[line[0]] = [complex(values[2 * i], values[2 * i + 1]) for i in range(3)]
        result.append(entry)
    return result


def run_parallel(command, rows, prefixes, jobs):
    chunk = max(1, math.ceil(len(rows) / jobs))
    pieces = [rows[i : i + chunk] for i in range(0, len(rows), chunk)]
    results = [None] * len(pieces)

    def work(index):
        results[index] = run_chunk(command, pieces[index], prefixes)

    threads = [threading.Thread(target=work, args=(i,)) for i in range(len(pieces))]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    return [entry for piece in results for entry in piece]


def distance(actual, expected):
    """Largest coefficient error relative to the largest reference magnitude."""
    if actual is None or expected is None:
        return math.inf
    norm = max(abs(v) for v in expected)
    error = 0.0
    for a, e in zip(actual, expected):
        if not (math.isfinite(a.real) and math.isfinite(a.imag)):
            return math.inf
        error = max(error, abs(a - e))
    if norm == 0:
        return error
    return error / norm


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--oracle", required=True)
    parser.add_argument("--native", required=True)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--scale", type=float, default=1.0)
    parser.add_argument("--count", type=int, default=400, help="base points per family")
    parser.add_argument("--families", default="A0,B0,dB0,C0,D0")
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--failures", help="write failing rows in fixture format")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument(
        "--generic", action="store_true", help="no exact zero, on-shell or threshold invariants"
    )
    parser.add_argument("--rows", help="evaluate explicit sweep rows from this file instead")
    args = parser.parse_args()

    global GENERIC
    GENERIC = args.generic
    rng = random.Random(args.seed)
    points = []
    for family in args.families.split(","):
        kind = KINDS[family]
        if kind == 1:
            points += a0_points(rng, args.count, args.scale)
        elif kind in (2, -2):
            points += bubble_points(kind, rng, args.count, args.scale)
        elif kind == 3:
            points += triangle_points(rng, args.count, args.scale)
        else:
            points += box_points(rng, args.count, args.scale)
    texts = [row_text(*point) for point in points]
    if args.rows:
        with open(args.rows) as handle:
            texts = [line.strip() for line in handle if line.strip() and line[0] != "#"]
    seen = set()
    rows = []
    for text in texts:
        if text not in seen:
            seen.add(text)
            rows.append(text)

    oracle = run_parallel([args.oracle], rows, ("dp", "qp"), args.jobs)
    native = run_parallel([args.native], rows, ("f64", "f256"), args.jobs)

    counts = collections.Counter()
    worst = collections.defaultdict(float)
    failing = []
    for row, reference, value in zip(rows, oracle, native):
        family = NAMES[int(row.split()[0])]
        counts[family, "rows"] += 1
        dp, qp, f64, f256 = reference["dp"], reference["qp"], value["f64"], value["f256"]
        if not all(math.isfinite(v.real) and math.isfinite(v.imag) for v in qp) or not any(qp):
            # OneLOop reports a singular or unsupported point (non-finite or
            # its exact-zero error convention); no reference is available.
            counts[family, "oracle-singular"] += 1
            continue
        oracle_error = distance(dp, qp)
        stable = oracle_error <= ORACLE_STABLE
        precise_error = distance(f256, qp)
        if stable:
            error = distance(f64, dp)
            worst[family, "f64"] = max(worst[family, "f64"], error)
            if error <= TOLERANCE:
                counts[family, "pass"] += 1
                continue
            status = "analytic" if precise_error > 1e-9 else "binary64"
        else:
            counts[family, "oracle-unstable"] += 1
            error = distance(f64, qp)
            if precise_error <= 1e-9 and error <= TOLERANCE:
                continue
            status = "unstable-analytic" if precise_error > 1e-9 else "unstable-binary64"
        counts[family, status] += 1
        worst[family, status] = max(worst[family, status], error)
        failing.append((family, status, row, qp, dp, f64, f256, error, precise_error, oracle_error))

    for family in sorted({NAMES[int(row.split()[0])] for row in rows}):
        summary = ", ".join(
            f"{key[1]}={value}" for key, value in sorted(counts.items()) if key[0] == family
        )
        print(f"{family}: {summary}; max stable f64 error {worst[family, 'f64']:.3e}")
    for family, status, row, qp, dp, f64, f256, error, precise_error, oracle_error in failing:
        if args.verbose or status in ("analytic", "binary64"):
            print(
                f"{status:18s} {row}\n    f64 err {error:.2e} f256 err {precise_error:.2e} "
                f"olo dp/qp {oracle_error:.2e}\n    qp   {qp}\n    f64  {f64}\n    f256 {f256}"
            )
    if args.failures:
        with open(args.failures, "w") as handle:
            for family, status, row, qp, *_ in failing:
                if status in ("analytic", "binary64"):
                    references = " ".join(f"{v.real:.17e} {v.imag:.17e}" for v in qp)
                    handle.write(f"# {status}\n{row} {references}\n")
    return 1 if any(f[1] in ("analytic", "binary64") for f in failing) else 0


if __name__ == "__main__":
    sys.exit(main())
