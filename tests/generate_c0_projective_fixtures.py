"""Independent massive Euclidean C0 controls using Feynman parameters.

Run with Python and mpmath >=1.3. No oneloop formula, continuation helper,
finite-difference stencil, or symbolic differentiation is used. All input
strings are exact rational numbers; 60/90-digit quadratures must converge.
"""
from fractions import Fraction
from pathlib import Path
import mpmath as mp

BASE = list(map(Fraction, [-1, -2, -3, 2, 3, 5, 7]))
POINTS = [("projective_zero", BASE)]
for parameter in range(6):
    for sign in [-1, 1]:
        point = BASE.copy()
        point[parameter] += Fraction(sign, 64)
        POINTS.append((f"parameter_{parameter}_{sign:+}", point))


def integrate(point, digits):
    with mp.workdps(digits):
        p0, p1, p2, m0, m1, m2, _ = [
            mp.mpf(v.numerator) / v.denominator for v in point
        ]
        def integrand(x, t):
            y = (1 - x) * t
            z = (1 - x) * (1 - t)
            denominator = m0*x + m1*y + m2*z - p0*x*y - p1*y*z - p2*x*z
            return -(1 - x) / denominator
        return mp.quadgl(integrand, [0, 1], [0, 1])


lines = ["# name p0 p1 p2 m0_squared m1_squared m2_squared scale_squared finite"]
for name, point in POINTS:
    low = integrate(point, 60)
    high = integrate(point, 90)
    assert abs(high - low) < mp.mpf("1e-58"), (name, high, low)
    lines.append(" ".join([name, *(str(v) for v in point), mp.nstr(high, 85)]))
Path(__file__).with_name("data").joinpath("c0_projective.tsv").write_text("\n".join(lines) + "\n")
print(f"Generated {len(POINTS)} C0 controls with mpmath {mp.__version__}")
