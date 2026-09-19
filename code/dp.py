"""Fixed-step RK4 for the planar double pendulum in arbitrary binary precision.

Units. Equal masses, equal rod lengths l, and time measured in units of
sqrt(l/g), so the equations of motion carry no dimensional parameters. For a
pendulum with l = 1 m and g = 9.8 m/s^2 one time unit is 0.3194 s.

Precision. Every floating-point operation is carried out by MPFR (through
gmpy2) with a p-bit significand and round-to-nearest, so p = 24, 53, 64 and 113
reproduce the significand widths of IEEE binary32, binary64, x87 extended and
binary128. The exponent range is MPFR's, which is never approached here, so
only the significand width is being varied.

Step sizes are powers of two, h = 2**-k. A power of two is represented exactly
at every precision, so lowering p changes the arithmetic and nothing else; a
decimal step such as 1e-4 would be rounded differently at each precision and
the runs would be integrating with slightly different steps. It also puts the
output times, which are multiples of 2**-4, exactly on the grid of every run.

State vector y = (theta1, theta2, omega1, omega2), angles measured from the
downward vertical.
"""

from __future__ import annotations

import gmpy2
from gmpy2 import mpfr, sin_cos, sin, cos

SAVE_K = 4  # output every 2**-SAVE_K time units


def rhs(t1, t2, w1, w2, two):
    """Angular accelerations for m1 = m2, l1 = l2, g = l = 1."""
    s, c = sin_cos(t2 - t1)
    s1 = sin(t1)
    s2 = sin(t2)
    den = two - c * c
    w1s = w1 * w1
    w2s = w2 * w2
    a1 = (w1s * s * c + s2 * c + w2s * s - two * s1) / den
    a2 = (-w2s * s * c + two * (s1 * c - w1s * s - s2)) / den
    return a1, a2


def energy(t1, t2, w1, w2):
    """Total energy in units of m g l. Zero for both rods horizontal at rest."""
    return (w1 * w1 + w2 * w2 / 2 + w1 * w2 * cos(t1 - t2)
            - 2 * cos(t1) - cos(t2))


def integrate(prec: int, k: int, t_end: int, y0=None, save_k: int = SAVE_K):
    """RK4 with step 2**-k and p-bit arithmetic from t = 0 to t_end.

    Returns a list of (theta1, theta2, omega1, omega2) tuples of mpfr values
    sampled every 2**-save_k time units, including t = 0. Conversion to double for plotting is left to the caller so that no
    precision is lost before differences are taken.
    """
    if k < save_k:
        raise ValueError("step must not exceed the output interval")
    ctx = gmpy2.get_context()
    ctx.precision = prec
    ctx.round = gmpy2.RoundToNearest

    two = mpfr(2)
    h = mpfr(1) / (1 << k)          # exact
    h2 = h / 2                      # exact
    h6 = h / 6                      # rounded once, at precision p

    if y0 is None:
        half_pi = gmpy2.const_pi() / 2
        t1, t2, w1, w2 = half_pi, +half_pi, mpfr(0), mpfr(0)
    else:
        t1, t2, w1, w2 = (mpfr(v) for v in y0)

    per_save = 1 << (k - save_k)
    n_save = t_end << save_k
    out = [(t1, t2, w1, w2)]

    for _ in range(n_save):
        for _ in range(per_save):
            a1, a2 = rhs(t1, t2, w1, w2, two)
            k1t1, k1t2, k1w1, k1w2 = w1, w2, a1, a2

            u1, u2 = w1 + h2 * k1w1, w2 + h2 * k1w2
            a1, a2 = rhs(t1 + h2 * k1t1, t2 + h2 * k1t2, u1, u2, two)
            k2t1, k2t2, k2w1, k2w2 = u1, u2, a1, a2

            u1, u2 = w1 + h2 * k2w1, w2 + h2 * k2w2
            a1, a2 = rhs(t1 + h2 * k2t1, t2 + h2 * k2t2, u1, u2, two)
            k3t1, k3t2, k3w1, k3w2 = u1, u2, a1, a2

            u1, u2 = w1 + h * k3w1, w2 + h * k3w2
            a1, a2 = rhs(t1 + h * k3t1, t2 + h * k3t2, u1, u2, two)

            t1 = t1 + h6 * (k1t1 + two * (k2t1 + k3t1) + u1)
            t2 = t2 + h6 * (k1t2 + two * (k2t2 + k3t2) + u2)
            w1 = w1 + h6 * (k1w1 + two * (k2w1 + k3w1) + a1)
            w2 = w2 + h6 * (k1w2 + two * (k2w2 + k3w2) + a2)
        out.append((t1, t2, w1, w2))
    return out
