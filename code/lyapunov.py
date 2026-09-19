"""Largest Lyapunov exponent by the two-trajectory method of Benettin et al.

A fiducial trajectory and a shadow displaced by d0 are advanced together. Every
tau time units the separation d is recorded, ln(d/d0) is added to a running sum,
and the shadow is pulled back along the line joining the two to distance d0.
The exponent is the running sum divided by the elapsed time.

Done in ordinary double precision with NumPy, vectorised over an ensemble of
initial conditions that share the energy of the run studied in the paper
(E = 0, both rods released from rest). Shadowing is the justification for
using double precision here: the computed orbit is not the true orbit from the
stated initial condition, but it stays close to some true orbit on the same
energy surface, and the exponent is a property of the surface.

Writes data/lyapunov.npz.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "data" / "lyapunov.npz"

H = 1.0 / 128.0
TAU_STEPS = 128            # renormalise once per time unit
T_TOTAL = 20_000           # time units
D0 = 1e-8
N_ENSEMBLE = 48
SEED = 165


def f(y):
    t1, t2, w1, w2 = y
    d = t2 - t1
    s, c = np.sin(d), np.cos(d)
    s1, s2 = np.sin(t1), np.sin(t2)
    den = 2.0 - c * c
    a1 = (w1 * w1 * s * c + s2 * c + w2 * w2 * s - 2.0 * s1) / den
    a2 = (-w2 * w2 * s * c + 2.0 * (s1 * c - w1 * w1 * s - s2)) / den
    return np.array([w1, w2, a1, a2])


def rk4(y, h):
    k1 = f(y)
    k2 = f(y + 0.5 * h * k1)
    k3 = f(y + 0.5 * h * k2)
    k4 = f(y + h * k3)
    return y + (h / 6.0) * (k1 + 2.0 * (k2 + k3) + k4)


def initial_conditions(n, rng):
    """Released from rest with E = -2 cos(theta1) - cos(theta2) = 0.

    The first member is the run studied in the paper, theta1 = theta2 = pi/2.
    The rest take theta1 uniform where a solution exists and solve for theta2.
    """
    t1 = np.empty(n)
    t2 = np.empty(n)
    t1[0] = t2[0] = np.pi / 2
    i = 1
    while i < n:
        a = rng.uniform(np.pi / 3, 2 * np.pi / 3)      # |2 cos a| <= 1
        t1[i] = a * rng.choice([-1, 1])
        t2[i] = np.arccos(-2.0 * np.cos(a)) * rng.choice([-1, 1])
        i += 1
    z = np.zeros(n)
    return np.array([t1, t2, z, z])


def main():
    rng = np.random.default_rng(SEED)
    y = initial_conditions(N_ENSEMBLE, rng)
    u = rng.normal(size=y.shape)
    u /= np.linalg.norm(u, axis=0)
    ys = y + D0 * u

    n_renorm = T_TOTAL
    running = np.zeros((n_renorm, N_ENSEMBLE))
    total = np.zeros(N_ENSEMBLE)
    for i in range(n_renorm):
        for _ in range(TAU_STEPS):
            y = rk4(y, H)
            ys = rk4(ys, H)
        delta = ys - y
        d = np.linalg.norm(delta, axis=0)
        total += np.log(d / D0)
        ys = y + delta * (D0 / d)
        running[i] = total / ((i + 1) * TAU_STEPS * H)
        if (i + 1) % 2000 == 0:
            print(f"t = {i + 1:6d}   lambda[0] = {running[i, 0]:.4f}   "
                  f"ensemble median = {np.median(running[i]):.4f}", flush=True)

    t = (np.arange(n_renorm) + 1) * TAU_STEPS * H
    np.savez_compressed(OUT, t=t[::10], running=running[::10], final=running[-1],
                        h=H, d0=D0, tau=TAU_STEPS * H)
    fin = running[-1]
    print(f"paper IC: {fin[0]:.4f};  ensemble median {np.median(fin):.4f}, "
          f"mean {fin.mean():.4f}, sd {fin.std(ddof=1):.4f}, "
          f"min {fin.min():.4f}, max {fin.max():.4f}")


if __name__ == "__main__":
    main()
