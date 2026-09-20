"""Check for the shadowing problem: long-time averages agree across precisions.

Time-averaged kinetic energy over 1e4 time units in single and in double
precision, RK4 with h = 2**-7, from the initial condition of the paper. The two
trajectories part company after about fifty time units; the averages should not.
"""
import json
from pathlib import Path
import numpy as np

def f(y, two):
    t1, t2, w1, w2 = y
    d = t2 - t1
    s, c = np.sin(d), np.cos(d)
    s1, s2 = np.sin(t1), np.sin(t2)
    den = two - c * c
    return np.array([w1, w2,
                     (w1*w1*s*c + s2*c + w2*w2*s - two*s1) / den,
                     (-w2*w2*s*c + two*(s1*c - w1*w1*s - s2)) / den], dtype=y.dtype)

out = {}
for dt in (np.float32, np.float64):
    two, half, six = dt(2), dt(0.5), dt(6)
    h = dt(2.0 ** -7)
    y = np.array([np.pi/2, np.pi/2, 0, 0], dtype=dt)
    acc, n = 0.0, 0
    blocks = []
    for i in range(10_000 * 128):
        k1 = f(y, two); k2 = f(y + half*h*k1, two); k3 = f(y + half*h*k2, two); k4 = f(y + h*k3, two)
        y = y + (h/six)*(k1 + two*(k2 + k3) + k4)
        ke = float(y[2]*y[2] + 0.5*y[3]*y[3] + y[2]*y[3]*np.cos(y[0]-y[1]))
        acc += ke; n += 1
        if n % (1000*128) == 0:
            blocks.append(acc / n)
    e = float(y[2]**2 + 0.5*y[3]**2 + y[2]*y[3]*np.cos(y[0]-y[1]) - 2*np.cos(y[0]) - np.cos(y[1]))
    out[dt.__name__] = {"mean_KE": acc / n, "running": blocks, "final_energy_error": e}
    print(dt.__name__, acc / n, e, flush=True)
Path(__file__).resolve().parents[1].joinpath("data", "shadow_check.json").write_text(json.dumps(out, indent=1))
