"""Step-by-step record of the 11-bit runs, to pin down why the pendulum stalls.

At p = 11 and h = 2**-12 the angles never change. This prints what each state
variable does at every step until nothing changes any more, and the same for
h = 2**-8, where the pendulum does fall, up to the first step at which theta1
moves. Writes data/lowprec.json.
"""
import json
from pathlib import Path

import gmpy2
from gmpy2 import mpfr

from dp import rhs

OUT = Path(__file__).resolve().parents[1] / "data" / "lowprec.json"


def step(y, h, h2, h6, two):
    t1, t2, w1, w2 = y
    a1, a2 = rhs(t1, t2, w1, w2, two)
    k1 = (w1, w2, a1, a2)
    u1, u2 = w1 + h2 * k1[2], w2 + h2 * k1[3]
    a1, a2 = rhs(t1 + h2 * k1[0], t2 + h2 * k1[1], u1, u2, two)
    k2 = (u1, u2, a1, a2)
    u1, u2 = w1 + h2 * k2[2], w2 + h2 * k2[3]
    a1, a2 = rhs(t1 + h2 * k2[0], t2 + h2 * k2[1], u1, u2, two)
    k3 = (u1, u2, a1, a2)
    u1, u2 = w1 + h * k3[2], w2 + h * k3[3]
    a1, a2 = rhs(t1 + h * k3[0], t2 + h * k3[1], u1, u2, two)
    k4 = (u1, u2, a1, a2)
    return tuple(y[i] + h6 * (k1[i] + two * (k2[i] + k3[i]) + k4[i]) for i in range(4))


def run(p, k, t_max):
    ctx = gmpy2.get_context()
    ctx.precision = p
    ctx.round = gmpy2.RoundToNearest
    two = mpfr(2)
    h = mpfr(1) / (1 << k)
    h2, h6 = h / 2, h / 6
    hp = gmpy2.const_pi() / 2
    y = (hp, hp, mpfr(0), mpfr(0))
    y0 = y
    first_theta_move = None
    last_change = 0
    n_steps = int(t_max * (1 << k))
    trace = []
    for n in range(1, n_steps + 1):
        new = step(y, h, h2, h6, two)
        if new != y:
            last_change = n
        if first_theta_move is None and (new[0] != y0[0] or new[1] != y0[1]):
            first_theta_move = n
            trace.append({"n": n, "t": n / (1 << k), "omega1_before": float(y[2]),
                          "theta1_before": float(y[0]), "theta1_after": float(new[0])})
        y = new
    a1, a2 = rhs(*y, two)
    return {
        "p": p, "k": k, "t_max": t_max,
        "theta1_0": float(y0[0]), "ulp_theta": float(mpfr(2) ** (gmpy2.floor(gmpy2.log2(y0[0])) - (p - 1))),
        "first_theta_move_step": first_theta_move,
        "first_theta_move_t": None if first_theta_move is None else first_theta_move / (1 << k),
        "last_step_at_which_anything_changed": last_change,
        "last_change_t": last_change / (1 << k),
        "final": [float(v) for v in y], "final_accel": [float(a1), float(a2)],
        "trace": trace,
    }


if __name__ == "__main__":
    res = {"k12": run(11, 12, 8), "k08": run(11, 8, 2), "k12_p16": run(16, 12, 2)}
    OUT.write_text(json.dumps(res, indent=1))
    for name, r in res.items():
        print(name, json.dumps({k: v for k, v in r.items() if k != "trace"}))
        print("   trace:", r["trace"])
