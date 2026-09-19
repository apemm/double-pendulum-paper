"""Check of the suggested problem on integrator order.

If the global error of an order-n method grows as h**n exp(lambda t), halving
the step should extend the horizon by n ln2 / lambda. RK4 gives n = 4 in the
main text. This script repeats the 237-bit step sweep with forward Euler
(n = 1) and the explicit midpoint rule (n = 2) and measures the same slope.

Reference: the RK4 run at h = 2**-16 from runs.py, whose error is smaller than
that of any run here by more than ten orders of magnitude. Writes
data/orders.json.
"""

from __future__ import annotations

import json
import pickle
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
T_END = 100
T_FIXED = 10     # time at which the error is read off to measure the order
KS = list(range(6, 17))


def integrate(method, k, t_end, prec=237, save_k=4):
    import gmpy2
    from gmpy2 import mpfr
    from dp import rhs

    ctx = gmpy2.get_context()
    ctx.precision = prec
    two = mpfr(2)
    h = mpfr(1) / (1 << k)
    h2 = h / 2
    half_pi = gmpy2.const_pi() / 2
    t1, t2, w1, w2 = half_pi, half_pi, mpfr(0), mpfr(0)
    out = [(t1, t2, w1, w2)]
    for _ in range(t_end << save_k):
        for _ in range(1 << (k - save_k)):
            a1, a2 = rhs(t1, t2, w1, w2, two)
            if method == "euler":
                t1, t2, w1, w2 = t1 + h * w1, t2 + h * w2, w1 + h * a1, w2 + h * a2
            else:  # explicit midpoint
                m1, m2 = w1 + h2 * a1, w2 + h2 * a2
                b1, b2 = rhs(t1 + h2 * w1, t2 + h2 * w2, m1, m2, two)
                t1, t2, w1, w2 = t1 + h * m1, t2 + h * m2, w1 + h * b1, w2 + h * b2
        out.append((t1, t2, w1, w2))
    return out


def job(args):
    method, k = args
    return method, k, integrate(method, k, T_END)


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).parent))
    import runs
    from analysis import log10_bob_distance, horizon, linfit

    with open(ROOT / "data" / "raw" / f"{runs.tag(('h', runs.P_REF, 16, runs.T_MAIN, None))}.pkl", "rb") as fh:
        ref = pickle.load(fh)["traj"]
    t = np.arange(T_END * 16 + 1) / 16.0

    result, early = {}, {}
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 20) as pool:
        for method, k, traj in pool.imap_unordered(job, [(m, k) for m in ("euler", "midpoint") for k in KS]):
            ld = log10_bob_distance(traj, ref)
            result.setdefault(method, {})[k] = horizon(t, ld, 1e-2)
            early.setdefault(method, {})[k] = float(ld[T_FIXED * 16])
            print(method, k, result[method][k], flush=True)

    out = {}
    for method, order in (("euler", 1), ("midpoint", 2)):
        ks = sorted(result[method])
        T = [result[method][k] for k in ks]
        use = [i for i, k in enumerate(ks) if k >= 9 and np.isfinite(T[i])]
        slope, icpt, se = linfit([ks[i] for i in use], [T[i] for i in use])
        # order actually observed: the error at a fixed early time against the step
        ke = [k for k in ks if k >= 9]
        obs, _, obs_se = linfit(ke, [early[method][k] / np.log10(2.0) for k in ke])
        out[method] = {"order": order, "k": ks, "horizon": T, "slope": slope, "se": se,
                       "lambda": order * float(np.log(2)) / slope,
                       "log10_d_fixed_t": [early[method][k] for k in ks],
                       "observed_order": -obs, "observed_order_se": obs_se}
    (ROOT / "data" / "orders.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
