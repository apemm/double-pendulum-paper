"""Run every simulation the paper uses and cache the raw trajectories.

Each job is one fixed-step RK4 integration at a given significand width p and
step 2**-k. Trajectories are pickled as tuples of mpfr values so that the
differences between runs, which start out near 1e-70, can be formed at full
precision later. A job whose output file already exists is skipped, so the
script can be interrupted and restarted.

Usage:  python runs.py [n_workers]
"""

from __future__ import annotations

import pickle
import sys
import time
from multiprocessing import Pool
from pathlib import Path

RAW = Path(__file__).resolve().parents[1] / "data" / "raw"

T_MAIN = 160      # time units; 51 s for a 1 m pendulum
T_PREC = 240      # the precision sweep at fixed step can be followed for longer
P_REF = 237       # significand width of IEEE binary256
K_REF = 18        # reference step 2**-18 = 3.8e-6

# A perturbation of 2**-170 = 6.7e-52 in one coordinate stays in the linear
# regime for the whole of T_MAIN, so no renormalisation is needed to read the
# Lyapunov exponent off the separation.
PERTURB_EXP = 170


def jobs():
    out = []
    # (a) step-size sweep at 237 bits: truncation error only
    for k in range(4, K_REF + 1):
        out.append(("h", P_REF, k, T_MAIN, None))
    # (b) precision sweep at fixed step: roundoff only, reference is p = 237
    for k in (8, 12, 16):
        for p in (11, 16, 20, 24, 28, 32, 36, 40, 48, 53, 64, 80, 96, 113, 128, P_REF):
            out.append(("p", p, k, T_PREC, None))
    # (c) step-size sweep at working precisions: both errors at once
    for p in (24, 32, 40, 53, 64):
        for k in range(4, 18):
            out.append(("h", p, k, T_MAIN, None))
    # (d) perturbed initial conditions for a direct measurement of lambda
    for which in (0, 3):
        out.append(("ly", P_REF, 14, T_MAIN, which))
    # longest first, so the pool is never waiting on one straggler at the end
    out.sort(key=lambda j: -(j[3] << j[2]) * (1.0 + j[1] / 200.0))
    return out


def tag(job):
    kind, p, k, t_end, which = job
    base = f"{kind}_p{p:03d}_k{k:02d}_T{t_end}"
    return base if which is None else f"{base}_c{which}"


def run(job):
    import gmpy2
    from gmpy2 import mpfr
    import dp

    kind, p, k, t_end, which = job
    path = RAW / f"{tag(job)}.pkl"
    if path.exists():
        return tag(job), 0.0
    y0 = None
    if which is not None:
        gmpy2.get_context().precision = p
        half_pi = gmpy2.const_pi() / 2
        y0 = [half_pi, half_pi, mpfr(0), mpfr(0)]
        y0[which] = y0[which] + mpfr(1) / (mpfr(2) ** PERTURB_EXP)
    t0 = time.perf_counter()
    traj = dp.integrate(p, k, t_end, y0=y0)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "wb") as fh:
        pickle.dump({"kind": kind, "p": p, "k": k, "t_end": t_end,
                     "which": which, "save_k": dp.SAVE_K, "traj": traj}, fh)
    tmp.replace(path)
    return tag(job), time.perf_counter() - t0


if __name__ == "__main__":
    RAW.mkdir(parents=True, exist_ok=True)
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    todo = [j for j in jobs() if not (RAW / f"{tag(j)}.pkl").exists()]
    print(f"{len(todo)} jobs on {workers} workers", flush=True)
    t0 = time.perf_counter()
    with Pool(workers) as pool:
        for name, secs in pool.imap_unordered(run, todo):
            print(f"{time.perf_counter() - t0:8.0f}s  {name}  ({secs:.0f}s)", flush=True)
    print("done", flush=True)
