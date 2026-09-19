import json, sys
from multiprocessing import Pool
sys.path.insert(0, ".")
def job(a):
    from lowprec_check import run
    p, k = a
    r = run(p, k, 4)
    return p, k, r["first_theta_move_t"], r["final"][2], r["last_change_t"]
if __name__ == "__main__":
    grid = [(p, k) for p in (10, 11, 12, 13, 14, 16) for k in (p - 2, p - 1, p, p + 1)]
    with Pool(12) as pool:
        res = sorted(pool.map(job, grid))
    out = []
    for p, k, tmove, w1, tlast in res:
        print(f"p={p:2d} k={k:2d} (k-p={k-p:+d}): theta first moves at t={tmove}; omega1 at t=4: {w1}; last change t={tlast}")
        out.append({"p": p, "k": k, "theta_first_moves_t": tmove, "omega1_final": w1})
    json.dump(out, open("../data/stall_grid.json", "w"), indent=1)
