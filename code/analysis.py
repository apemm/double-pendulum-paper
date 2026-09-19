"""Turn cached trajectories into error curves, horizons and fitted exponents.

Differences between runs are formed with 320-bit arithmetic before anything is
converted to double, because the early separations are of order 1e-70 and would
vanish in a double-precision subtraction.

The error measure is the distance between the lower bobs of two pendulums, in
units of the rod length, so it saturates near the maximum possible value of 4.

Writes data/analysis.npz (curves) and data/summary.json (every number quoted in
the papers).
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import gmpy2
import numpy as np
from gmpy2 import mpfr, sin, cos, sqrt, log10

import runs

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
LN2 = float(np.log(2.0))
BS, NL = chr(92), chr(10)

THRESHOLDS = (1e-6, 1e-4, 1e-2, 1e-1)
MAIN_THRESHOLD = 1e-2


def load(kind, p, k, t_end, which=None):
    path = RAW / f"{runs.tag((kind, p, k, t_end, which))}.pkl"
    with open(path, "rb") as fh:
        return pickle.load(fh)["traj"]


def log10_bob_distance(a, b):
    """log10 of the lower-bob separation between two trajectories."""
    gmpy2.get_context().precision = 320
    n = min(len(a), len(b))
    out = np.full(n, -np.inf)
    for i in range(n):
        a1, a2 = mpfr(a[i][0]), mpfr(a[i][1])
        b1, b2 = mpfr(b[i][0]), mpfr(b[i][1])
        dx = (sin(a1) + sin(a2)) - (sin(b1) + sin(b2))
        dy = (cos(a1) + cos(a2)) - (cos(b1) + cos(b2))
        d = sqrt(dx * dx + dy * dy)
        if d > 0:
            out[i] = float(log10(d))
    return out


def log10_state_distance(a, b):
    gmpy2.get_context().precision = 320
    n = min(len(a), len(b))
    out = np.full(n, -np.inf)
    for i in range(n):
        s = mpfr(0)
        for u, v in zip(a[i], b[i]):
            e = mpfr(u) - mpfr(v)
            s += e * e
        if s > 0:
            out[i] = float(log10(sqrt(s)))
    return out


def energy_error(traj):
    import dp
    gmpy2.get_context().precision = 320
    e0 = dp.energy(*[mpfr(v) for v in traj[0]])
    return np.array([float(abs(dp.energy(*[mpfr(v) for v in y]) - e0)) for y in traj])


def running_max(a):
    return np.maximum.accumulate(np.where(np.isfinite(a), a, -400.0))


def horizon(t, log10_d, threshold):
    """First time the separation reaches the threshold; nan if it never does."""
    hit = np.flatnonzero(log10_d >= np.log10(threshold))
    return float(t[hit[0]]) if hit.size else float("nan")


def linfit(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    A = np.vstack([x, np.ones_like(x)]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    n = len(x)
    s2 = float(res[0]) / (n - 2) if n > 2 and res.size else 0.0
    se = float(np.sqrt(s2 / ((x - x.mean()) ** 2).sum())) if n > 2 else float("nan")
    return float(coef[0]), float(coef[1]), se


# Ranges used in the straight-line fits. Coarse steps are left out because RK4
# is not yet in its asymptotic regime there; the finest is left out because its
# error is no longer small compared with the reference's. Low precisions are
# left out because the arithmetic stalls rather than perturbs (see text); high
# ones because their horizon lies beyond the end of the run.
K_FIT = list(range(7, 17))
P_FIT = [24, 28, 32, 36, 40, 48, 53, 64, 80]


def fits(S):
    th = str(MAIN_THRESHOLD)
    F = {}
    k = np.array(S["step_sweep"]["k"])
    T = np.array(S["step_sweep"]["horizon"][th])
    m = np.isin(k, K_FIT) & np.isfinite(T)
    slope, icpt, se = linfit(k[m], T[m])
    F["step"] = {"k_used": k[m].tolist(), "slope": slope, "intercept": icpt, "se": se,
                 "lambda": 4 * LN2 / slope, "lambda_se": 4 * LN2 / slope ** 2 * se}
    # the same fit at the other thresholds, to show the slope does not depend on it
    F["step_by_threshold"] = {}
    for other in S["step_sweep"]["horizon"]:
        To = np.array(S["step_sweep"]["horizon"][other])
        mo = np.isin(k, K_FIT) & np.isfinite(To)
        so, _, _ = linfit(k[mo], To[mo])
        F["step_by_threshold"][other] = 4 * LN2 / so

    F["precision_by_k"] = {}
    for kfix, d in S["precision_sweep"].items():
        p = np.array(d["p"])
        Tp = np.array(d["horizon"][th])
        mp = np.isin(p, P_FIT) & np.isfinite(Tp)
        sl, ic, sep = linfit(p[mp], Tp[mp])
        F["precision_by_k"][kfix] = {"p_used": p[mp].tolist(), "slope": sl, "intercept": ic,
                                     "se": sep, "lambda": LN2 / sl, "lambda_se": LN2 / sl ** 2 * sep}
    F["precision"] = F["precision_by_k"]["12"]

    F["best"] = {}
    for p, d in S["mixed"].items():
        Tm_ = np.array(d["horizon"], float)
        i = int(np.nanargmax(Tm_))
        F["best"][p] = {"k": d["k"][i], "T": float(Tm_[i])}
    S["fits"] = F


def main():
    import dp
    dt_save = 2.0 ** -dp.SAVE_K
    curves, summary = {}, {"units": "time in sqrt(l/g); distance in l",
                           "seconds_per_unit_for_1m": float(1 / np.sqrt(9.8))}

    # ---- (a) truncation only: step sweep at 237 bits ---------------------
    ref = load("h", runs.P_REF, runs.K_REF, runs.T_MAIN)
    t_main = np.arange(len(ref)) * dt_save
    curves["t_main"] = t_main
    ks = list(range(4, runs.K_REF))
    hz = {th: [] for th in THRESHOLDS}
    for k in ks:
        ld = log10_bob_distance(load("h", runs.P_REF, k, runs.T_MAIN), ref)
        curves[f"h_p237_k{k:02d}"] = ld
        for th in THRESHOLDS:
            hz[th].append(horizon(t_main, ld, th))
    summary["step_sweep"] = {"k": ks, "horizon": {str(th): hz[th] for th in THRESHOLDS}}

    # energy error for a few of them, to contrast with trajectory error
    for k in (6, 8, 10):
        curves[f"dE_p237_k{k:02d}"] = energy_error(load("h", runs.P_REF, k, runs.T_MAIN))
    curves["dE_p053_k10"] = energy_error(load("h", 53, 10, runs.T_MAIN))

    # ---- (b) roundoff only: precision sweep at fixed step -----------------
    summary["precision_sweep"] = {}
    for k in (8, 12, 16):
        refp = load("p", runs.P_REF, k, runs.T_PREC)
        t_prec = np.arange(len(refp)) * dt_save
        curves["t_prec"] = t_prec
        ps = [p for (kind, p, kk, _, _) in runs.jobs() if kind == "p" and kk == k and p != runs.P_REF]
        ps = sorted(set(ps))
        hz = {th: [] for th in THRESHOLDS}
        for p in ps:
            ld = log10_bob_distance(load("p", p, k, runs.T_PREC), refp)
            curves[f"p_p{p:03d}_k{k:02d}"] = ld
            for th in THRESHOLDS:
                hz[th].append(horizon(t_prec, ld, th))
        summary["precision_sweep"][str(k)] = {"p": ps, "horizon": {str(th): hz[th] for th in THRESHOLDS}}

    # ---- (c) both at once: step sweep at working precisions ---------------
    summary["mixed"] = {}
    for p in (24, 32, 40, 53, 64):
        kk = list(range(4, 18))
        hzp = []
        for k in kk:
            ld = log10_bob_distance(load("h", p, k, runs.T_MAIN), ref)
            curves[f"m_p{p:03d}_k{k:02d}"] = ld
            hzp.append(horizon(t_main, ld, MAIN_THRESHOLD))
        summary["mixed"][str(p)] = {"k": kk, "horizon": hzp}

    # ---- (d) direct measurement of the growth of a perturbation -----------
    base = load("h", runs.P_REF, 14, runs.T_MAIN)
    summary["direct"] = {}
    for which in (0, 3):
        ld = log10_state_distance(load("ly", runs.P_REF, 14, runs.T_MAIN, which), base)
        curves[f"ly_c{which}"] = ld
        ln_d = ld * np.log(10.0)
        for lo, hi in ((0, 160), (20, 160), (0, 80), (80, 160)):
            m = (t_main >= lo) & (t_main <= hi) & np.isfinite(ln_d)
            slope, _, se = linfit(t_main[m], ln_d[m])
            summary["direct"][f"c{which}_{lo}_{hi}"] = {"lambda": slope, "se": se}

    # ---- how well the step-sweep curves collapse under h**-4 ----------------
    # Running maxima, as plotted, and steps fine enough for RK4 to be in its
    # asymptotic regime.
    raw = np.array([running_max(curves[f"h_p237_k{k:02d}"]) for k in range(8, 17)])
    stack = raw + 4 * np.log10(2.0) * np.arange(8, 17)[:, None]
    ok = (t_main > 2) & np.all(raw < -3, axis=0)
    summary["collapse"] = {"decades_spanned": float(stack[:, ok].max() - stack[:, ok].min()),
                           "max_spread_dex": float(np.max(np.ptp(stack[:, ok], axis=0))),
                           "median_spread_dex": float(np.median(np.ptp(stack[:, ok], axis=0))),
                           "t_max": float(t_main[ok].max())}
    # observed order of RK4 from the error at t = 10
    kk = np.arange(8, 17)
    sl, _, se = linfit(kk, [curves[f"h_p237_k{k:02d}"][160] / np.log10(2.0) for k in kk])
    summary["rk4_observed_order"] = {"order": -sl, "se": se}

    # energy error of the ordinary double-precision run used in the energy figure
    summary["energy"] = {"max_dE_p053_k10": float(np.max(curves["dE_p053_k10"])),
                         "horizon_p053_k10": summary["mixed"]["53"]["horizon"][summary["mixed"]["53"]["k"].index(10)]}

    # exponent gamma of the roundoff term: at fixed precision and fixed early
    # time, how the error changes between steps 2**-8 and 2**-16
    early = (t_prec > 5) & (t_prec < 40)
    g = []
    for p in (32, 36, 40, 48, 53, 64, 80, 96, 113, 128):
        a = running_max(curves[f"p_p{p:03d}_k08"])[early]
        b = running_max(curves[f"p_p{p:03d}_k16"])[early]
        g.append(float(np.median(b - a) / (8 * np.log10(2.0))))
    summary["gamma"] = {"median": float(np.median(g)), "min": min(g), "max": max(g)}

    fits(summary)
    np.savez_compressed(ROOT / "data" / "analysis.npz", **curves)
    (ROOT / "data" / "summary.json").write_text(json.dumps(summary, indent=1))
    write_numbers(summary)
    print(json.dumps(summary["fits"], indent=1))
    print(json.dumps(summary["direct"], indent=1))
    print(json.dumps(summary["collapse"], indent=1), summary["rk4_observed_order"], summary["energy"], summary["gamma"])


def write_numbers(S):
    """numbers.tex: every number quoted in either manuscript, as a macro."""
    F = S["fits"]
    sec = S["seconds_per_unit_for_1m"]
    L = np.load(ROOT / "data" / "lyapunov.npz")
    fin = L["final"]
    others = fin[1:]
    chaotic = others[others > 0.05]
    O = json.loads((ROOT / "data" / "orders.json").read_text())
    lam = float(fin[0])

    def T_of(kind, key, k):
        d = S[kind] if key is None else S[kind][key]
        return d["horizon"]["0.01"][d["k"].index(k)] if "0.01" in d["horizon"] else d["horizon"][d["k"].index(k)]

    best = F["best"]
    rows = []
    names = {"24": "24 (single)", "32": "32", "40": "40", "53": "53 (double)", "64": "64 (extended)"}
    for p in ("24", "32", "40", "53", "64"):
        b = best[p]
        rows.append(f"{names[p]} & $2^{{-{b['k']}}}$ & {b['T']:.0f} & {b['T'] * sec:.0f} \\\\")
    Tdbl = best["53"]["T"]
    m = {
        "secPerUnit": f"{sec:.3f}",
        "lamBen": f"{lam:.3f}",
        "lamBenEns": f"{np.median(chaotic):.3f}",
        "nEns": f"{len(others)}",
        "nChaotic": f"{len(chaotic)}",
        "nReg": f"{len(others) - len(chaotic)}",
        "lamDir": f"{S['direct']['c0_0_160']['lambda']:.2f}",
        "lamDirEarly": f"{S['direct']['c0_0_80']['lambda']:.2f}",
        "lamDirLate": f"{S['direct']['c0_80_160']['lambda']:.2f}",
        "lamDirOther": f"{S['direct']['c3_0_160']['lambda']:.3f}",
        "slopeStep": f"{F['step']['slope']:.1f}",
        "slopeStepSE": f"{F['step']['se']:.1f}",
        "slopeStepPred": f"{4 * LN2 / lam:.1f}",
        "lamStep": f"{F['step']['lambda']:.2f}",
        "lamStepSE": f"{F['step']['lambda_se']:.2f}",
        "slopePrec": f"{F['precision']['slope']:.2f}",
        "slopePrecSE": f"{F['precision']['se']:.2f}",
        "slopePrecPred": f"{LN2 / lam:.2f}",
        "lamPrec": f"{F['precision']['lambda']:.2f}",
        "lamPrecSE": f"{F['precision']['lambda_se']:.2f}",
        "gammaRound": f"{S['gamma']['median']:.1f}",
        "gammaLo": f"{S['gamma']['min']:.1f}",
        "gammaHi": f"{S['gamma']['max']:.1f}",
        "lamDirWin": f"{S['direct']['c0_20_160']['lambda']:.2f}",
        "lamPrecCoarse": f"{F['precision_by_k']['8']['lambda']:.2f}",
        "lamPrecFine": f"{F['precision_by_k']['16']['lambda']:.2f}",
        "lamExcess": f"{100 * (F['step']['lambda'] / lam - 1):.0f}",
        "TsingleFine": f"{S['mixed']['24']['horizon'][S['mixed']['24']['k'].index(17)]:.0f}",
        "collapseSpan": f"{S['collapse']['decades_spanned']:.0f}",
        "collapse": f"{S['collapse']['max_spread_dex']:.2f}",
        "rkOrder": f"{S['rk4_observed_order']['order']:.3f}",
        "ordEuler": f"{O['euler']['observed_order']:.3f}",
        "ordMid": f"{O['midpoint']['observed_order']:.3f}",
        "TEulerFine": f"{O['euler']['horizon'][-1]:.0f}",
        "TRKcoarse": f"{T_of('step_sweep', None, 4):.0f}",
        "TkEight": f"{T_of('step_sweep', None, 8):.0f}",
        "TkEleven": f"{T_of('step_sweep', None, 11):.0f}",
        "TBestDouble": f"{Tdbl:.0f}",
        "TBestDoubleSec": f"{Tdbl * sec:.0f}",
        "kBestDouble": f"{best['53']['k']}",
        "TBestSingle": f"{best['24']['T']:.0f}",
        "TBestSingleSec": f"{best['24']['T'] * sec:.0f}",
        "kBestSingle": f"{best['24']['k']}",
        "kRatio": f"{2 ** (17 - best['53']['k'])}",
        "TdoubleFine": f"{S['mixed']['53']['horizon'][S['mixed']['53']['k'].index(17)]:.0f}",
        "TdoubleKten": f"{S['energy']['horizon_p053_k10']:.0f}",
        "dEexp": f"{int(np.floor(-np.log10(S['energy']['max_dE_p053_k10'])))}",
        "costFactor": f"{float(f'{np.exp(lam * Tdbl / 4):.1g}'):.0f}",
        "extraBits": f"{lam * Tdbl * (4 + S['gamma']['median']) / (4 * LN2):.0f}",
        "lamSI": f"{lam / sec:.2f}",
        "decadeSI": f"{4 * np.log(10) / (lam / sec):.0f}",
        "bestRows": NL.join(rows),
    }
    lines = ["% Generated by code/analysis.py from data/summary.json. Do not edit by hand."]
    for k, v in m.items():
        lines.append(BS + "newcommand{" + BS + k + "}{" + v + "}")
    (ROOT / "numbers.tex").write_text(NL.join(lines) + NL, encoding="utf-8")



if __name__ == "__main__":
    main()
