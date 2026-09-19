"""Static figures for both manuscripts.

Sized for a single journal column (3.37 in) or the full width (6.69 in), with
no text smaller than 8 pt, written as vector PDF and as 600 dpi PNG. Colours are
ordered light to dark within a sweep so the figures still read in greyscale, and
line styles alternate for the same reason.
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import runs

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
COL, FULL = 3.37, 6.69

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "font.family": "serif", "mathtext.fontset": "cm",
    "axes.linewidth": 0.6, "lines.linewidth": 1.0,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def save(fig, name):
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png", dpi=600)
    plt.close(fig)


def shades(n, cmap="viridis", lo=0.0, hi=0.88):
    return [plt.get_cmap(cmap)(x) for x in np.linspace(hi, lo, n)]


def running_max(a):
    return np.maximum.accumulate(np.where(np.isfinite(a), a, -400.0))


def load_traj(kind, p, k, t_end, which=None):
    with open(ROOT / "data" / "raw" / f"{runs.tag((kind, p, k, t_end, which))}.pkl", "rb") as fh:
        tr = pickle.load(fh)["traj"]
    return np.array([[float(v) for v in y] for y in tr])


def bob_xy(y):
    x1, y1 = np.sin(y[:, 0]), -np.cos(y[:, 0])
    return x1, y1, x1 + np.sin(y[:, 1]), y1 - np.cos(y[:, 1])


def main():
    A = np.load(ROOT / "data" / "analysis.npz")
    S = json.loads((ROOT / "data" / "summary.json").read_text())
    L = np.load(ROOT / "data" / "lyapunov.npz")
    t = A["t_main"]
    tp = A["t_prec"]
    ln10 = np.log(10.0)

    # ------------------------------------------------------------------ 1
    # Snapshots: three pendulums that differ only in step size.
    ref = load_traj("h", runs.P_REF, runs.K_REF, runs.T_MAIN)
    show = [(8, "#d95f02", "--"), (11, "#1b9e77", "-.")]
    trajs = {k: load_traj("h", runs.P_REF, k, runs.T_MAIN) for k, _, _ in show}
    T8 = S["step_sweep"]["horizon"]["0.01"][S["step_sweep"]["k"].index(8)]
    T11 = S["step_sweep"]["horizon"]["0.01"][S["step_sweep"]["k"].index(11)]
    times = [0.5 * T8, T8 + 6.0, T11 + 6.0]
    fig, axes = plt.subplots(1, 3, figsize=(FULL, 2.35), sharey=True)
    for ax, tt in zip(axes, times):
        i = int(round(tt * 16))
        for k, colr, ls in show:
            x1, y1, x2, y2 = bob_xy(trajs[k][i:i + 1])
            ax.plot([0, x1[0], x2[0]], [0, y1[0], y2[0]], ls, color=colr, lw=1.3,
                    marker="o", ms=3.5, label=rf"$h=2^{{-{k}}}$")
        x1, y1, x2, y2 = bob_xy(ref[i:i + 1])
        ax.plot([0, x1[0], x2[0]], [0, y1[0], y2[0]], "-", color="k", lw=1.3,
                marker="o", ms=3.5, label=rf"$h=2^{{-{runs.K_REF}}}$")
        ax.set_aspect("equal")
        ax.set_xlim(-2.15, 2.15)
        ax.set_ylim(-2.15, 2.15)
        ax.set_title(rf"$t={tt:.1f}$", fontsize=8, pad=3)
        ax.set_xlabel(r"$x/l$")
    axes[0].set_ylabel(r"$y/l$")
    axes[0].legend(loc="upper left", frameon=False, handlelength=2.2)
    save(fig, "fig1_snapshots")

    # ------------------------------------------------------------------ 2
    # Error curves for the step sweep, and their collapse under h**-4.
    ks = [6, 8, 10, 12, 14, 16]
    cols = shades(len(ks))
    fig, (a, b) = plt.subplots(2, 1, figsize=(COL, 4.3), sharex=True)
    for j, (k, c) in enumerate(zip(ks, cols)):
        ld = A[f"h_p237_k{k:02d}"]
        ls = "-" if j % 2 == 0 else "--"
        a.plot(t, running_max(ld), ls, color=c, label=rf"$2^{{-{k}}}$")
        b.plot(t, running_max(ld) + 4 * k * np.log10(2.0), ls, color=c)
    a.axhline(np.log10(4.0), color="0.5", lw=0.6, ls=":")
    a.set_ylim(-24, 1.5)
    a.set_ylabel(r"$\log_{10}\,(d/l)$")
    a.legend(title=r"step $h$", ncol=2, frameon=False, loc="lower right",
             columnspacing=1.0, handlelength=1.8)
    lam = S["direct"]["c0_0_160"]["lambda"]
    gl = A["ly_c0"]
    off = np.nanmedian((running_max(A["h_p237_k12"]) + 48 * np.log10(2.0) - gl)[(t > 20) & (t < 100)])
    b.plot(t, gl + off, color="k", lw=0.6, label="growth of a perturbation")
    b.set_ylim(-3, 24)
    b.set_xlim(0, 160)
    b.set_ylabel(r"$\log_{10}\,(d/l\,h^{4})$")
    b.set_xlabel(r"time $t\,\sqrt{g/l}$")
    b.legend(frameon=False, loc="lower right")
    a.text(0.04, 0.86, "(a)", transform=a.transAxes)
    b.text(0.04, 0.86, "(b)", transform=b.transAxes)
    fig.subplots_adjust(hspace=0.06)
    save(fig, "fig2_step_sweep")

    # ------------------------------------------------------------------ 3
    # Error curves for the precision sweep at fixed step.
    ps = [16, 24, 32, 40, 53, 64, 80, 96]
    cols = shades(len(ps), "plasma")
    fig, ax = plt.subplots(figsize=(COL, 2.5))
    for j, (p, c) in enumerate(zip(ps, cols)):
        ax.plot(tp, running_max(A[f"p_p{p:03d}_k12"]), "-" if j % 2 == 0 else "--",
                color=c, label=str(p))
    ax.axhline(np.log10(4.0), color="0.5", lw=0.6, ls=":")
    ax.set_xlim(0, 240)
    ax.set_ylim(-32, 1.5)
    ax.set_xlabel(r"time $t\,\sqrt{g/l}$")
    ax.set_ylabel(r"$\log_{10}\,(d/l)$")
    ax.legend(title="significand bits $p$", ncol=4, frameon=False, loc="lower right",
              columnspacing=0.9, handlelength=1.6)
    save(fig, "fig3_precision_sweep")

    # ------------------------------------------------------------------ 4
    # Horizons against log2(1/h) and against p, with straight-line fits.
    F = S["fits"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 2.5))
    k_all = np.array(S["step_sweep"]["k"])
    T_all = np.array(S["step_sweep"]["horizon"]["0.01"])
    # What the horizons should be if the error is A h^4 G(t), with G(t) the
    # measured growth of a perturbation and A fixed once from the collapse.
    G = running_max(A["ly_c0"])
    unsat = (t > 5) & (running_max(A["h_p237_k06"]) < -3)
    offs = np.median([np.median((running_max(A[f"h_p237_k{k:02d}"]) + 4 * k * np.log10(2.0) - G)[unsat])
                      for k in range(6, 17)])
    kc = np.linspace(4, 17.4, 500)
    Tpred = [t[np.argmax(G + offs - 4 * x * np.log10(2.0) >= -2.0)] for x in kc]
    a.plot(kc, Tpred, "-", color="0.6", lw=0.8, zorder=0)
    a.plot(k_all, T_all, "o", ms=3.5, mfc="w", color="k")
    kk = np.array(F["step"]["k_used"])
    a.plot(kk, F["step"]["slope"] * kk + F["step"]["intercept"], "-", color="#d95f02")
    a.set_xlabel(r"$k=\log_2(1/h)$")
    a.set_xticks(range(4, 19, 2))
    a.set_ylabel(r"horizon $T\,\sqrt{g/l}$")
    a.text(0.05, 0.88, "(a) 237-bit arithmetic", transform=a.transAxes)
    for kfix, mk, c in ((8, "s", "#7570b3"), (12, "o", "k"), (16, "^", "#1b9e77")):
        d = S["precision_sweep"][str(kfix)]
        b.plot(d["p"], d["horizon"]["0.01"], mk, ms=3.5, mfc="w", color=c,
               label=rf"$h=2^{{-{kfix}}}$")
    Gp = G[: len(t)]
    unsat_p = (t > 5) & (running_max(A["p_p024_k12"])[: len(t)] < -3)
    offp = np.median([np.median((running_max(A[f"p_p{q:03d}_k12"])[: len(t)] + q * np.log10(2.0) - Gp)[unsat_p])
                      for q in (24, 32, 40, 53, 64, 80)])
    pc = np.linspace(20, 82, 500)
    Tp_pred = np.array([t[np.argmax(Gp + offp - x * np.log10(2.0) >= -2.0)] for x in pc])
    good = Tp_pred > 0
    b.plot(pc[good], Tp_pred[good], "-", color="0.6", lw=0.8, zorder=0)
    pp = np.array(F["precision"]["p_used"])
    b.plot(pp, F["precision"]["slope"] * pp + F["precision"]["intercept"], "-", color="#d95f02")
    b.set_xlabel(r"significand bits $p$")
    b.set_ylabel(r"horizon $T\,\sqrt{g/l}$")
    b.legend(frameon=False, loc="lower right")
    b.text(0.05, 0.88, "(b) fixed step", transform=b.transAxes)
    fig.subplots_adjust(wspace=0.28)
    save(fig, "fig4_horizons")

    # ------------------------------------------------------------------ 5
    # Horizon against step at working precisions: the turnover.
    fig, ax = plt.subplots(figsize=(COL, 2.7))
    ax.plot(k_all, T_all, "-", color="k", lw=0.8, label="237")
    style = {24: ("v", "#e7298a"), 32: ("s", "#7570b3"), 40: ("D", "#1b9e77"),
             53: ("o", "#d95f02"), 64: ("^", "#66a61e")}
    for p, (mk, c) in style.items():
        d = S["mixed"][str(p)]
        ax.plot(d["k"], d["horizon"], "-", marker=mk, ms=3.2, mfc="w", color=c, lw=0.8,
                label=str(p))
    ax.set_xlabel(r"$k=\log_2(1/h)$")
    ax.set_xticks(range(4, 19, 2))
    ax.set_ylabel(r"horizon $T\,\sqrt{g/l}$")
    ax.legend(title="bits $p$", frameon=False, loc="upper left", ncol=2,
              columnspacing=1.0, handlelength=1.8)
    save(fig, "fig5_turnover")

    # ------------------------------------------------------------------ 6
    # Energy is conserved long after the trajectory is wrong.
    fig, ax = plt.subplots(figsize=(COL, 2.5))
    k = 10
    ax.plot(t, running_max(A[f"m_p053_k{k:02d}"]), "-", color="k", label=r"trajectory error $d/l$")
    dE = A["dE_p053_k10"]
    ax.plot(t, np.log10(np.maximum(dE, 1e-300)), "-", color="#d95f02", lw=0.7,
            label=r"energy error $|\Delta E|/mgl$")
    ax.set_xlim(0, 160)
    ax.set_ylim(-17, 1.5)
    ax.set_xlabel(r"time $t\,\sqrt{g/l}$")
    ax.set_ylabel(r"$\log_{10}$ error")
    ax.legend(frameon=False, loc="center right")
    save(fig, "fig6_energy")

    # ------------------------------------------------------------------ 7
    # Convergence of the Benettin estimate.
    fig, ax = plt.subplots(figsize=(COL, 2.3))
    run = L["running"]
    ax.plot(L["t"], run[:, 1:], color="0.75", lw=0.4)
    ax.plot(L["t"], run[:, 0], color="k", lw=1.0)
    ax.set_xscale("log")
    ax.set_xlim(10, L["t"][-1])
    ax.set_ylim(0.1, 0.6)
    ax.set_xlabel(r"averaging time $t\,\sqrt{g/l}$")
    ax.set_ylabel(r"$\lambda\,\sqrt{l/g}$")
    save(fig, "fig7_benettin")


if __name__ == "__main__":
    main()
