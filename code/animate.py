"""Supplementary animations.

Each animation overlays several pendulums that start from the same state and
obey the same equations, and differ only in how the arithmetic was done. The
left panel shows the pendulums with a short trail on the lower bob. The right
panel shows the separation of each from the reference on a logarithmic scale,
drawn up to the current time, so the viewer can watch the error climb for tens
of time units while the pendulums on the left still look identical, and then
see them peel away one at a time in the order the right panel predicts.

  python animate.py step        pendulums differing in step size, 237 bits
  python animate.py precision   pendulums differing in precision, h = 2**-12
  python animate.py both

Writes GIF always and MP4 when an ffmpeg binary is available (through the
imageio-ffmpeg package if it is installed).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter

import runs
from figures import load_traj, bob_xy, running_max

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "animations"

FPS = 25
STRIDE = 2          # samples are 1/16 apart; two per frame = 1/8 time unit
TRAIL = 40          # frames of trail on the lower bob

try:
    import imageio_ffmpeg
    matplotlib.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    HAVE_FFMPEG = True
except Exception:
    HAVE_FFMPEG = False


def build(name, members, reference, t, curves, t_end, title):
    """members: list of (label, colour, trajectory); reference drawn last, in black."""
    plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm"})
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.0, 4.4), dpi=100,
                                 gridspec_kw={"width_ratios": [1.0, 1.15]})
    ax.set_aspect("equal")
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-2.2, 2.2)
    ax.set_xticks([])
    ax.set_yticks([])
    clock = ax.text(0.03, 0.95, "", transform=ax.transAxes, va="top")
    ax.plot([0], [0], marker="+", color="0.35", ms=9, mew=1.2, zorder=5)   # the pivot
    fig.suptitle(title, fontsize=10, y=0.975)

    bx.set_xlim(0, t_end)
    bx.set_ylim(-26, 1.5)
    bx.axhline(np.log10(4.0), color="0.6", lw=0.6, ls=":")
    bx.set_xlabel(r"time $t\sqrt{g/l}$")
    bx.set_ylabel(r"$\log_{10}$ (separation of lower bobs / $l$)")

    everyone = members + [reference]
    xy = [bob_xy(tr) for _, _, tr in everyone]
    rods, trails, errs = [], [], []
    for label, colr, _ in everyone:
        lw = 2.0 if colr == "k" else 1.4
        (rod,) = ax.plot([], [], "-o", color=colr, lw=lw, ms=5, label=label)
        (trl,) = ax.plot([], [], "-", color=colr, lw=0.6, alpha=0.6)
        rods.append(rod)
        trails.append(trl)
    for (label, colr, _), c in zip(members, curves):
        (e,) = bx.plot([], [], "-", color=colr, lw=1.2, label=label)
        errs.append((e, running_max(c)))
    (cursor,) = bx.plot([], [], "-", color="0.4", lw=0.6)
    ax.legend(loc="lower left", frameon=False, fontsize=8)
    bx.legend(loc="lower right", frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    n_frames = int(t_end * 16) // STRIDE

    def draw(f):
        i = f * STRIDE
        for (x1, y1, x2, y2), rod, trl in zip(xy, rods, trails):
            rod.set_data([0, x1[i], x2[i]], [0, y1[i], y2[i]])
            j = max(0, i - TRAIL * STRIDE)
            trl.set_data(x2[j:i + 1], y2[j:i + 1])
        for e, c in errs:
            e.set_data(t[:i + 1], c[:i + 1])
        cursor.set_data([t[i], t[i]], [-26, 1.5])
        clock.set_text(rf"$t\sqrt{{g/l}}={t[i]:6.1f}$")
        return rods + trails + [e for e, _ in errs] + [cursor, clock]

    anim = FuncAnimation(fig, draw, frames=n_frames, blit=True)
    OUT.mkdir(exist_ok=True)
    anim.save(OUT / f"{name}.gif", writer=PillowWriter(fps=FPS))
    if HAVE_FFMPEG:
        anim.save(OUT / f"{name}.mp4", writer=FFMpegWriter(fps=FPS, bitrate=2400))
    plt.close(fig)
    print("wrote", name, n_frames, "frames", "(gif+mp4)" if HAVE_FFMPEG else "(gif)")


def step():
    A = np.load(ROOT / "data" / "analysis.npz")
    ks = [(6, "#e7298a"), (8, "#d95f02"), (10, "#7570b3"), (12, "#1b9e77")]
    members = [(rf"$h=2^{{-{k}}}$", c, load_traj("h", runs.P_REF, k, runs.T_MAIN)) for k, c in ks]
    ref = (rf"$h=2^{{-{runs.K_REF}}}$ (reference)", "k",
           load_traj("h", runs.P_REF, runs.K_REF, runs.T_MAIN))
    curves = [A[f"h_p237_k{k:02d}"] for k, _ in ks]
    build("step_size", members, ref, A["t_main"], curves, 140,
          "Same pendulum, same 237-bit arithmetic, different step $h$")


def precision():
    A = np.load(ROOT / "data" / "analysis.npz")
    ps = [(16, "#e7298a"), (24, "#d95f02"), (32, "#7570b3"), (53, "#1b9e77")]
    members = [(rf"$p={p}$ bits", c, load_traj("p", p, 12, runs.T_PREC)) for p, c in ps]
    ref = (rf"$p={runs.P_REF}$ bits (reference)", "k", load_traj("p", runs.P_REF, 12, runs.T_PREC))
    curves = [A[f"p_p{p:03d}_k12"] for p, _ in ps]
    build("precision", members, ref, A["t_prec"], curves, 160,
          "Same pendulum, same step $h=2^{-12}$, different precision")


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "both"
    if what in ("step", "both"):
        step()
    if what in ("precision", "both"):
        precision()
