"""Assemble the supplementary material for anonymous review, and prove it is anonymous.

AJP reviews without author names. The repository README names both authors, their
college and the course, so it cannot go to referees. This script builds a separate
bundle, supplement_anonymous/, with its own README that identifies nobody, then
scans every file in the bundle, the anonymous manuscript PDFs and the alt text for anything
that would identify the authors: in the text, in PDF metadata, and in the raw
bytes of binary files. It exits with an error, and writes no zip, if it finds
anything.

  python make_supplement.py
"""

from __future__ import annotations

import re
import shutil
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "supplement_anonymous"

# Strings that would identify the authors, their institution, the course the
# project began in, or the repository. Matched without regard to case.
IDENTIFYING = ["pemmasani", "arjun", "luke", "wang", "harvey", "mudd", "hmc.edu", "claremont",
               "apemm", "math 165", "math165", "math~165", "gmail", "github.com"]

CODE = ["dp.py", "runs.py", "lyapunov.py", "orders.py", "shadow_check.py", "lowprec_check.py",
        "stall_grid.py", "analysis.py", "figures.py", "animate.py"]
DATA = ["summary.json", "analysis.npz", "lyapunov.npz", "orders.json", "shadow_check.json",
        "lowprec.json", "stall_grid.json"]
ANONYMOUS_FILES = ["ajp/manuscript.pdf", "ajp/manuscript_long.pdf", "ajp/alt_text.md",
                   "ajp/upload/fig1.pdf", "ajp/upload/fig2.pdf", "ajp/upload/fig3.pdf"]

# The anonymous manuscript cites a paper by P. Wang and J. Li, who are not the
# authors. Those citation strings are removed before the scan, so that the bare
# surname still counts as a leak anywhere else.
CITED_NAMESAKES = ["wang and li", "p. wang", "p. f. wang", "pengfei wang"]

README = """# Supplementary material

Code, data and animations for the manuscript on step size, floating-point
precision and the predictability horizon of a simulated double pendulum.

## Animations

- `animations/step_size.mp4` (also `.gif`): pendulums that start from the same
  state and are integrated with the same 237-bit arithmetic, differing only in
  the step h. The right-hand panel draws the separation of each from the
  reference run as it grows.
- `animations/precision.mp4` (also `.gif`): the same at a fixed step h = 2^-12,
  with the pendulums differing only in the number of significand bits.

## Code

Python 3 with `numpy`, `matplotlib` and `gmpy2` (an interface to the MPFR
library). `imageio-ffmpeg` is needed only to write MP4 files.

| File | Purpose |
|---|---|
| `code/dp.py` | Equations of motion, energy, and fixed-step RK4 at a chosen precision |
| `code/runs.py` | Every simulation used in the paper; caches trajectories in `data/raw/` |
| `code/lyapunov.py` | Largest Lyapunov exponent by the two-trajectory method |
| `code/orders.py` | The same step sweep with forward Euler and the midpoint rule |
| `code/shadow_check.py` | Long-time averages in single and in double precision |
| `code/lowprec_check.py`, `code/stall_grid.py` | Step-by-step record of the low-precision stall |
| `code/analysis.py` | Error curves, horizons and fits; writes every number quoted in the paper |
| `code/figures.py`, `code/animate.py` | Figures and animations |

To reproduce everything:

```
cd code
python runs.py 20        # about 35 minutes on 20 cores; writes 53 MB to data/raw
python lyapunov.py
python orders.py 16
python shadow_check.py
python lowprec_check.py
python stall_grid.py
python analysis.py
python figures.py
python animate.py both
```

`data/raw/` is not included because of its size; `runs.py` regenerates it, and
the arithmetic is correctly rounded, so the regenerated trajectories are
identical bit for bit.

## Data

`data/analysis.npz` holds every curve plotted, as the base-ten logarithm of the
separation sampled every 1/16 time unit. `data/summary.json` holds the horizons
and fitted slopes. `numbers.tex` holds every number quoted in the text, as
written by `analysis.py`, except the two long-time averages of the shadowing
check, which are in `data/shadow_check.json`.

Units: time in sqrt(l/g), length in l, energy in m g l.
"""


def build():
    # The figure files the submission uses, named by the number they carry in it.
    up = ROOT / "ajp" / "upload"
    up.mkdir(exist_ok=True)
    for n, f in enumerate(("fig2_step_sweep", "fig4_horizons", "fig5_turnover"), start=1):
        shutil.copy2(ROOT / "figures" / f"{f}.pdf", up / f"fig{n}.pdf")
    if OUT.exists():
        shutil.rmtree(OUT)
    for sub in ("code", "data", "animations"):
        (OUT / sub).mkdir(parents=True)
    for f in CODE:
        shutil.copy2(ROOT / "code" / f, OUT / "code" / f)
    for f in DATA:
        shutil.copy2(ROOT / "data" / f, OUT / "data" / f)
    for f in sorted((ROOT / "animations").iterdir()):
        shutil.copy2(f, OUT / "animations" / f.name)
    shutil.copy2(ROOT / "numbers.tex", OUT / "numbers.tex")
    (OUT / "README.md").write_text(README, encoding="utf-8")


def hits_in_bytes(blob: bytes):
    low = blob.lower()
    found = []
    for s in IDENTIFYING:
        for enc in ("utf-8", "utf-16-le", "utf-16-be"):
            if s.encode(enc) in low:
                found.append(s)
                break
    return found


def scan():
    problems = []
    targets = [p for p in OUT.rglob("*") if p.is_file()] + [ROOT / p for p in ANONYMOUS_FILES]
    for path in targets:
        rel = path.relative_to(ROOT)
        blob = path.read_bytes()
        if path.suffix == ".npz":
            # zip members are compressed, so look inside them
            with zipfile.ZipFile(path) as z:
                blob = b"".join(z.read(n) for n in z.namelist()) + " ".join(z.namelist()).encode()
        if path.suffix == ".pdf":
            # PDF streams are compressed, so the raw bytes say little; read the
            # metadata and the extracted text instead.
            from pypdf import PdfReader
            reader = PdfReader(str(path))
            meta = " ".join(f"{k}={v}" for k, v in (reader.metadata or {}).items())
            text = " ".join((pg.extract_text() or "") for pg in reader.pages)
            hay = " ".join((meta + " " + text).lower().split())
            for ok in CITED_NAMESAKES:
                hay = hay.replace(ok, " ")
            found = [s for s in IDENTIFYING if s in hay]
        else:
            found = hits_in_bytes(blob)
        for s in found:
            problems.append(f"{rel}: contains '{s}'")
    return targets, problems


if __name__ == "__main__":
    build()
    targets, problems = scan()
    print(f"scanned {len(targets)} files for {len(IDENTIFYING)} identifying strings")
    if problems:
        print("NOT ANONYMOUS:")
        for p in problems:
            print("  " + p)
        sys.exit(1)
    zpath = ROOT / "supplement_anonymous.zip"
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for p in sorted(OUT.rglob("*")):
            if p.is_file():
                z.write(p, p.relative_to(OUT.parent))
    print(f"clean. wrote {zpath.name} ({zpath.stat().st_size / 1e6:.1f} MB)")
