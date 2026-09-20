# How long is a simulated double pendulum right?

Two manuscripts from one body of work, on how step size and floating-point
precision limit the time for which a numerical simulation of the chaotic double
pendulum is accurate, and how that time is governed by the largest Lyapunov
exponent.

The work began as the Math 165 (Numerical Analysis) final project of Luke Wang
and I (Arjun Pemmasani) at Harvey Mudd College in fall 2025. This repository is the
version rebuilt for publication.

## The two manuscripts

They are alternatives, not companions. Journals do not allow the same work to be
under review at two of them at once, so the second is submitted only if the
first declines.

| | First choice | Second choice |
|---|---|---|
| Journal | *American Journal of Physics*, Computational Physics section | *European Journal of Physics*, Paper (Mechanics) |
| Folder | `ajp/` | `ejp/` |
| Framing | A computational exercise for instructors, with four suggested problems (six in the long version) | A student computational project, with its educational use and level stated |
| Format | REVTeX 4.2 preprint, 16 pages. See "Length" below | 12 pt article, about 2900 words of body text against a limit of 4000 |
| Review | Anonymous: `manuscript.pdf` carries no names; `manuscript_named.pdf` does | Single-anonymous: names on the manuscript |
| Supplement | `supplement_anonymous.zip`, built and scanned by `code/make_supplement.py`. Never this README, which names the authors | The repository contents |
| Cover letter | `ajp/cover_letter.pdf` | `ejp/cover_letter.pdf`, which EJP requires, stating level and usefulness |
| Also needed | `ajp/alt_text.md`, `ajp/statement_for_submission_form.md`, `ajp/upload/fig1.pdf` to `fig3.pdf` | |

## Length of the AJP manuscript

AJP asks that papers normally fit on 6 journal pages. Two estimates disagree, and
the truth for AJP's own layout is probably between them.

| Build | Preprint pages | AJP's rule of thumb (pages / 3) | Two-column REVTeX build |
|---|---|---|---|
| `manuscript.pdf`, the one to submit | 16 | 5.3 | about 6.2, the last 0.2 being references |
| `manuscript_long.pdf` | 19 | 6.3 | about 7.0 |

The submitted version leaves out the snapshot figure, the precision-sweep curves
(their content is in the horizon figure and the second animation), the energy
figure (its numbers stay in the text), suggested problems 4 and 5, and the
general form of the stall rule, and runs the list of five lessons into a
paragraph. All of it is in the source behind `\ifdefined\longversion`, so
compiling `manuscript_long.tex` restores it if an editor allows the space.
`manuscript_twocolumn.tex` exists only to measure length.

## What was found

All in units where time is measured in sqrt(l/g); one unit is 0.319 s for rods
one metre long.

- With RK4 in 237-bit arithmetic the error at fixed time scales as h^4.000, and
  the error curves for all steps collapse onto one curve when divided by h^4.
  That curve is the growth of a perturbation of the initial condition, measured
  separately by perturbing the start by 2^-170.
- Halving the step extends the accurate time by 4 ln2 / lambda; each extra bit of
  precision extends it by ln2 / lambda. The slopes give lambda = 0.35 +/- 0.02 and
  0.35 +/- 0.01, equal to the directly measured growth rate over the same interval
  and 18 percent above the long-time exponent 0.299 from the two-trajectory
  method, because the first 160 time units of this trajectory are more unstable
  than average.
- In double precision no step keeps the pendulum accurate beyond about 101 time
  units (32 s for one-metre rods). The best step is 2^-11; going to 2^-17 costs 64
  times as much and gives 98. In single precision the best is 53, at 2^-7, and it
  falls to 30 by 2^-17.
- With 11-bit arithmetic and h = 2^-12 the pendulum never moves. The velocity falls
  by exactly h per step for 2048 steps, reaches -1/2 at t = 1/2, where h is half a
  unit in its last place, and round-to-even freezes it. The increment to the angle
  is then 2^-13 for ever, a quarter of the 2^-11 needed to change the angle. With
  h = 2^-8 the angle first moves at t = 0.129. On a grid of 24 cases every
  h < 2^-p stalled and every h > 2^-p moved, and at h = 2^-p the increment is a tie
  decided by the last bit of pi/2 (`code/lowprec_check.py`, `code/stall_grid.py`).
- The cost of a longer horizon depends on which exponent is used. Doubling the
  double-precision horizon multiplies the steps by about 2000 and needs about 50
  more bits at the long-time 0.299, or 7000 and 59 at the finite-time 0.35.
- Energy stays conserved to 4e-11 in a double-precision run whose trajectory is
  lost at t = 92.
- The original observation, that a step ten times smaller bought only five to
  ten more seconds, is 4 ln10 / lambda = 10 s.

## Prior work that had to be cited

A literature search on 19 September 2026 (ten searches by different routes,
every result checked against Crossref or arXiv) found that the main results exist
in the research literature for other systems, and both manuscripts now say so.

- Li, Zeng and Chou (2001): a given precision has a best step and a longest
  reliable time.
- Kehlet and Logg (2017): an error estimate of the form used here, with
  accumulated roundoff proportional to h^-1/2. The fitted exponent here is 0.6.
- Wang and Li (2014, arXiv): reliable time is linear in the number of digits with
  slope ln B / lambda, which is the one-bit-buys-ln2/lambda law.
- Mendes and Nepomuceno (2016): the Lyapunov exponent estimated from the
  divergence of two computations that round differently.
- Hayes (2004) and Faux and Godolphin (2021), both in AJP: shadowing, and
  floating-point pitfalls for students.
- Calvao and Penna (2015) in EJP and Rafat, Wheatland and Bedding (2009) in AJP:
  numerical studies of the double pendulum.
- Wild (2019): an undergraduate thesis at James Madison University that varies
  the word size in double pendulum simulations.

What remains ours is the classroom treatment: separating the two errors with the
precision as a program parameter, measuring each, the collapse onto a perturbed
trajectory, the stall, and doing all of it on the double pendulum.

Wild's thesis has now been read in full, and the papers describe it from the
text. The thesis uses the same ingredients (RK4, MPFR, the IEEE significand widths 11,
24, 53, 64 and 113, comparison with the 113-bit run at the same step) and
observed three things found here as well: wider formats stay accurate for
longer, the 11-bit solution is constant at h = 1e-4 (attributed there to
underflow; the mechanism is the rounding stall described above), and replacing
h = 1e-4 by 2^-13 keeps the 53- and 64-bit runs accurate for longer. The results
are qualitative, judged by eye from the l2 norm of the state, and relating the
divergence time to the Lyapunov exponent is listed as future work. Both papers
cite the thesis at each of those three points.

One thing still to do. Calvao and Penna is not available online and has been
requested. What the papers say about it (that a step of 1e-4 served RK4 best on
this system and smaller steps gave less precision) rests on two secondary
sources that agree: Wild's thesis, which took its step size from that paper and
says so, and one literature search that reported the same. Check that sentence
against the paper itself when it arrives, and see whether they judged accuracy by
energy error, which would deserve a sentence of contrast with the energy result
here.

## Changes from the original project, and why

- Steps are powers of two, not powers of ten. No binary format represents 10^-4
  exactly, so runs at different precisions were also using slightly different
  steps.
- Natural units (g = l = 1) remove every parameter, so no constant has to be
  rounded differently at each precision. The original used g = 9.8.
- Precisions follow the IEEE significand widths (24, 53, 64, 113, 237). The
  original list (11, 23, 52, 113, 237) mixed significand and fraction widths.
- In the precision sweep each run is compared with the 237-bit run at the *same*
  step, so truncation error cancels exactly and only roundoff is seen.
- Python with gmpy2 in place of Julia's BigFloat. Both wrap MPFR. The initial
  condition, the integrator and the 237-bit reference precision are unchanged.

## Reproducing everything

Requires Python 3 with `numpy`, `matplotlib`, `gmpy2`, `pypdf` (for
`make_supplement.py`) and optionally `imageio-ffmpeg` for MP4 output; and a TeX
distribution with REVTeX 4.1, which is what AJP's own sample manuscript uses.

```
cd code
python runs.py 20          # all simulations, about 35 minutes on 20 cores, 53 MB in data/raw
python lyapunov.py         # two-trajectory exponent, about 10 minutes
python orders.py 16        # Euler and midpoint check
python shadow_check.py     # long-time averages in single and double precision
python lowprec_check.py    # the 11-bit stall, step by step
python stall_grid.py       # the stall rule on a grid of precisions and steps
python analysis.py         # error curves, horizons, fits -> data/*.json, numbers.tex
python figures.py          # figures/*.pdf and *.png
python animate.py both     # animations/*.gif and *.mp4
python make_supplement.py  # anonymous bundle for AJP review, scanned for names
```

`data/raw/` is not committed; it is regenerated by `runs.py`. Everything
downstream of it is committed, so the three submitted figures, the animations'
inputs and both PDFs can be rebuilt without re-running the simulations. The one
exception is the snapshot figure of the long version, which reads the raw
trajectories and is skipped when they are absent. Build the PDFs with
`bash ajp/build.sh`, which fixes the clock so the PDF metadata carries no time
zone.

Every number quoted in either manuscript is a macro in `numbers.tex`, which
`analysis.py` writes from the data. Neither manuscript contains a number typed by
hand, other than the two long-time averages in the shadowing check.

## Layout

```
ajp/                 manuscript.tex (one source); _named, _long and _twocolumn wrappers; cover letter; alt text
ejp/                 manuscript, cover letter
figures/             fig1 to fig7, PDF and 600 dpi PNG
animations/          step_size and precision, GIF and MP4 (supplementary material)
code/                dp.py (integrator), runs.py, analysis.py, figures.py, animate.py, ...
data/                summary.json, analysis.npz, lyapunov.npz, orders.json, shadow_check.json
refs.bib             shared bibliography; every entry checked against Crossref
numbers.tex          generated; every quoted number
```

## Before submitting

- Authorship: Arjun Pemmasani is first and corresponding author, Luke Wang second.
  Confirm Luke's affiliation line.
- AJP encourages a proposal to the Computational Physics section editors before
  submission, summarizing the physics, the algorithm and the intended level. Send
  it, and then say in the cover letter that it was sent and what the reply was.
- AJP: upload `ajp/manuscript.pdf` (anonymous), `ajp/upload/fig1.pdf` to
  `fig3.pdf`, `ajp/alt_text.md`, `supplement_anonymous.zip`, and paste
  `ajp/statement_for_submission_form.md` into the form. Do not upload this README
  or link this repository, since both name the authors.
- AJP requires disclosure to the editor, in the cover letter, if artificial
  intelligence was used to generate any of the text. That decision is the
  authors'.
- EJP: upload `ejp/manuscript.pdf` and the cover letter.
- Run a fresh literature search on the day. AJP rejects without review for
  missing closely related work, and the search here is current to 19 September
  2026.
