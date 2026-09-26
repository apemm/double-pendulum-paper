# How long is a simulated double pendulum right?

A student computational project on how the step size and the floating-point precision limit the
time for which a chaotic trajectory is accurate, and how the Lyapunov exponent sets that time.

Arjun Pemmasani and Luke Wang, Harvey Mudd College. This started as our Math 165 final project in
fall 2025.

## Manuscript

Target: European Journal of Physics, Paper (Mechanics). The manuscript and cover letter are in
ejp/, as TeX and PDF. The manuscript has six figures and two tables, and TeXcount gives 3603 words
including headings and captions (3803 with math units), not counting references. Both authors have
approved submission, and it has not been uploaded yet. ejp/revision_notes.md lists the editorial
changes and numerical checks.

## Results

Time is in units of sqrt(l/g), which is 0.319 s for rods 1 m long.

- With RK4 in 237-bit arithmetic the error at a fixed time scales as h^4.000, and the error curves
  for all steps collapse onto one curve when divided by h^4. That curve is the growth of a
  perturbation of the initial condition, which we measured separately by perturbing the start
  by 2^-170.
- Halving the step extends the accurate time by 4 ln2 / lambda, and each extra bit extends it by
  ln2 / lambda. The slopes give lambda = 0.35 +/- 0.02 and 0.35 +/- 0.01. These equal the growth
  rate measured directly over the same interval and are 18 percent above the long-time exponent
  0.299 from the two-trajectory method, since the first 160 time units of this trajectory are more
  unstable than average.
- In double precision no step keeps the pendulum accurate past about 101 time units (32 s for 1 m
  rods). The best step is 2^-11; 2^-17 costs 64 times as much and gives 98. In single precision the
  best is 53, at 2^-7, and it falls to 30 by 2^-17.
- With 11-bit arithmetic and h = 2^-12 the pendulum never moves. The velocity falls by exactly h per
  step for 2048 steps and reaches -1/2 at t = 1/2, where h is half a unit in its last place, so
  round-to-even freezes it. The increment to the angle is then 2^-13 forever, a quarter of the
  2^-11 needed to change the angle. With h = 2^-8 the angle first moves at t = 0.129. On a grid of
  24 cases every h < 2^-p stalled and every h > 2^-p moved, and at h = 2^-p the increment is a tie
  decided by the last bit of pi/2 (`code/lowprec_check.py`, `code/stall_grid.py`).
- The cost of a longer horizon depends on which exponent is used. Doubling the double-precision
  horizon takes about 2000 times the steps and 50 more bits at the long-time 0.299, or 7000 and 59
  at the finite-time 0.35.
- Energy is conserved to 4e-11 in a double-precision run whose trajectory is lost at t = 92.
- The original observation, that a step 10 times smaller bought only 5 to 10 more seconds, is
  4 ln10 / lambda = 10 s.

## Prior work

A literature search on 19 September 2026 found the main results already in the research
literature for other systems, and the manuscript cites them:

- Li, Zeng and Chou (2001): a given precision has a best step and a longest reliable time.
- Kehlet and Logg (2017): an error estimate of the form used here, with accumulated roundoff
  proportional to h^-1/2. Our fitted exponent is 0.6.
- Wang and Li (2014, arXiv): reliable time is linear in the number of digits with slope
  ln B / lambda, which is the ln2 / lambda per bit law.
- Mendes and Nepomuceno (2016): the Lyapunov exponent from the divergence of two computations that
  round differently.
- Hayes (2004) and Faux and Godolphin (2021), both in AJP: shadowing, and floating-point pitfalls
  for students.
- Calvao and Penna (2015) in EJP and Rafat, Wheatland and Bedding (2009) in AJP: numerical studies
  of the double pendulum.
- Wild (2019): an undergraduate thesis at James Madison University that varies the word size in
  double pendulum simulations.

What is ours is the classroom treatment: making the precision a program parameter so the two errors
can be separated and measured, the collapse onto a perturbed trajectory, the stall, and doing all
of it on the double pendulum.

Wild's thesis uses the same ingredients (RK4, MPFR, the IEEE significand widths 11, 24, 53, 64 and
113, and comparison with the 113-bit run at the same step) and saw 3 of the same things: wider
formats stay accurate longer, the 11-bit solution is constant at h = 1e-4 (put down to underflow
there; the cause is the rounding stall above), and replacing h = 1e-4 by 2^-13 keeps the 53- and
64-bit runs accurate longer. Its results are qualitative, judged by eye from the l2 norm of the
state, and it lists relating the divergence time to the Lyapunov exponent as future work. The
manuscript cites the thesis at each of these 3 points.

We have requested the full Calvao and Penna paper but have not read it yet. The manuscript cites it
as a review of numerical methods and says nothing about an optimum step in it.

## Changes from the original project

- Steps are powers of two, not powers of ten. No binary format represents 10^-4 exactly, so runs at
  different precisions were using slightly different steps.
- Natural units (g = l = 1) remove every parameter, so no constant is rounded differently at each
  precision. The original used g = 9.8.
- Precisions follow the IEEE significand widths (24, 53, 64, 113, 237). The original list
  (11, 23, 52, 113, 237) mixed significand and fraction widths.
- In the precision sweep each run is compared with the 237-bit run at the *same* step, so the
  truncation error cancels exactly and only roundoff is left.
- Python with gmpy2 instead of Julia's BigFloat. Both wrap MPFR. The initial condition, the
  integrator and the 237-bit reference are unchanged.

## Reproducing everything

Needs Python 3 with numpy, matplotlib and gmpy2 (imageio-ffmpeg is optional, for MP4 output), and a
TeX distribution with lmodern, natbib, amsmath, graphicx, booktabs, setspace and hyperref.

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
```

The raw trajectories in data/raw/ are not committed; runs.py regenerates them. The processed data,
figures and animations are committed. The snapshot figure needs the raw trajectories.

Build the manuscript from ejp/:

```sh
pdflatex manuscript.tex
bibtex manuscript
pdflatex manuscript.tex
pdflatex manuscript.tex
pdflatex cover_letter.tex
texcount -inc -sum manuscript.tex
```

Most numbers in the paper are macros in numbers.tex, generated by analysis.py. The paired
horizon-shift comparison is described in ejp/revision_notes.md.

## Layout

- ejp/: manuscript, cover letter and revision notes
- figures/: figures in PDF and PNG
- animations/: step-size and precision animations in GIF and MP4
- code/: integrator, simulation drivers, analysis and plotting
- data/: processed curves, horizons, fitted slopes and check results
- refs.bib: bibliography
- numbers.tex: generated numerical macros

## Before submitting

- Read the final manuscript and cover letter in ejp/.
- Include the code, processed data and animations as supplementary material.
- Check the Acknowledgments and the affiliation details.
- Redo the literature search on submission day (the last one was 19 September 2026).
