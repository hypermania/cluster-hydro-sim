# Moving-fluid spatial resolution study

This study changes only the radial cell count, retaining the initial radial
domain, physical coefficients, timestep tolerance, boundary convention and
stopping criteria. The finite-volume initialization adjusts particle masses
and the capture coefficient very slightly to retain the prescribed stellar
number; the 500/1000-cell coefficient difference is about 1e-7 relative.
This is not a timestep-convergence test: adaptive accepted timesteps may
differ when the grid changes.

## Existing Statler oscillations

The plotting code reads the central dispersion history directly, without
smoothing. Over 50–10000 initial half-mass relaxation times, the 500-cell
history has three rebounds with log prominence exceeding log(1.03), at
approximately 173, 539 and 746 trh. The corresponding fractional prominences
are .073, .815 and .589. The 1000-cell history is monotonically decreasing
over this interval. These are features of the saved solution, not different
profile snapshot schedules or plot interpolation settings.

This strongly supports a resolution-dependent numerical origin, but does
not isolate an operator or prove that the 1000-cell solution is converged.
The diagnostic uses 10000 common logarithmic-time samples and no smoothing;
peak prominence is not a physical oscillation-amplitude definition.

One concrete grid-dependent candidate is the frozen initial reference-force
correction. At r=0.1 Plummer scales its ratio to initial self-gravity is
-.0944, -.0462, -.0228 for 500, 1000, 2000 cells. At r=10 the ratios are
.0842, .0436, .0222. These numbers follow directly from the solver's initial
pressure/source quadrature, reconstructed in the analysis script. As gravity
weakens during expansion, a fixed correction can become relatively more
important. This is a hypothesis to test independently, not a demonstrated
explanation of all observed features. No correction or transport operator
is changed in this study.

## 2000-cell runs

The eight Heggie runs at 2000 cells completed. The .001/.002 and segregation
models hit the density ceiling at 16.691, 17.865 and 13.476 initial
relaxation times; .005, .01, .02 and BS+BB reached their requested endpoints;
BS-only hit the ceiling at 13.744. Mass drift remains below 3e-15 and the
step energy-ledger residual below 4e-17. The PDFs are in
`output/heggie2000/plots/`; the 400-to-2000 BS+BB comparison reports relative
differences of about .92 (single density), .97 (binary density), 2.17 (core
radius), and .80 (core-mass fraction). Thus increasing resolution does not
yet establish convergence.

A 2000-cell Statler run reached approximately 1728.160 trh before the long
process was stopped. The runner serializes the full history only when
`evolve()` returns, so this interrupted run has a budget log but no
reconstructable profile/history files and cannot honestly be plotted. A
second attempt was stopped during its slow initial phase. The complete
Statler 2000-cell overlay remains pending; the partial run is retained as a
timing record, not used as scientific output.
The endpoint is recovered from the last complete budget row; the earlier
1309-trh report used the less frequent progress messages. Its final recorded
timestep is .04633 code units (.1125 trh), so the timestep is allowed to
increase. Small early timesteps, including approximately 5.29e-13 code
units when a trace binary density changes rapidly, contribute to runtime.

Six bounded startup timing trials per build at 2000 cells measured median
step times of 24.96 ms with O3/native, 25.19 ms with O3/LTO/no-math-errno,
and 25.70 ms with Ofast/LTO. These overlapping timings show no reliable
aggressive-optimization speedup; they are not late-time benchmarks.

## Commands

Use fresh output directories and the strict-floating-point executable:

```sh
make -j2 main-strict
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./main-strict moving-statler output/moving_statler_2000 2000 10000 2000000
python3 script/run_heggie_suite.py output/heggie2000 --zones 2000 --jobs 3
python3 script/plot_moving_statler.py path/to/Statler_Ostriker_Cohn.pdf output/moving_statler_2000 --compare output/moving_statler_1000
python3 script/analyze_statler_resolution.py output/moving_statler_matched500 output/moving_statler_1000 --output output/moving_resolution_500_1000
python3 script/plot_heggie_comparison.py path/to/Heggie_1992.pdf output/heggie2000
```

The Heggie suite includes all eight cases used for Figs. 1–4. Density-ceiling
stops remain labelled; a successful density-ceiling stop does not mean the
reference post-collapse solution was reproduced. The Statler run retains
the approximate reference snapshot epochs used by the matched500 run.
Generated PDFs/data are not committed.

The 2000-cell jobs were launched from code revision `f65379d`; the subsequent
audit used `2e49bce`. See [the Fig. 4 density audit](heggie-reproduction.md#2000-cell-fig-4-density-audit)
for the normalization check and remaining discrepancy. The initial 500/1000
comparison is not a substitute for a completed 2000-cell Statler run.
