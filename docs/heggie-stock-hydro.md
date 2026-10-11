# Heggie Figs. 1–4 with unchanged ThreeFluidSim

This is a separate, self-contained hydrostatic comparison based on upstream
`82d09b8`. It does not require the moving-binary PR or the earlier experimental
hydrostatic solver. `three_fluid.cpp`, `three_fluid.hpp`, `evolution.cpp`, and
the moving solver are unchanged. Only Heggie initialization, observation,
runners, plots, tests, and their build/CLI registration are added.

## The three experimental changes that are NOT included

1. **Relative-dispersion binary heating.** Heggie & Aarseth (1992), Eq. (3),
   uses a relative velocity dispersion for binary–single (BS) encounters.
   The earlier experiment added `c4[i,j]*rho[j]/sqrt(U[i]+U[j])`, with
   derivatives in both new energies and the corresponding constant RHS.
   The public hydrostatic solver uses a single-component denominator.
   This PR does not add a switch, modify that stencil, or disguise a copy
   of the modified solver as reproduction code.
2. **Reflecting sphere.** Section 3.1 of the paper specifies a reflecting
   outer sphere. The earlier experiment fixed the last radius, set the last
   conduction row to `U[N-1]-U[N-2]=0`, and retained the wall's entropy-based
   pressure through realignment. The stock solver instead solves the
   one-sided outer hydrostatic equation, allows the last radius to move,
   and keeps the last U fixed in each conduction solve. This PR retains
   those conditions. A distant low-density boundary is not proof that
   the two closures are equivalent, particularly during late expansion.
3. **Configurable Plummer endpoints.** The earlier experiment extended
   `initPlummer` to accept 1e-6 and 1e4. These were numerical choices, not
   prescribed Heggie parameters. This PR uses the canonical initializer's
   built-in 1e-2..1e3 grid. A profile homology scales the whole native grid
   along with the profile; no custom-grid initializer is copied in.

None of these options exists in upstream `ThreeFluidSim` at `82d09b8`.
The earlier experimental version defaulted both new physics switches off,
but it is not the public version used by this PR.

## Physical settings and the stock heating approximation

All variables inside the solver remain dimensionless. The one-dimensional
dispersion obeys `sigma_i^2 = 2 U_i/3`; both fluids have gamma=5/3. Binary
energy is centre-of-mass kinetic energy, not internal binding energy.
We set C=.104, m_b/m_s=2 for the two-component cases, and ln Lambda=10.
The exact wall radius and Fig.4 Coulomb logarithm are not uniquely specified
by the paper. No formation, ejection, disruption or tidal removal is used.

In relative-U notation the desired BS coefficient is
`B=3.807/(16*pi*lnLambda)` and BB coefficient
`D=5*sqrt(3)/(24*pi*lnLambda)` (paper Eqs. 3 and 10).
The initializer sets the stock coefficients to

```
c4[s,b] = (2/3) B/sqrt(2)
c4[b,s] = (1/3) B/sqrt(2)
c4[b,b] = D/sqrt(2)        # zero in the BS-only case
```

Thus the implemented BS terms use `1/sqrt(2 U_s)`, whereas Heggie uses
`1/sqrt(U_s+U_b)`. They agree when U_s=U_b. Their local rate ratio is
`sqrt((U_s+U_b)/(2 U_s))`; at equipartition, U_b=U_s/2, it is sqrt(3)/2,
about 13.4% smaller. No constant coefficient can make them identical for
arbitrary U_s/U_b. BB self-heating has the correct denominator exactly.
The 2:1 BS heat partition is preserved. Conduction and exchange use
`c2_i=2 sqrt(2/3) C m_i/m_s` and `c1_ij=1/sqrt(3*pi)`.
The companion Mathematica script checks the conversions and approximation.

Figs.1–2 are single-effective-fluid models, not segregating binary systems.
Their self-heating is proportional to `3.807 f_star/(8*pi*sqrt(2 U))` times
the total density. The stock solver supports self-heating in slot FB, so
that slot holds the effective fluid; FS/FD are negligible identical clones,
all particle masses are equal, and exchange is off. This is a storage
choice, not a physical binary fraction of unity. The observer records the
dominant component as column 2. Fig.3 uses M_b/M_s=.01 and no heating.
Fig.4 uses f=N_b/(N_s+N_b)=.03 and BS or BS+BB heating.

## Units and optional profile calibration

The baseline is the canonical hydrostatic Plummer initialization. We use
`L=3*pi/16`, `r_plot=L*r_code`, and
`rho_plot=rho_code/(4*pi*L^3*M_initial_code)` so G=M(0)=1 and the baseline
continuum Plummer energy is -1/4. With each run's initial half-mass radius,
`t/trh(0)=t_code/(.18*ln(10)*rh_code(0)^1.5*sqrt(M_initial_code))`.
Paper Eq.15 uses log10 Lambda; the evolution coefficients use ln Lambda.
Core mass/radius and central relaxation times follow paper Eqs.5, 7, 20.

The optional Fig.4 fit varies the singles' Plummer scale a and relative
binary scale b, preserving the continuum total mass and binary fraction.
Canonical `initPlummer` first uses central density ratio `(M_b/M_s)/b^3`
and binary radius b, then the homology multiplies R by a and divides rho,
P, U by a^3, a^4, a. It preserves the discrete hydrostatic equilibrium.
The observer length unit remains fixed; initial energy is NOT reset to
-1/4. Different component scales also imply different initial dispersions.
The fit therefore deliberately changes the paper's equal-profile setup.

Only the approximate solid-curve readings (rho_s,rho_b)=(3.6,.6) at
t/trh(0)=1 are fitted (scan uncertainty roughly 5–10%). The BS-only case
uses the same fitted profile. Later times are predictions. This PR reruns
the fit using the unchanged solver; do not use previous experimental data
as though they were generated by this runner.

## Build, checks, runs and PDF plots

Initialize the repository's Eigen and Boost.PFR submodules first. The
Python scripts require NumPy, SciPy, Matplotlib and PyMuPDF; Mathematica
is needed only for the symbolic checks.

```sh
git submodule update --init --recursive
make -j2 main-strict
make -j2 check
wolframscript -file script/mathematica/heggie_stock_units.wls
wolframscript -file script/mathematica/heggie_profile.wls
python3 script/run_heggie_hydro_suite.py output/heggie_stock2000 --zones 2000 --jobs 3
python3 script/calibrate_heggie_hydro.py output/heggie_stock_calibrated2000
python3 script/calibrate_heggie_hydro.py output/heggie_stock_calibrated2000 --full
python3 script/plot_heggie_hydro.py /path/to/Heggie_1992.pdf output/heggie_stock2000
python3 script/plot_heggie_calibrated.py /path/to/Heggie_1992.pdf \
  output/heggie_stock2000 output/heggie_stock_calibrated2000
```

The user supplies the paper PDF; it is not redistributed. PDFs are the
default output, and the reference panel aspect ratios and limits are
preserved. The suite caps independent jobs and uses one BLAS/OpenMP thread
per simulation. `main-strict` uses O3/native without fast-math.
The optional single-run CLI is

```
main-strict hydro-heggie OUTPUT MODEL [FSTAR] [ZONES] [FINAL_TRH]
                        [U_TOLERANCE] [MAX_STEPS] [A] [B]
```

Model 0 produces a Fig.1/2 case; 1 produces Fig.3; 2/3 produce Fig.4 BS/BS+BB.
All runs use `evolve()` and a read-only observer. Effective evolution,
initialization, observer and profile settings are saved separately. Outputs
must be fresh directories, including when a previous attempt failed. Final
steps may overshoot. Defaults are U tolerance .001, density ceiling 1e9,
initial timestep 1e-7, maximum timestep 1e4 and one million steps. Density
ceiling is a termination label, not evidence of a bounce or a crash.

## Data format and interpretation

`history.dat` contains little-endian doubles, 13 columns per accepted state:
code time, t/trh(0), plotted central singles/effective density, plotted
binary/trace density, positive central potential depth, plotted core radius,
plotted half-mass radius, Mc/M(0), code total mass, singles central relaxation
time, binary central relaxation time, total mechanical code energy, timestep.
The relaxation-time columns are used only for the two-component Fig.3.
Column 11 (zero-based) is NOT a conservation error or a moving-fluid ledger.

`diagnostics.dat` has 8 columns: t/trh(0), central U_s/U_b, the three outer
radii in code units, outer U_s and U_b, mass fraction outside code radius 100.
It exposes the BS approximation and boundary motion. No claim of exactly
conservative heat flux or mass realignment is made. Binary heating draws on
an unresolved binding reservoir, so mechanical energy alone is not conserved.
Summary JSON files distinguish time, density and failed/step-limit stops.

## Recorded 2000-shell results (unchanged solver)

All eight baseline runs exit without exceptions. The .001/.002 single-fluid
cases stop at the density ceiling at 16.7107/17.9365 trh; .005/.01/.02 reach
1001.19/1000.69/1000.15. Segregation stops at 13.2274 and BS-only at 11.5874.
BS+BB reaches 150.0403, with singles' peak 5.4992 at 10.5656, and maximum
mass drift .542%. Its central BS rate/Heggie rate reaches .8689. The outer
radius changes by 7.47%; maximum mass fraction outside radius 100 is .0728%.
These measurements do not establish that the boundary error is zero.

The new fit is a=.7172313183564625, b=1.5692616385707612. At t/trh=1 it
gives densities 3.5999879/.6000001. Initial energy becomes -.339463, not
-1/4. The full fitted BS+BB run reaches 150.216, with peak singles density
14.7736 at 10.7820 and mass drift .708%. At t/trh=10, fitted densities
14.618/15.752 exceed the paper's approximate 10.8/9.5 by 35%/66%.
At t/trh=100 they are .4082/.5621 versus .35/.39. Matching the left edge
does not reproduce the later evolution exactly. Fitted BS-only still
has no bounce and reaches the density ceiling at 12.0356.

A repeat with 1000 shells and the same fitted parameters (no refit) changes
central singles/binary densities by at most 1.53%/1.31% on 1000 common
logarithmic times from 1 to 150 trh. This is smaller than the remaining
paper mismatch, but is not proof of full convergence or exact conservation.
The baseline and fitted overlays, summaries and density comparison table
are generated by the commands above. The standalone TeX note explains each
curve in Figs.1–4; compile twice from docs with
`pdflatex -output-directory=../output heggie-stock-hydro.tex`.
