# Statler comparison through the shared evolution driver

This example migrates the existing Statler–Ostriker–Cohn comparison to the
same pattern as the AB examples: parameter setup, initialization, observer
construction, `sim.evolve(observer)`, then observer saving. There is no
study-specific time loop or generated copy of solver source.

## Run

Initialize the repository's pinned dependencies and build:

```sh
git submodule update --init --recursive
make -j3 main-strict
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./main-strict statler direct output/statler/direct
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./main-strict statler control output/statler/control
```

`main-strict` uses the normal source with `-ffast-math` removed. The existing
`main` target also exposes the same example; no arguments retains the
previous binary-formation examples. Optional positional limits are
`max_steps` and `final_time_trh`:

```sh
./main-strict statler direct output/statler/smoke 100 10000
```

The full defaults are 500 zones, two million maximum accepted steps,
`1e4 trh`, and a maximum-species central-density ceiling of `1e8` in code units.
The control normally reaches the density ceiling first. Saved histories
include the initial and terminal accepted states, not rejected trials.
Rerunning in a directory replaces its named result files.

## Configuration and physical scope

`statler_reproduction.cpp/hpp` contains the external initializer and preset.
Configuration has three independent owners:

- `ThreeFluidSim::param` (`ThreeFluidParam`): dimensionless evolution
  settings only, saved in the run directory. There is no copying
  `parameters()` getter or duplicate set of simulator settings.
- `StatlerInitParam`: physical scales and Plummer inputs, saved under
  `initialization/`. `initializeCaptureCluster()` calls the canonical
  `initPlummer` and `initCoeffs`. The Yiming initializers remain available
  for old examples but are deprecated, not alternative initialization modes.
- `StatlerObserver::param` (`StatlerObserverParam`): reference-unit
  conversions and snapshot schedule, saved under `observer/`. The observer
  can be passed alone or composed with other observers through `ObserverPack`.
  No observer-selection flag belongs to the simulator.

Each record uses the same field-name/type/offset serializer. The external
initializer performs physical-unit conversions; physical constants are
`constexpr` there, not adjustable evolution options. `param.Deltat` is the
initial timestep and `sim.Deltat` is the current adaptive timestep. The
former is not overwritten during evolution. Source exponents, capture
partition, two relaxation passes, retry safety and growth bounds are fixed
parts of the implemented method, not runtime switches.

Physical metadata: `N_star=3e5`, stellar mass `0.7 Msun`, stellar radius
`0.57 Rsun`, Plummer scale `1.13 pc`, and reference half-mass relaxation time
`225 Myr`. Binary and DM density seeds are `1e-12` and `1e-20` relative to
singles. Binaries have twice the stellar mass; DM interactions are disabled.
Physical scales are metadata only: all evolving physical state remains
dimensionless, including time and the cumulative formed-binary count.

The Press–Teukolsky capture operator implements

\[
\widehat S=A\widehat\rho_s^2\widehat u_s^{-p},\qquad
\Delta\widehat\rho=\widehat m_b\widehat S\Delta\widehat t,
\qquad
\widehat u_b'=
\frac{\widehat\rho_b\widehat u_b+
\Delta\widehat\rho\,\eta\widehat u_s}{\widehat\rho_b+\Delta\widehat\rho}.
\]

For this example `p=0.6`, `eta=0.5`, and the Press–Teukolsky coefficient is

\[
A=\frac{107\Gamma(0.9)}{50\,2^{0.7}3^{0.4}\pi^{3/2}}
\left(\frac{R_*}{r_0}\right)^{0.9}
\frac{\widehat m_s^{-1.9}}{\ln(0.8N_*)}.
\]

The applied source is tidal capture only. The three-body diagnostic is
recorded as a counterfactual rate, not applied formation. Capture conserves
local mass but removes resolved random energy
`(1-eta) * transferred_density * U_s`; it is not the energy-preserving
legacy formation mode. `direct` retains the existing binary heating
coefficients; `control` sets them to zero. Neither includes binary
disruption, binding-energy evolution or encounter-driven ejection.

Thus `direct` is the closest available comparison to the paper's
direct-heating-only experiment. `control` is a capture/cooling control,
**not** its ejection-only experiment. The atlas labels missing observables
and physics rather than treating numerical mass loss as physical escape.

## Shared numerical machinery

`evolve()` only validates, observes and calls `advanceAcceptedStep()`.
Separate functions prepare capture sources, apply formation, project
hydrostatic equilibrium, compute the conduction change and select the next
timestep. The accepted-step order remains conduction, optional tidal sink,
formation, two corrections per fluid, and realignment.

The bounded controller compares fractional conduction changes to `thres`
and logarithmic central-density changes to `density_change_tolerance`.
Central density for control and stopping is the maximum over species.
The controller clamps timestep growth and caps the step at the remaining integration
interval. The capture donor limit retries conduction at a smaller step
before any density changes; both U and P are restored. This is not a full
post-projection error-control or rollback scheme. Growth factors are fixed
at `[0.5,1.5]`, the error floor is `1e-12`, retry safety is `0.8`, and the
retry limit is 100. These constants are recorded by the source revision.

Intentional shared-driver corrections: stop at time-limit equality, clip
the last step, bound a zero-error timestep proposal, and throw on numerical
failure instead of `exit(0)`. There is now one bounded controller for all
examples, not a legacy/new dispatcher. This deliberately changes timestep
histories for examples that previously used the unbounded U-only controller.
Existing
interpolation-based realignment and the outer fixed-temperature conduction
row are unchanged, including their conservation limitations.

`sim.param.runtime_validation=0` disables automatic full-state finiteness,
positivity and shell-order scans; it defaults to 1. Cheap configuration,
solver-error, donor-depletion and time-advance guards remain mandatory.
`sanityCheck()` can always be invoked explicitly. Use `main-strict` when
relying on IEEE finiteness checks; fast-math can invalidate those checks.

## Observers and plotting

`StatlerObserver` is read-only with respect to the simulator. It computes
the legacy `history_*.dat` and `snapshot_*.dat` diagnostics, preserves the
atlas's unit conventions, and estimates mean binary age from accepted
formation increments. That age estimate assumes numerical remapping loss
is age-independent; it is not an evolved age distribution. The snapshot
schedule belongs to the observer parameter record. Every accepted state is
recorded; the peak snapshot is always retained and output is binary. These
are fixed conventions, not extra options. History and snapshot arrays are
little-endian float64; the zone-count file is int64. Existing offset-aware
Python readers handle each independent parameter record.

Supply your own copy of the reference paper; it is not bundled:

```sh
python3 script/plot_statler_comparison_atlas.py Statler_Ostriker_Cohn.pdf \
  output/statler/control output/statler/direct \
  --output output/statler/comparison.pdf
```

The Python analysis uses NumPy, Matplotlib and PyMuPDF. The generated atlas
contains one page per numbered figure (23 figures). It does not imply that
the fluid closure reproduces every plotted physical quantity.

## Validation

`make check` runs the hydrostatic test, `check_statler`, bounded checks of
the eight existing examples, and Python tests.
The new checks cover local formation mass/energy balance, legacy 100-step
direct/control values, observer lengths and saved parameter readback,
time-limit equality and clipping, rejected-trial isolation/restoration,
zero-error timestep control, observer composition, and runtime validation
on/off. The retry test compares the accepted state bitwise with a fresh
run using the reduced initial timestep.

Historical validation of commit `683fdfb`, before the parameter simplification:
the migrated 500-zone direct run reaches `1e4 trh` in 12,141 accepted steps;
the control stops at its density ceiling after 11,237. All 16 full-history
diagnostics agree with the saved legacy calculation to a maximum difference
below `3e-10`, normalized by each diagnostic's maximum absolute legacy
value. The direct density maximum is at `14.2753387104 trh`. These checks
establish preservation of the existing comparison, not improved physical
agreement with the paper or improved global energy conservation.

At that revision, the eight existing sample configurations matched the preceding public
revision bit for bit over 100 steps when built with identical compiler
flags. Do not mix Eigen objects compiled with different architecture or
alignment flags in one executable; the comparison builds used the same
`-O3 -march=native -DNDEBUG` strict-IEEE settings throughout.

The current checks retain the 100-step Statler reference values (relative
tolerance `2e-10`). The previous full-history and eight-example bitwise
claims above do not extend to the intentionally unified timestep controller.
Both command-line Statler configurations were also run for 1,000 steps after
the simplification: all 16 history arrays match the corresponding prefixes
of the previous full runs bitwise.
