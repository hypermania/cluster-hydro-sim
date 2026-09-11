# Cluster Hydro Simulation

`cluster-hydro-sim` is a research code for the secular evolution of a
spherically symmetric, self-gravitating system with three components:

- single stars (`s`),
- stellar binaries (`b`), and
- dark matter (`d`).

The code evolves dimensionless density, enclosed mass, pressure, and specific
kinetic energy profiles on radial Lagrangian grids. Its numerical structure
extends the two-fluid conduction method of Zhong and Shapiro to three fluids
and includes experimental binary-formation, binary-heating, and tidal-removal
terms.

The current executable is a research prototype rather than a general-purpose
simulation package. In particular, experiment selection and parameters are
set in `src/main.cpp`; the Statler example additionally has a command-line
selector. Tests cover hydrostatic projection, the shared Statler evolution
and output readers. The mode-2
binary-formation example has also been audited operator by operator; formation
is conservative to roundoff, but the complete trajectory is currently
dominated by nonconservative grid realignment.

## Numerical method

For each fluid `f`, the principal arrays are

```text
R[f][i]     radial shell boundary
Rho[f][i]   shell-averaged mass density
Menc[f][i]  enclosed mass
P[f][i]     pressure
U[f][i]     specific kinetic energy
```

The code uses

```text
P = (2/3) Rho U,        U = (3/2) P/Rho.
```

With fiducial density `rho0` and radius `r0`, the mass unit is

```text
M0 = 4 pi rho0 r0^3.
```

The factor `4 pi` is absorbed into this unit, so the dimensionless enclosed
mass is `Menc = integral Rho r^2 dr`. The particle-mass fields `ms`, `mb`, and
`md` are stored in units of `M0`; coefficient initialization separately takes
the ratios `Mtot/ms`, `mb/ms`, and `md/ms`.

One evolution step consists of:

1. A semi-implicit conduction and inter-fluid energy-exchange solve. The
   resulting block-banded system is solved with `LAPACKE_dgbsv`. This solve
   also includes binary-heating source terms.
2. Optional tidal mass removal.
3. Optional binary formation. Mode 2 transfers density locally from single
   stars to binaries and chooses the updated binary specific energy to
   preserve `Rho_s U_s + Rho_b U_b` pointwise.
4. Two linearized hydrostatic-relaxation solves per component, using
   `LAPACKE_dgtsv` while preserving shell mass and specific entropy within
   each relaxation solve. The solve includes the finite-mass outer shell and
   imposes the one-sided equation
   `-P[N-1]/(R[N-1]-R[N-2]) + M[N-1] Rho[N-1]/R[N-1]^2 = 0`.
5. Logarithmic-grid realignment, followed by reconstruction of enclosed mass
   and hydrostatic pressure.

Formation is therefore applied before hydrostatic projection. The accepted
state at the end of the step has passed through relaxation and the common-grid
pressure reconstruction.

The center uses a regularity boundary condition, `U[f][0] = U[f][1]`, and the
outer value of `U` is held fixed during the conduction solve. Entry `N-1` is
an active shell with positive density and pressure; the adjacent conduction
stencil uses its density, while the fixed outer-`U` row acts as a prescribed
thermal boundary. The timestep starts at `1e-3` and is adjusted using a target
maximum relative change in `U` of `1e-3`. The measured change determines the
next step through `dt_next = dt * u_change_tolerance / maxChange`, capped
by `max_timestep` (also used when the measured change is zero). Density and
formation changes are not included in this estimate; there are no retries
or multiplicative growth clamps. Invalid capture transfers terminate the run.
Evolution stops at the configured density, step or time limit;
the final step is taken in full and may overshoot `maxTime`. Observers save
the actual time, which both Python plotters use without endpoint clipping.

### Conduction and heating linearization

The conduction solve linearizes `sqrt(U)` and the binary-heating factor
`1/sqrt(U)` about the old state. For a source-fluid energy `U_g`,

```text
1/sqrt(U_g_new) ~= 3/(2 sqrt(U_g_old))
                   - U_g_new/(2 U_g_old^(3/2)).
```

Consequently, a heating source `c4[f][g] Rho[g]/sqrt(U[g])` contributes

```text
A[f,g] += dt c4[f][g] Rho[g] / (2 U[g]^(3/2))
b[f]   += 3 dt c4[f][g] Rho[g] / (2 sqrt(U[g])).
```

An independent dense-system assembly agrees with the production banded solve
exactly for the `c1` and `c4` terms and to `5.6e-17` for the `c2` terms. The
combined test differs by at most `1.7e-16`.

### Three-body convention

The three-body rate follows Spitzer Eq. 6-37. Spitzer's `v_m` is the
three-dimensional velocity dispersion, whereas `sigma` is one-dimensional:

```text
v_m = sqrt(3) sigma,       U = (3/2) sigma^2,       v_m^2 = 2 U.
```

Thus `3^(9/2)/v_m^9 = 1/sigma^9`. Do not retain the explicit `3^(9/2)` while
also substituting `sigma = sqrt(2 U/3)`, because that counts the conversion
twice. With the code's mass normalization, this gives the literal
`0.0009373511756007407` used in `applyBinaryFormation()`.

The implementation and notation follow:

- Yi-Ming Zhong and Stuart L. Shapiro, [*Dynamical Evolutions in Globular
  Clusters and Dwarf Galaxies: Conduction Fluid
  Simulations*](https://arxiv.org/abs/2505.18251).

## Repository layout

```text
.
|-- Makefile             GCC/C++20 release build
|-- docs/                numerical validation notes
|-- external/            pinned Eigen and Boost.PFR submodules
|-- plot_packed.nb       Mathematica import and plotting notebook
|-- script/              Python analysis and plotting utilities
|-- test/                standalone numerical and plotting checks
`-- src/
    |-- main.cpp         experiment definitions and active entry point
    |-- three_fluid.*    simulation state and numerical operators
    |-- io.hpp           binary array I/O
    |-- param.hpp        parameter serialization for Mathematica
    `-- utility.hpp      output-directory and timing helpers
```

## Requirements

The tested Ubuntu build uses:

- GCC with C++20 support,
- GNU Make,
- Eigen,
- Boost.PFR,
- LAPACKE and LAPACK, and
- OpenBLAS.

Eigen and Boost.PFR are pinned as Git submodules. On Ubuntu, the remaining
development packages can be installed with

```bash
sudo apt-get install build-essential liblapacke-dev liblapack-dev libopenblas-dev
```

## Clone and build

Clone the repository together with its submodules:

```bash
git clone --recurse-submodules https://github.com/hypermania/cluster-hydro-sim.git
cd cluster-hydro-sim
make -j6
```

If the repository was cloned without submodules, initialize them before
building:

```bash
git submodule update --init --recursive
make -j6
```

The compiler can be selected through Make's standard `CXX` variable, for
example `make CXX=g++ -j6`. The default build enables `-O3`, native CPU
instructions, OpenMP, `-ffast-math`, and `NDEBUG`, so its executable is not
portable across all CPU architectures and does not retain Eigen's debug
assertions.

Run the strict-IEEE numerical checks and Python tests with

```bash
make check
```

The check constructs discrete equilibria at several resolutions and verifies
that projection changes them only at roundoff level. It also exercises both
production initializers, identical-grid realignment, the positive outer-shell
state, and the fixed outer conduction value. The Statler checks additionally
cover the capture update, bounded evolution, rejected trials and serialization.

## Run

Run from the repository root so the relative output paths resolve correctly:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./main
```

Legacy experiment selection is made by editing calls at the bottom of
`src/main.cpp`. The file defines one-fluid split,
single-fluid tidal, two-component tidal, and binary-formation experiments.
The checked-in entry point currently selects `binary_formation()`. That source
function first runs a no-formation baseline and then runs mode-2 formation.
For the Statler command-line selector, see the final section below.
The non-formation functions provide the baseline examples used to
validate conduction, projection, realignment, and tidal removal.

The current 150-zone mode-2 example reaches the central-density stop at
`t=5.75102` after 42,567 steps. This is a numerically finite trajectory, but
it is not a converged conservative result; see
[`docs/mode2-validation.md`](docs/mode2-validation.md).

`make clean` removes the executable and object files. It deliberately leaves
simulation output in place.

## Output format

A run directory contains parameter metadata and the packed, headerless binary
arrays produced by its attached observers:

| File | Contents | Shape |
| --- | --- | --- |
| `t.dat` | profile snapshot times | `n_snapshot` |
| `R.dat` | common radial grid | `n_snapshot x N` |
| `Rho_s.dat`, `Rho_b.dat`, `Rho_d.dat` | component densities | `n_snapshot x N` |
| `U_s.dat`, `U_b.dat`, `U_d.dat` | component specific energies | `n_snapshot x N` |
| `radii_t.dat` | Lagrangian-radius sample times | `n_radius_sample` |
| `radii_fraction.dat` | enclosed-mass fractions | `n_fraction` |
| `radii_s.dat`, `radii_b.dat`, `radii_d.dat` | component Lagrangian radii | `n_radius_sample x n_fraction` |
| `central_t.dat` | central-value sample times | `n_state` |
| `central_Rho_s.dat`, `central_Rho_b.dat`, `central_Rho_d.dat` | component central densities | `n_state` |
| `central_U_s.dat`, `central_U_b.dat`, `central_U_d.dat` | component central specific energies | `n_state` |

Numeric arrays are native binary `double` values; on the tested x86-64 system
they are little-endian IEEE-754 `float64`. Snapshot times are approximate:
the observer records the first evolved state whose time reaches each requested
value. If a run ends before a requested time, that snapshot is absent. The
Lagrangian-radius observer, by contrast, records every evolution step and can
therefore consume substantial memory in a long run.

The profile observer does not automatically save the terminal state. In the
checked-in mode-2 example, the requested snapshots end near `t=5.567`, while
the run ends near `t=5.751`.

The `central_*` files are currently written only by the no-tidal run. Its
central-value observer is called at the same cadence as the Lagrangian-radius
observer, including the initial and terminal states, so `central_t.dat` and
`radii_t.dat` have the same length and timestamps for that run.

The parameter metadata consists of `param.dat`, `paramNames.txt`,
`paramTypes.txt`, and `paramOffsets.txt`. `param.dat` is a raw C++
`ThreeFluidParam` object, including implementation-dependent layout and
padding; use the accompanying type and offset files when importing it.
Evolution settings are owned directly by `sim.param`; `sim.param.Deltat`
is the initial timestep, while `sim.Deltat` is adaptive state. Initialization
and observer settings are separate records, not fields in `ThreeFluidParam`.
The Statler runner saves those under `initialization/` and `observer/`.
See [the reproduction notes](docs/statler-reproduction.md) for the parameter
split and the optional `param.runtime_validation` full-state checks.

### Plotting packed output

`script/plot_output.py` is a headless Python equivalent of the per-run parts
of `plot_packed.nb`. It reads `N` from the parameter-offset metadata, validates
the packed-array shapes, and makes the following figures when the corresponding
observer files are present:

- log-log density profiles,
- log-log one-dimensional velocity-dispersion profiles, using
  `sigma = sqrt(2 U / 3)`,
- Lagrangian-radius evolution with linear time and logarithmic radius, and
- central-density and central-velocity-dispersion evolution.

In each profile figure, snapshot time is encoded by a rainbow gradient from
the earliest to the latest saved state. Curve style identifies the component:
asterisk markers for single stars, dashed lines for binaries, and solid lines
for dark matter. Long curves automatically subsample the single-star markers
to keep the profiles legible.

The central-value and Lagrangian-radius families are optional, but a partially
written family is reported as an error. Non-finite and nonpositive points are
omitted from logarithmic plots with a warning, which allows the script to be
used for diagnosing a failed trajectory.

NumPy and Matplotlib are required. From the repository root, run

```bash
python3 script/plot_output.py output/one_fluid_binary_formation
```

PDF files are written to the run's `plots/` subdirectory by default. Use
`--formats png` or `--formats pdf png` when raster output is also needed. The
`--plot-directory`, `--title`, and `--max-profile-snapshots` options control
the remaining output choices. The plotting regression test is independently
executable:

```bash
python3 test/test_plot_output.py
```

For example, NumPy can load one run as follows:

```python
from pathlib import Path
import numpy as np

run = Path("output/baseline_AB")
N = int(np.fromfile(run / "param.dat", dtype="<i8", count=1)[0])
t = np.fromfile(run / "t.dat", dtype="<f8")
r = np.fromfile(run / "R.dat", dtype="<f8").reshape(len(t), N)
rho_s = np.fromfile(run / "Rho_s.dat", dtype="<f8").reshape(len(t), N)
u_s = np.fromfile(run / "U_s.dat", dtype="<f8").reshape(len(t), N)
```

`plot_packed.nb` imports the parameter metadata and packed arrays and contains
sections for the one-fluid split, single-fluid tidal comparison, and AB
comparison. It requires an activated Mathematica installation.

## Current validation limits

Before using the results for scientific inference, add invariant and
convergence checks appropriate to the experiment. Important current limits
include:

- the log-density interpolation used during realignment is not
  mass-conservative; in a previously audited mode-2 run, it removed 17.462% of the
  initial mass and dominates the resolved-energy error;
- the standalone hydrostatic check does not yet cover mass conservation,
  energy conservation, complete-run regression, or grid convergence;
- mode-2 formation is locally mass- and resolved-energy-conservative to
  roundoff, but its donor depletion and recipient growth are not included in
  timestep acceptance;
- binary binding energy is not evolved as a separate reservoir, so
  random-plus-gravitational energy is not the complete physical energy when
  binary heating is active;
- the fixed outer-`U` conduction row is a prescribed thermal boundary rather
  than an explicitly conservative zero-flux condition;
- the tidal prescription is a phenomenological explicit density sink rather
  than an external gravitational potential;
- the profile observer does not guarantee a terminal snapshot;
- output writes do not currently report short writes or file-open failures.

The repository does not currently declare a software license.

## Statler comparison

The Statler comparison now uses the shared `evolve()` driver and saved
`ThreeFluidParam` settings, with observer-only output. See
[the setup, numerical scope and reproduction commands](docs/statler-reproduction.md).
