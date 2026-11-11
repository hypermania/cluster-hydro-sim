# Mass-scaled tidal cutoff

Select the hydrostatic cutoff through `ThreeFluidParam::tidal_cutoff`:

| Flag | Behavior |
| --- | --- |
| `TIDAL_CUTOFF_OFF` (0) | No density sink |
| `TIDAL_CUTOFF_ON` (1) | Existing fixed-radius sink |
| `TIDAL_CUTOFF_MASS_SCALED` (2) | Radius follows total cluster mass |

For the new mode, `param.tidal_radius` is the **initial** radius `r_J0`.
At each call to `applyTidalCutoff()`, before any density is modified,

```text
M_cluster = sum_f Menc[f][N-1]        // M_tot[N-1]
r_J       = r_J0 * cbrt(M_cluster / M_cluster_initial)
```

The initial reference includes stars, binaries, and DM, and is captured on
the first cutoff call at time zero. In the existing step ordering, the
preceding conduction operation does not change mass. Every later call uses
the updated enclosed masses, including the previous step's drain and
hydrostatic realignment. The same radius is used for every component; the
reference radius is never overwritten or multiplied cumulatively.

The density update remains the existing explicit prescription:
`rho *= 1 - tidal_cutoff_factor * Deltat` for fully exterior shells,
with the existing radial-width weight for a shell straddling the radius.
The innermost zone remains excluded, as in the fixed-radius implementation.
Specific energy and radii are unchanged by the sink itself; pressure and
enclosed mass are recomputed before the usual hydrostatic projection.
This adds no external force or radial-velocity equation.

After initialization (which may reset physics flags), configure e.g.:

```cpp
sim.param.tidal_cutoff = TIDAL_CUTOFF_MASS_SCALED;
sim.param.tidal_radius = 20.0;
sim.param.tidal_cutoff_factor = 10.0;
sim.evolve(observer);
```

All quantities are dimensionless. The existing settings serializer records
the flag, initial radius, and drain rate without changing the parameter
layout. `sim.tidal_initial_mass` is runtime state, not another tunable
parameter. It persists when continuing the same simulation object and is
cleared by `initSolver()`. A restart must restore it together with the
profiles and time; parameter files alone are not checkpoints. Enabling the
new mode at nonzero time without supplying the original reference mass
throws instead of silently redefining the initial mass.

Invalid modes, nonpositive radii or masses, invalid drain rates, and steps
with `lambda * Deltat >= 1` fail rather than retry. The current mass used
in the scaling includes **numerical mass drift**: this change does not
repair non-conservative realignment or add an escaped-energy ledger.

Run the operator and evolution regression checks independently with:

```sh
make output/check_tidal_cutoff
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./output/check_tidal_cutoff
```

They are also included in `make check` and use strict floating-point
semantics (no `-ffast-math`).
