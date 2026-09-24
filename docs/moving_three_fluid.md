# MovingThreeFluidSim: implementation and validation status

This is an experimental radial-fluid solver, not a replacement for the
hydrostatic solver. The alternating acoustic flux removes the observed
checkerboard failure in the tested single, split, canonical and concentrated-DM
runs. Physical convergence and unequal-mass agreement remain separate,
incomplete validation requirements. POWER_LAW tidal capture is supported;
other formation modes and explicit stripping are rejected. See
[the formation notes](moving-binary-formation.md) for the source and budgets.

## State and step

Each fixed radial cell stores

```
[u_s,v_s,rho_s,u_b,v_b,rho_b,u_d,v_d,rho_d,M_s,M_b,M_d]
```

All entries are dimensionless. Density, velocity and internal energy are
cell-centred; enclosed mass is at the outer face. `v` is in `r0/t0` units,
not `sqrt(u0)`. Thus `epsilon = r0^2/(u0*t0^2)` multiplies radial inertia
and bulk kinetic energy. The one-dimensional dispersion is `sqrt(2*u/3)`.
The evolution parameters, initialization adapter and output observer are
separate. The existing ThreeFluidSim equations have not been changed.

The solve advances continuity, radial momentum, random-plus-bulk energy,
and the enclosed-mass constraint together. It includes self-gravity,
outward `q*r`, conservative shared-face conduction, POWER_LAW capture and the
existing c1/c4 thermal sources. The square-root/inverse-square-root Taylor expansions
match the conduction derivation. The c1/c4 closures do not include a
separate intercomponent bulk-velocity drag model.

One pivoted LAPACK band factorization solves for primitive increments.
The scalar half-bands are `kl=13, ku=14`, with `ldab=41`. Row/column
equilibration improves scaling without changing the equations. Assembly,
solve, recovery, validation and ledgers are O(N) per step. Work arrays are
reused. Local conservative recovery is required after the linear solve;
simply adding increments to velocity/u does not preserve the conservative
storage update.

The physical centre has even rho/u, odd v, zero face flux/area and the
regular mass relation `M_0 = V_0*rho_0`. It does not delete the first cell's
evolution equations. The outer face uses the gas-to-vacuum Riemann flux;
heat uses an absorbing Robin condition on sqrt(u), with extrapolation
length equal to the domain radius by default. The tidal force-balance
radius is not a fixed wall or a mass-deletion rule.

Mass and advected/thermal energy exports use the **solved affine fluxes**,
not fluxes recomputed from recovered snapshots. The energy/work check
tests random-plus-bulk energy against the included gravity work, heating
and exported energy. It is not a claim of exact total gravitational-energy
conservation. Invalid states raise errors; there are no retries or density
floors. The adaptive controller monitors fractional rho/u changes and
velocity change in thermal-speed units. Final-step overshoot is allowed.

## Corrections discovered during implementation

1. Scalar sound-speed Rusanov viscosity nearly froze the slow contraction.
   The current test scheme uses advective mass/energy diffusion plus a
   frozen velocity viscosity in momentum, including its mechanical energy
   flux. Centred acoustic coupling still produced checkerboarding. It is now
   replaced by left-velocity/right-pressure alternating traces, including
   matching pressure work. The resulting gradient and divergence are an
   adjoint pair without the centred operator's odd/even null space. This
   is first-order accurate, not a proven nonlinear positivity-preserving flux.
2. A frozen reference force density omitted mechanical work and drove
   exterior cooling to zero. The reference correction is now a frozen
   acceleration, with its density-dependent force and corresponding work.
   It balances the supplied q=0 reference but never cancels the q*r force.
3. Frozen primitive reconstruction also corrupted exterior thermal transport
   in the alternating-flux trial. It has been removed completely: advective
   diffusion vanishes at rest and needs no such correction. Reference balance
   now uses the actual right-pressure trace in the source quadrature.
4. Hitting a comparison's safety step limit is reported as failure to reach
   its requested time, not as a completed comparison.

Robust treatment of true empty cells still needs additional work. Preserving the initial
equilibrium and matching a symbolic Jacobian do not establish that property.

## Reproducing checks and plots

Use strict floating-point compilation, including for production experiments:

```sh
make -j2 main-strict check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./main-strict moving-comparison \
    output/moving_single 0 5
python3 script/plot_moving_comparison.py output/moving_single
```

The runner arguments are

```
moving-comparison output_directory sample final_time [epsilon] [max_dt] [canonical_zones] [max_steps]
```

`epsilon=0` requests the initializer's relaxation-time conversion; a
positive override is a numerical experiment, not the same physical system.
`sample` selects existing main.cpp configurations:

- 0: Yiming single-fluid sample (150 cells, trace other components).
- 1: the same fluid split between stars and equal-mass DM (150 cells).
- 2: AB baseline with unequal particle masses and no explicit stripping.
- 3: canonical no-formation baseline (500 cells by default).

The Yiming preset's coefficients and stored particle masses are preserved.
The reference Coulomb logarithm uses the actual total enclosed mass and
stored stellar/DM particle masses. The initializer's `stars` setting sets
the particle count only for the canonical sample. Formation/stripping runs are not
silently relabelled as equivalent source-free runs.

All run/initializer/observer parameters are saved using the shared binary
serializer. Profiles are binary float64 arrays of shape `(zones,3,5)`:
`(cell-centre radius,rho,u,enclosed mass,radial velocity)`. CSV files contain
snapshot metadata and scalar diagnostic histories. Observers alone save
output. Existing output directories containing snapshots are not overwritten.

Plots are PDF. Both solvers use the same time/length/density/dispersion
units, actual accepted times and common physical radii. Snapshot curves
are interpolated in time without extrapolation; no collapse-time shift or
amplitude rescaling is applied. Colours encode time; solid/dashed curves
encode hydrostatic/moving solvers. Components have separate panels.
Partial runs retain their actual stopping time and must not be described
as having reached the requested endpoint.

`check_moving` covers band-vs-dense solution, mass and energy/work ledgers,
equilibrium residuals, tidal acceleration/evolution, equal-fluid symmetry,
reference depletion, fail-early checks, local thermal agreement with the
existing conduction operator and parameter serialization. The separate
Mathematica audit compares every assembled entry and RHS in a five-cell
fixture, including all three outer Riemann branches. These local checks
do not substitute for the long-run comparison gate above.

## Historical centred-flux results (17285f4)

On upstream baseline `49dd74e`, the 150-cell single and equal-split cases
reach code time 5. Their central density agrees with the hydrostatic solver
within 10% for about 62% of the interval; dispersion agrees within 10%
throughout. The two moving calculations agree with each other to about
1e-9 in total central density. The unequal-mass AB-like case differs
substantially by time 2.8. Canonical 500/1000-cell runs track the dominant
stellar fluid closely but reach only 3.1371/2.7291 of the requested time 5
within 12000 steps, owing to trace-fluid oscillations. Refinement does not
resolve this failure.

To reproduce the five reported comparisons, run these commands from Code
with `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` exported. Use fresh
output directories. The last two commands deliberately return failure at
the step limit, while preserving their partial histories.

```sh
./main-strict moving-comparison output/moving_comparison/matched_single 0 5
./main-strict moving-comparison output/moving_comparison/matched_split 1 5
./main-strict moving-comparison output/moving_comparison/matched_ab 2 2.8
./main-strict moving-comparison output/moving_comparison/review_canonical 3 5 0 0.001 500 12000
./main-strict moving-comparison output/moving_comparison/review_fine 3 5 0 0.001 1000 12000
python3 script/plot_moving_comparison.py output/moving_comparison/matched_single output/moving_comparison/matched_split output/moving_comparison/matched_ab output/moving_comparison/review_canonical output/moving_comparison/review_fine
python3 script/plot_moving_comparison.py --split-check output/moving_comparison/matched_single output/moving_comparison/matched_split
```

The existing solver is not a conservative ground truth: `realign()`
interpolates density and rebuilds pressure/energy. The independent diagnostic
`OPENBLAS_NUM_THREADS=1 ./check_moving --audit-hydro` measures its mass,
random-energy and entropy-proxy changes, including a nonzero boundary
modification even at zero timestep. This does not explain away the new
solver's instability; the two issues need separate convergence tests.

## Alternating-flux checks (22 September 2026)

The new single and split runs reach t=5 in 5537 steps; canonical 500/1000-cell
runs finish in 7237/7317 steps, and the AB-like run reaches t=2.8. The single
fluid's core second-difference velocity estimator falls from 0.119 to
1.92e-4. At 500/1000 cells the final stellar estimators are 1.33e-5/3.05e-6;
the remaining signal includes smooth curvature, not an alternating mode.
All components are inspected, not only those above a mass-fraction cutoff.

The heating-off concentrated-DM regression (rho0=1, xi2=1, zeta2=0.3,
ms=1e-6, md/ms=1e-10, other inputs saved by its initializer) finishes t=5
at 150 and 300 cells. Its estimator drops from about 0.53/1.26 in the old
scan to 1.93e-4/3.05e-5. Central densities still depend on resolution:
removing the instability does not establish converged core collapse.
The coarse unequal-mass comparison also still differs from the hydrostatic
solver. At 1000 cells the dominant stellar density and dispersion in the
canonical comparison agree within 10% over the tested interval.

New short tests perturb pressure and velocity separately with conduction
off, at epsilon=0.03 and 1e-9. At acoustic Courant number 10, their measured
one-step alternating amplitudes are below 0.061 of the initial amplitude.
A separate no-conduction/no-heating evolution reaches t=0.12. Tests retain
the mass and gas-energy/work ledgers. No state floor, retry or velocity
filter is used. The outer vacuum and thermal boundary formulas are unchanged.

Example reproduction commands (fresh directories, single BLAS/OMP thread):

```sh
./main-strict moving-comparison output/new_single 0 5
./main-strict moving-comparison output/new_canonical 3 5 0 0.001 500 16000
./main-strict moving-comparison output/new_fine 3 5 0 0.001 1000 16000
./check_moving --contraction 150 output/new_contraction150
./check_moving --contraction 300 output/new_contraction300
python3 script/plot_moving_oscillation_check.py output/old_single/moving output/new_single/moving output/velocity_comparison
```

The last command requires saved old (17285f4) output. Its JSON/PDF compare
raw signed velocities at common time without extrapolation. Velocity is in
r0/t0 units; the physical ratio to dispersion is sqrt(epsilon)*v/sqrt(2u/3),
not v/sqrt(2u/3). The estimator is grid-dependent and must not be presented
as a fixed-wavenumber spectral convergence test.
