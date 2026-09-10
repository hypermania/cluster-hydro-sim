# HydroSim codebase rules

These rules apply to subsequent codebase changes, including examples,
reproduction runners, diagnostics and tests. Read this file before making
changes. Existing code is not necessarily compliant: violations are work
to address, not precedents to copy. Do not perform unrelated rewrites merely
to enforce these rules.

These rules have been reviewed by the project owner. Changes require the
owner's agreement; do not silently introduce exceptions.

## 1. Save every run setting in `ThreeFluidParam`

- Every configurable choice needed to reproduce a run must be represented
  in `ThreeFluidParam` and included in the saved parameter record.
- Include physical scales and units, particle masses, initial-profile type
  and parameters, grid settings, enabled processes, source coefficients,
  boundary conditions, tidal settings, solver tolerances, timestep controls,
  stopping conditions and random seeds when applicable.
- Include observer configuration that determines the saved result: sampling
  times or cadence, selected diagnostics, Lagrangian mass fractions and
  output format. Do not leave reproducibility-critical settings only in
  runner literals, comments, command-line arguments or observer constructors.
- Command-line arguments and presets may populate the configuration; they
  must not bypass its serialization. Save the effective configuration after
  initialization and coefficient overrides, before evolution changes state.
- Avoid independently maintained copies of the same setting. If simulator
  members mirror parameter fields, use a single explicit import/export path
  and test that values agree. Update serialization and readers when adding
  a field; never silently omit a value because its type is inconvenient.
- Separate configuration from evolving state. For example, record the
  initial timestep as a setting; the current adaptive timestep belongs to
  state/history. A restart needs both configuration and the required state,
  not just a parameter file.

## 2. Every evolution runner calls `evolve()`

- Examples, benchmarks and reproduction runners configure the simulation,
  initialize it, construct observers, and call `evolve()`.
- Do not implement another time loop in a runner or manually chain the
  physical operators as an alternative production evolution path.
- If a study requires a new process, stopping rule, integration mode or
  accepted-step diagnostic, extend or generalize the shared `evolve()`
  infrastructure. Do not copy it into a study-specific driver.
- Single-operator unit tests may call the operator being tested directly.
  Comparative evolution tests still use `evolve()`; this exception is not
  permission to maintain a second solver in a test or diagnostic.

## 3. Keep `evolve()` as orchestration

- `evolve()` coordinates the numerical steps and their ordering. Put
  logically separate operations into named functions, not inline blocks
  inside the evolution loop.
- This includes conduction/exchange, binary formation, tidal loss,
  hydrostatic projection, realignment, timestep selection, stopping checks,
  validation, and any future trial-step acceptance or rollback.
- Functions must state their inputs, modified state, preconditions and
  intended invariants. Names and comments must describe the actual operation.
- Keep the physical step order explicit in one place. Treat a change in
  ordering, iteration count or acceptance policy as a numerical change,
  not merely a refactor.
- Helpers must operate on the shared simulation state and configuration;
  do not hide study-specific constants or alternate physics in them.

## 4. Pass saving and output through observers

- All run-output behaviour supplied to `evolve()` must enter through its
  observer interface. Do not add filenames, output streams, plotting logic,
  study-specific diagnostics or save calls to the numerical operators.
- Compose observers for profiles, scalar diagnostics, Lagrangian radii,
  checkpoints and progress reporting rather than duplicating evolution.
- Observers inspect a read-only simulation state. They must not change the
  physics, timestep or stopping policy as a side effect of observation.
- Define whether observations represent the initial state, accepted steps
  or the terminal state. If retries are introduced, rejected trial states
  must not silently enter accepted-history output.
- The runner may create output directories and invoke observer saving or
  finalization outside `evolve()`. Parameter saving must use the shared
  serializer, not a hand-written runner-specific copy of the settings.
- Error reporting must remain possible without a saving observer. Return
  errors or throw exceptions for the runner to report; diagnostic failure
  reporting is not a reason to add scientific output to the solver.

## 5. Keep runners declarative and features reusable

- A runner should explain a physical setup through configuration, not
  implement numerical methods. Keep initialization in reusable functions.
- Add optional physics through explicit parameters and named operations.
  Do not introduce paper-name conditionals into generic solver routines.
- Do not edit or generate temporary copies of production source merely to
  change rates, coefficients, timestep logic or output for a study.
- Keep the public repository sufficient to build and reproduce published
  examples. Do not depend on absolute server paths, untracked helper code,
  private patches or undocumented local environments.
- Describe implemented behaviour separately from proposed methods. Pin the
  code revision used for scientific comparisons and record deliberate
  differences from a reference calculation.

## 6. Make numerical contracts explicit

- All physical quantities in the simulation's internal state are
  dimensionless. This convention is not formally enforced by the code or
  type system; review it explicitly whenever state variables, physical
  processes or numerical operations are added or changed.
- Document the reference scales used to nondimensionalize each quantity.
  Convert dimensional inputs before they enter the internal state; keep
  conversions back to physical units at the analysis/output boundary.
  Dimensional reference scales may be saved as configuration metadata,
  but must not be confused with dimensionless evolving state.
- State units and normalizations for parameters, rates and output fields.
  Distinguish particle mass from enclosed mass, specific energy from energy
  density, and binary centre-of-mass energy from internal binding energy.
- Keep initialization, projection and pressure reconstruction consistent
  with the same discrete equations and boundary conventions.
- State conservation claims per operation. Exact shell mass and entropy
  preservation during relaxation does not imply conservative realignment
  or full-step energy conservation.
- Account separately for physical sinks/sources, boundary fluxes and
  numerical drift. Do not interpret remapping losses as physical escape.
- Validate finite values, positivity and shell ordering at appropriate
  boundaries. A numerical failure must not return success status.
- Make stopping and restart semantics explicit. Guard timestep selection
  against zero error estimates, invalid proposals and unintended endpoint
  overshoot. Configuration validation belongs in reusable code.

## 7. Require tests and reproducibility with each change

- Add or update tests for changed behaviour and each fixed bug. Include
  control-flow and output tests, not only matrix or equilibrium tests.
- For a nominal refactor, compare bounded numerical histories against the
  previous revision under identical settings. Require bitwise agreement
  where arithmetic is unchanged, otherwise document justified tolerances.
- Changes to parameter layout or observers need round-trip/readback tests.
  New physics needs operator-level checks and an integration check through
  `evolve()`.
- Run the applicable C++ and Python suites, including after each commit.
  Use strict floating-point semantics for finiteness and conservation tests;
  an optimized build passing is not a substitute for those checks.
- Keep short deterministic checks practical to run routinely. Distinguish
  smoke tests from long-run convergence and scientific validation.
- Record the tested revision, build options, configuration and commands.
  Keep generated data and binaries out of Git; retain essential source,
  tests and reconstruction instructions.

## Change checklist

Before handing off a code change, check:

1. Are all new run settings saved in `ThreeFluidParam`?
2. Does every evolution runner still use the shared `evolve()`?
3. Are distinct operations implemented as named functions?
4. Is saving/output observer-based and observational only?
5. Are all internal physical state quantities dimensionless, with their
   normalizations, invariants, boundary conditions and failure semantics
   explicitly reviewed?
6. Do the relevant tests pass, with evidence for any behaviour change?
7. Can the public version reproduce the result without server-local code?
