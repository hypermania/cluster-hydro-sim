# Heggie & Aarseth (1992): primordial-binary gas tests

The target is *Dynamical effects of primordial binaries in star clusters – I.
Equal masses*, MNRAS **257**, 513–536, specifically the gas models of Section 3
and Figs. 1–4, **not** the later N-body models. Supply your own reference PDF.
The companion [TeX note](heggie-reproduction.tex) derives the coefficients;
[Mathematica](../script/mathematica/heggie_units.wls) checks their units.

## Physics and assumptions

- Figs. 1–2: single-component gas, no segregation; prescribed heating strength
  `f_star = f(1-f)/(ln(Lambda)*(1+f)^2)` takes .001, .002, .005, .01, .02.
- Fig. 3: two components with particle-mass ratio 2 and total mass ratio .01;
  conduction and equipartition exchange, no encounter heating.
- Fig. 4: initial binary **number** fraction .03, with binary–single heating
  alone or binary–single plus binary–binary heating.
- No formation, disruption or encounter ejection in any of these models.
  All start with a Plummer profile and equal velocity dispersions.
- Reflecting, insulated outer sphere, as in the gas models. Default radius
  is 10000 Plummer scales: the paper only specifies a large radius.
- Conductivity C=.104 (Heggie–Ramamani calibration); ln(Lambda)=10 for the
  two-component heating comparison. The latter is an explicit assumption:
  Section 3.2 does not provide an unambiguous value for the plotted runs.
  N_star=300000 sets only the small dynamical/relaxation time ratio here;
  it is not the 2500-star N-body setup later in the paper.
- Gas figure captions do not explicitly restate absolute density/length
  units. We use G=M=-4E0=1. Heggie & Ramamani (1989), p. 760, explicitly
  adopts this convention for the earlier gas implementation cited in the
  1992 paper. This supports, but does not independently establish, the
  normalization of the original Fig. 4 data. No normalization is fitted;
  Fig. 3 ratios are independent of it.

Evolution settings remain dimensionless in `MovingThreeFluidParam`.
`HeggieInitParam` and `HeggieObserverParam` separately own initialization and
observable conversions. The reusable solver additions are an insulated
reflecting boundary and heating proportional to `(U_f+U_h)^(-1/2)`.
Defaults, including the Statler donor-dispersion heating convention, are
unchanged. No paper-specific branches or time loops enter the solver.

## Reproduce

```sh
make -j2 main-strict
python3 script/run_heggie_suite.py output/heggie400 --zones 400 --jobs 3
python3 script/plot_heggie_comparison.py path/to/Heggie_1992.pdf output/heggie400
```

Choose fresh directories. The suite records all eight run logs; a numerical
failure is nonzero exit, not a successful endpoint. Unheated collapse and
insufficiently heated models may reach the configured central-density ceiling
of 1e9; this is recorded as `DENSITY_LIMIT`, not extrapolated to later times.
The individual command is
`./main-strict moving-heggie OUTPUT MODEL [FSTAR] [ZONES]`, with model 0
for Figs. 1–2, 1 for Fig. 3, 2 for BS only, 3 for BS+BB.

`history.dat` is little-endian float64, 13 entries per accepted state:
code time, time/trh0, central singles density, central binaries density,
positive central potential depth, core radius, half-mass radius, core mass
fraction, total code mass, singles central relaxation time, binaries central
relaxation time, gas-energy ledger error, current code timestep.
Densities/radii use G=M=1 and initial energy -1/4; potential and relaxation
times retain code units (only their ratios/products are plotted).
All effective settings and conversion factors are saved independently.

`plots/heggie_1_4_overlays.pdf` and individual PDFs preserve the scanned
axis-frame aspect ratio and limits, with no fitted coordinate shifts.
Fig. 3 differentiates log density in **code time** and multiplies by the
corresponding **code-unit** relaxation time. The gray background is the
reference scan, not digitized reference measurements. `summary.json`
records endpoints, stopping conditions, and discrete budget errors.

Short tests cover wall fluxes, heating matrix/RHS, initialization, observable
normalization, and conservation; they do not establish long-run convergence.

## Initial 400-cell results

All eight runs terminate without a numerical exception. The .001 and .002
single-component cases reach the density ceiling at 16.684 and 17.627 trh;
the .005, .01 and .02 cases reach 1000 trh. The unheated segregation case
reaches its density ceiling at 13.644 trh (its Fig. 3 comparison interval is
covered). The BS-only case instead reaches the ceiling at 12.620 trh, **not
the bounce shown in the reference**. BS+BB reaches 150 trh.

Across the suite, relative mass drift is below 1.6e-15 and accepted-step
gas-energy ledger residual below 5.1e-17. Fig. 3 is qualitatively close;
Figs. 1–2 develop late behavior absent from the reference, and Fig. 4 is not
quantitatively reproduced. These discrepancies must not be hidden by fitted
axis scalings or described as validated physics. Core transport and the
reference-force correction require convergence checks; reference absolute
density normalization is also an unresolved input for Fig. 4.

For a resolution check, run the same model at 800 zones and use
`python3 script/compare_heggie_resolution.py COARSE FINE OUTPUT.pdf`.

The completed 800-cell checks also reach their endpoints, but **do not
establish convergence**. For f*=.005, the peak singles density changes from
374.12 to 272.96 and the minimum core-mass fraction from .00433 to .00119.
For BS+BB, the peak changes from 4.506 to 5.083, and the final singles density
from .00824 to .03941. Thus the disagreement cannot yet be attributed solely
to differing physical models. In particular, the initial reference-force
correction is frozen at fixed radii; as the core expands and its gravity
weakens, its relative importance can grow. Auditing this correction and
discrete gravitational work is a concrete next numerical check, not a
demonstrated explanation of the entire discrepancy.

## 2000-cell Fig. 4 density audit

The audit uses revision `2e49bce` and the completed 2000-cell BS+BB history.
In standard N-body units the Plummer scale is a=3*pi/16 and the total
initial central density is 3/(4*pi*a^3)=1.16804. With binary number fraction
f=.03, the single and binary mass fractions are (1-f)/(1+f) and 2*f/(1+f).
The saved initial central densities, 1.09999946 and .06804120, agree with
these conventions. No missing constant conversion factor was demonstrated.

Fig. 4 begins at t/trh(0)=1, not at the unevolved initial state. Approximate
readings of the solid reference curves give:

| t/trh(0) | HydroSim singles | HydroSim binaries | Reference singles | Reference binaries |
| --- | ---: | ---: | ---: | ---: |
| 1 | 1.39560 | .243745 | 3.6 | .6 |
| 10 | 5.34491 | 5.49701 | 10.8 | 9.5 |

Reference values are scan estimates, not original numerical data. The
left-edge offset resembles a common factor of about 2.5, but the required
factors near the peak differ between components and are smaller. A common
vertical rescaling therefore cannot reconcile the entire evolution or the
component density ratio. The cause remains unresolved: reference settings,
transport/heating accuracy and spatial convergence need further checks.
Calling this a quantitatively matched reproduction would be premature.

The scanned horizontal decade ticks also suggest a right frame endpoint
near 158, whereas the overlay currently uses 150. This small horizontal
registration uncertainty cannot explain the factor-of-two density offset.
No rescaling, solver change or new evolution run was made for this audit.

Source: [Heggie & Ramamani (1989), MNRAS 237, 757, p. 760](https://adsabs.harvard.edu/pdf/1989MNRAS.237..757H).
