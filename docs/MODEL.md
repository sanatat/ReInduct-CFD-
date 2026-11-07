# Reduced-order reactive-bed model: governing equations and limitations

All state variables are 1D axial cell averages in a stationary porous bed, not resolved particle-scale flow or Navier–Stokes CFD. The advection of gas through the bed is prescribed by a constant **superficial** velocity U; porosity ε is constant. Gas-phase species are in mol/m³ **gas**, bulk NiO inventory n0 is mol/m³ **total reactor**.

### H2 and H2O conservation

ε ∂c_H2/∂t + U ∂c_H2/∂z = −r,
ε ∂c_H2O/∂t + U ∂c_H2O/∂z = +r,

where `r = n0 * A exp(-Ea/(R Ts)) * (c_H2 / c_ref) * (1 − X)` and is capped by available gas H2 and unconverted NiO per time step.

### Solid conversion

n0 ∂X/∂t = r,  0 ≤ X ≤ 1.

### Temperature

ε Cp_g ∂Tg/∂t + U Cp_g ∂Tg/∂z = h_v (Ts − Tg),

Cp_s ∂Ts/∂t = k_eff ∂²Ts/∂z² + q_ind(z) − h_v (Ts − Tg) − h_loss (Ts − Tamb) − ΔH_r r.

`Cp_g` is per gas volume; `Cp_s` is per bulk reactor volume. Assumed heat deposition `q_ind(z) = q_peak exp(-((z − 0.065)/0.04)²)`. The name is suggestive only; **no induction or plasma physics is solved**. End boundaries are adiabatic for solid; gas uses upwind inlet advection.

The reaction is only a simple pedagogical first-order apparent rate law. NiO reduction may exhibit nucleation, autocatalysis, particle-size dependence, diffusion limitations, water inhibition and other effects. The model uses fixed bulk heat capacities and ignores pressure, detailed gas thermodynamics, radiative transfer and porous permeability coupling.

### Electromagnetic scalar reference

Under simple homogeneous good-conductor, non-magnetic assumptions, classical field penetration depth is δ = √(2/(ω μ σ)), ω=2πf. Evaluated only as a standalone dimensional check in the C++ output. It is **not** connected to the prescribed axial heat source and does not demonstrate Maxwell or MHD coupling.

### Physical verification and validation distinctions

- **Conservation test**: hydrogen in minus out should equal gas inventory plus chemically reduced solid inventory. This is an algebraic numerical verification.
- **Grid and time-step check**: differences should decrease when resolution is improved, for fixed physical assumptions.
- **External validation**: requires an independent experiment or trusted benchmark with specified uncertainty. None is included.
- **Thermodynamic/kinetic calibration**: illustrative ΔH = −2500 J/mol and Arrhenius parameters are NOT fitted to NiO data. Replace using literature or experimental parameters with matching rate-law dimensions before quantitative conclusions.

### Parameters currently used

| Name | Value | Unit | Meaning |
|---|---:|---|---|
| L | 0.12 | m | model length |
| ε | 0.45 | – | void fraction |
| U | 0.016 | m/s | superficial gas speed |
| n0 | 150 | mol/m³ bulk | active NiO inventory |
| Cp_g | 450 | J/m³ gas/K | gas volumetric heat capacity |
| Cp_s | 1.188e6 | J/m³ bulk/K | solid volumetric heat capacity |
| h_v | 800 | W/m³/K | gas-solid heat exchange |
| k_eff | 1 | W/m/K | axial effective solid thermal conductivity |
| h_loss | 35 | W/m³/K | bulk heat leak factor |
| T0, Tin, Tamb | 850 | K | initial, inlet, ambient temperatures |
| A | 1500 | 1/s | assumed pre-exponential factor |
| Ea | 85000 | J/mol | assumed activation energy |
| c_ref | 15 | mol/m³ gas | reference gas concentration |
| q_peak | 600000 | W/m³ bulk | nominal peak prescribed heating |
| ΔH_r | −2500 | J/mol | nominal assumed reaction heat |

## Numerical method

First-order conservative upwind finite-volume advection, forward Euler for gas species, thermal exchange, solid conduction and reaction. Explicit reaction is positivity-limited. Mass conservation follows analytically from matched H2 source and solid conversion, including a conservative face flux. Physical accuracy is limited by the reduced-order constitutive assumptions.
