# ReInduct-CFD Pro v5 — evidence register, validation, and remaining skill gaps

**Generated from local execution: November 5, 2026.** This file distinguishes verification of code from validation of experimental reactor physics. The exact command `make v5all` regenerates the datasets and plots.

## Computed numerical evidence

| Check | Result | Meaning |
|---|---:|---|
| V5 unit tests | 13 passed | EM manufactured solution, q positivity, mesh, kinetic holdout, thermo, MHD, thermal/chemical invariants |
| V4+V5 tests combined | 23 passed | Includes C++17 kernel compilation in legacy v4 tests |
| EM manufactured solution relative L2, 18×12 | 1.2604e-3 | 2-D complex Maxwell-diffusion operator verification |
| EM manufactured solution relative L2, 36×25 | 3.1676e-4 | ~4-fold reduction |
| EM manufactured solution relative L2, 72×50 | 8.2275e-5 | ~4-fold reduction |
| Coil model integrated power, 51×25 | ~5702 W/m | Smooth conductor sources; integrated over out-of-plane metre |
| Coil model integrated power, 81×40 | ~5715 W/m | Mesh difference ~0.22% |
| NiO growth holdout: simple Arrhenius | ~0.267 µm/s MAE | 26% H2 held out; microkinetic experimental observable |
| NiO growth holdout: two-regime first order | ~0.091 µm/s MAE | Exploratory post-comparison choice; risks overfitting |
| NiO/H2 Kp at 800 K | ~335 | idealized reference thermochemistry, not measured in reactor |
| With EM: max bed Ts (1.2 s) | ~876.63 K | representative properties, NOT industrial prediction |
| With EM vs none: max Ts delta | +26.62 K | one controlled simulation contrast |
| With EM: mean X after 1.2 s | 0.00072443 | model fraction, NOT recovered battery yield |
| Coupled thermal H balance residual | ~4.6e-16 relative | numerical face-flux bookkeeping |
| Coupled thermal energy residual | ~1.4e-14 relative | numerical thermal source/face-flux bookkeeping |
| MHD proxy with sigma_g=5 S/m | mean velocity ~3.84% lower vs without | illustrative low-Rm conductivity assumption |
| OpenFOAM v2606 compilation | **not run** | No `wmake`/Docker in this runtime |
| Reactor experimental data calibration | **not done** | no reactor measurements provided |

These small numerical residuals demonstrate exact discrete bookkeeping. They do not establish chemical, pressure, thermodynamic or kinetic realism. The 2-D flow uses periodic axial hydrodynamics and inlet/outlet species boundary conditions, an intentionally simplified proxy. The NiO rate constant `1500` and apparent activation energy `85000 J/mol` used in the spatial reactor have **not** been calibrated against measured volumetric conversion data; measured dense NiO layer growth rates have different dimensions and morphology.

## Physics implemented and not implemented

### 2-D electromagnetic induction

Vector potential `A_z` phasor, long-conductor translational geometry:

`(-∇² + i ω μ σ) A_z = μ J_coil,z`

`E_z = -iω A_z; B_x = ∂A_z/∂y; B_y = -∂A_z/∂x`

`q_Joule = (σ/2) |E_z|²`, `σ(T) = σ_ref / (1 + 0.003 (T-850 K))` for demonstrator only.

This is a mathematically verified Cartesian long-conductor approximation, not an axisymmetric/helical coil model. Independent material conductivity, skin depth and measured coil power are absent.

### Thermochemical NiO/H2

`NiO(s) + H2(g) -> Ni(s) + H2O(g)`.

Gas thermochemistry uses NIST Shomate equations in the domain 500–1000 K. Solids use JAEA NiO/Ni sensible heat and entropy increments shifted to published standard anchors (NiO -239.7 kJ/mol at 298 K and S=37.99 J/mol/K; Ni S=29.87 J/mol/K). This **mixed-source approximation** must not be used outside its tabulated 500–1000 K range without a consistent reference-state audit. The driving term `1 - Q/Kp`, `Q=pH2O/pH2` is limited to ideal gas with unit activity solid phases; the volumetric kinetic rate remains a surrogate.

### Kinetics holdout

Hidayat et al. 2009, DOI 10.1007/s11663-008-9212-0, Table III gives initial **Ni product-layer growth [µm/s] in 0–5 seconds** on dense synthetic NiO. These are *not* conversion percentages or molar rates. A freely-fitted temperature-dependent Arrhenius model trained on 53%/75%/100% H2 is compared with a two-segment temperature law, including a *post-error-examination* hypothesis of first order in hydrogen. The holdout is a 26% H2 composition condition in the same paper (4 temperatures); not external replication or independent facility experiments. At a high bar, future validation should use another study, different oxide morphology and gas pressure.

### Low magnetic Reynolds number MHD proxy

Assumed out-of-plane current closure and prescribed complex `B` field give `J_z=σ_g (u_x B_y-u_y B_x)`, and `f_x= -0.5 Re(J_z conj(B_y))`, `f_y=0.5 Re(J_z conj(B_x))`. The flow includes the source divided by an assumed constant ρ_g=0.1 kg/m³ and returns dissipated mechanical energy into gas heat. Magnetic induction in gas and electric current continuity in finite 3-D electrodes are not solved. A constant conductivity σ_g=5 S/m is a *hypothetical partially conducting gas*, not a chemistry/ionization-derived plasma result.

### OpenFOAM

The v2606 C++ source contains velocity–pressure PISO, H2/H2O, gas/solid temperatures, X and imported `qEM` (360 cells). This is **static source preparation only**. The Docker launcher attempts to actually compile with `wmake`, generate `blockMesh`, `checkMesh`, and run. The blocked local attempt is in `v5/results/openfoam_env_check.log`. Upon obtaining `wmake`, compilation or solver errors may require API corrections. The .github workflow is committed as a configuration, not as an executed successful CI job.

### Real black mass chemistry

The literature-derived LiCoO2/C overall equation provides elemental and gravimetric balances. No measured NMC/LFP cell composition, measured carbon loss, solvent/binder properties, or electrode metal recovery has been used; rate laws and product separations are missing.

## Remaining gaps to claim the position's complete technical skillset

| Priority | Gap | Required external evidence |
|---|---|---|
| P0 | OpenFOAM full build & run | wmake compile, mesh and solver logs, field output, Courant/mesh refinement |
| P0 | Reactor/battery experimental validation | independently measured conversion, gas composition, sample geometry and uncertainty |
| P1 | Full plasma discharge or self-consistent magnetohydrodynamics | ionization and charge equations, electric potential, conductivity and species/energy closures, comparison to benchmark |
| P1 | Motion-resolving gas–solid multiphase (if intended) | Euler–Euler or Lagrangian particles with drag/exchange and verified phase balances |
| P1 | Real coil/reactor geometry | measured winding dimensions, skin effect, temperature-dependent property datasets, wall losses |
| P1 | Fully validated NMC/LFP chemistry | multi-step reduction and Li/Co/Ni/Mn recovery kinetics with experimental reference |

## Files reproducible without OpenFOAM

- `v5/em2d.py`, `v5/mhd2d.py`, `v5/thermochemistry.py`, `v5/kinetics_validation.py`
- `v4/coupled_reactor.py` updated with optional `em_backend='2d', thermo_backend='tabulated', mhd_sigma_g=5`
- `v5/results/*`: raw data, metrics, figures; `v5/run_study.py` regenerates evidence
- `v5/openfoam/export_foam_heat.py` writes 0/qEM; `v5/openfoam/run_local_or_docker.sh` runs real OpenFOAM when available

### Training-project CV statement

> Built and numerically verified a multiphysics reactor research prototype in Python/C++ with 2D quasistatic electromagnetic induction and temperature-dependent Joule heating coupled to reactive finite-volume gas–solid transport, checked a simplified MHD Lorentz-force formulation, and evaluated NiO growth kinetics on published concentration-holdout measurements; prepared a conditional OpenFOAM v2606 solver/CI port (runtime validation pending).
