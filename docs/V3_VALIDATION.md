# v3 Executed numerical verification report (2026-11-05)

**All values below are computed from synthetic model inputs; they are NOT experimental battery recycling data.**

## Verified outputs from `python scripts/run_v3_study.py`

| Diagnostic | Result | Interpretation |
|---|---:|---|
| 2-D MAC flow, mean U | 0.013729 m/s | local Brinkman resistance slows periodic channel |
| Flow maximum transverse velocity | 0.000141018 m/s | nonuniform drag drives a 2-D redistribution |
| Flow max discrete divergence | 1.38e-15 s⁻¹ | projection enforces discrete incompressibility |
| Hartmann benchmark Ha | 1.4142 | idealized imposed field |
| Hartmann profile rel. L2 error | 0.509% | agreement with analytic solution on this grid |
| Coupled model mean NiO conversion | 0.06560 | illustrative conversion at 12 s |
| Coupled model max solid T | 926.49 K | model prediction only |
| Coupled H₂ integral balance | 0.00e+00 relative | exact FV accounting (round-off) |
| Coupled sensible energy balance | 2.37e-15 relative | gas/solid exchange cancels discretely |
| Number of EM updates | 40 | feedback from mean solid T to σ and Joule heat |
| Fixed-σ max solid T | 921.65 K | shows approximated σ(T) feedback sensitivity |

## Mesh sensitivity

| Nx × Ny | Mean conversion | Maximum solid T (K) |
|---|---:|---:|
| 24 × 8 | 0.029719 | 880.42 |
| 48 × 16 | 0.030359 | 898.45 |
| 72 × 24 | 0.030528 | 902.92 |

**Interpretation:** mean conversion shifts less than 1% relative between the final two meshes; however the maximum solid temperature shifts by several kelvin. **The maximum temperature is NOT fully mesh-independent** on the tested range. Do not claim global mesh convergence.

## Time step sensitivity, fixed 36 × 12 mesh, 4 s

| dt (s) | Mean conversion | Maximum solid T (K) |
|---|---:|---:|
| 0.0040 | 0.016630 | 881.729 |
| 0.0020 | 0.016617 | 881.765 |
| 0.0010 | 0.016611 | 881.783 |

This modest numerical variation is a time-discretization check, **not** validation against experiment.

## Additional independent benchmarks and process screening

- `multiregion/conjugate1d.py`: two separate finite-volume material meshes with conservative interface flux; comparison to a series-thermal-resistance analytic solution.
- `plasma/radical0d.py`: H-atom and solid-oxygen stoichiometric tests.
- `scripts/doe_v3.py`: 3 × 3 matrix of inlet gas velocity and induction strength, a **hypothetical** maximum-temperature feasibility cut, Pareto screening, and 3-case activation-energy sensitivity. No industrial optimum is claimed.

## Known limitations

- The flow-to-chemistry link is one-way; no conversion/temperature feedback to fluid momentum.
- The gas–solid calculation is an **interpenetrating local-thermal-nonequilibrium continuum**, not explicit granular multiphase dynamics or separately meshed regions.
- `Ergun_pressure_drop_Pa_diagnostic_only` is post-processing, not a porous momentum sink determined from the particle-size data.
- The harmonic EM model is 1-D and conductivity uses temperature averaged over reactor length.
- Hartmann MHD is a **separate idealized electrically conducting liquid** benchmark, not plasma-gas MHD inside the reactor.
- The 0-D radical model is a prescribed-rate chemistry screen, not a plasma simulation.
- `OpenFOAM` scaffold has not been run here; no validated OpenFOAM reactive case is claimed.
- All kinetics, thermal properties, thermodynamics, electric conductivity, and reaction ΔH are illustrative and uncalibrated.

## Reproduce

```bash
python3 -m pip install -r requirements.txt
make all test emtest v3test
make v3study
make v3radicals
make v3multiregion
make v3doe
```
