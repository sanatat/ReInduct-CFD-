# ReInduct-CFD Pro v13 — Numerical verification portfolio

**Scope:** a reduced-order porous-gas reactor surrogate (1D axial finite volume + mapped axisymmetric induction solver); **not** a validated OpenFOAM/industrial plasma-MHD reactor.

### Verified changes
- Coupled gas pressure–species–enthalpy Picard correction with a measured EOS closure residual.
- Temperature-dependent NIST Shomate sensible enthalpies and heat capacities for CO, CO2, N2, H2; ideal monatomic H as a separate simple model.
- Nonlinear EM/solid-thermal feedback at field refreshes. Other sources remain lagged over a step.
- Python numerical tests and independent C++17 thermochemical regression.
- 10-scenario parameter/grid/time-step screening, with absolute elemental and energy-ledger checks.

### Reproduce

```bash
python3 -m pip install -r requirements.txt pytest
make v13all
# Optional longer 10-scenario sensitivity study:
make v13study
```

Check `v13/docs/ACCEPTANCE_AND_GAPS_FA.md` for the detailed scientific limitations and source URLs. The v11 OpenFOAM candidate is **not** equivalent to the v13 model and has not been compiled or executed here.
