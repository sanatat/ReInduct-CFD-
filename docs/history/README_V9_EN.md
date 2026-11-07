# ReInduct-CFD Pro v9 — Auditable simulation portfolio

**Scope:** reduced-order induction-heated reactive gas/solid bed, hydrogen chemistry energy ledger, and kinetic LiCoO2/C process modelling. Earlier v4–v8 modules remain included. **Not** a validated digital twin.

## Run

```bash
python3 -m pip install -r requirements.txt
bash v9/run_v9.sh
```

For OpenFOAM.com v2606 (separate, **unverified** locally):

```bash
bash v9/openfoam/Allrun
```

`v9/results/openfoam_attempt.log` records local exit code 2 for absent `wmake`. The GitHub Actions workflow `.github/workflows/v9-verification.yml` must be triggered in a GitHub repository to establish actual runtime compatibility.

## What is new

- `v9/plasma_energy.py`: simplified H₂ ⇌ 2H energy-coupled batch; bounded reaction rates, H atom and electron+neutral+chemical energy ledgers.
- `v9/lco_process.py`: carbon-constrained LiCoO₂ batch reaction with qualitative phase selection and mass conservation; predicted curves are illustrative.
- `v9/validation_protocol.py`: four held-out-temperature within-study kinetic diagnostics on approximate literature endpoint values; **not external validation**.
- `v9/external_validation.py`: refuses absent/nonindependent/mismatched reactor observations; requires provenance and uncertainty.
- `v9/tests/test_v9.py`: new numerical and provenance tests; prior tests retained.
- `v9/results/`: actual generated plots and logs.

See `v9/docs/SCIENTIFIC_VALIDATION_AND_LIMITATIONS.md` for explicit limitations and the pass/fail criteria for each missing research skill. Scientific sources include DOI `10.1007/s11663-022-02712-1` and the unrelated H₂ cathode kinetics study https://pmc.ncbi.nlm.nih.gov/articles/PMC9548595/.
