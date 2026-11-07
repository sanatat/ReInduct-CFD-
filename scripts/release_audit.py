"""Validate repository publication hygiene and fail-closed verification status.

Run from project root: python3 scripts/release_audit.py
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = [
    'LICENSE', 'README.md', 'README_FA.md', 'requirements.txt',
    '.github/workflows/numerical-verification.yml',
    'docs/REPRODUCIBILITY.md', 'docs/SCIENTIFIC_LIMITATIONS.md',
    'docs/figures/relative_sensitivity.png',
    'v13/reactor.py', 'v13/thermo_pressure.py', 'v13/em_feedback.py',
    'v13/tests/test_v13.py', 'v13/openfoam/Allrun',
]

def audit() -> dict:
    missing = [name for name in REQUIRED if not (ROOT / name).is_file()]
    if missing:
        raise AssertionError(f'Missing release essentials: {missing}')
    evidence_file = ROOT / 'v13/results/EVIDENCE_MANIFEST.json'
    if not evidence_file.exists():
        raise AssertionError('Numerical evidence has not been generated. Run: make v13all')
    manifest = json.loads(evidence_file.read_text(encoding='utf-8'))
    forbidden_claims = ['openfoam_v13_runtime_verified', 'full_plasma_mhd_verified',
                        'experimental_reactor_validated','experimental_black_mass_recovery_validated']
    if any(manifest.get(k) is not False for k in forbidden_claims):
        raise AssertionError('Evidence manifest must explicitly deny unsupported validations')
    if not manifest.get('source_sha256'):
        raise AssertionError('No source hashes in evidence manifest')
    if not (ROOT/'v13/results/final_5s/evidence.json').is_file():
        raise AssertionError('Numerical run evidence missing')
    if 'NOT VERIFIED' not in (ROOT/'v13/openfoam/Allrun').read_text():
        raise AssertionError('OpenFOAM run must fail closed until implementation exists')
    workflows = list((ROOT/'.github/workflows').glob('*.yml'))
    if len(workflows) != 1:
        raise AssertionError(f'Expected only one active workflow, got {len(workflows)}')
    return {'result': 'PASS', 'active_workflows': len(workflows),
            'openfoam_runtime_verified': False, 'experimental_validated': False,
            'hashed_active_source_files': len(manifest['source_sha256'])}

if __name__ == '__main__':
    print(json.dumps(audit(), indent=2))
