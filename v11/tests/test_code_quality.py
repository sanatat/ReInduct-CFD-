"""Static integration checks; NOT substitutes for wmake / OpenFOAM runtime."""
from pathlib import Path
import json,subprocess,shutil
import numpy as np
from v11.integrated_reactor import kinetic_rates_per_s
ROOT=Path(__file__).resolve().parents[2]

def test_cpp_kernel_build_and_matches_python(tmp_path):
    assert shutil.which('g++')
    b=tmp_path/'check'
    subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-pedantic',
                    str(ROOT/'v11/cpp/kinetic_stoichiometry.cpp'),'-o',str(b)],check=True)
    r=json.loads(subprocess.check_output([str(b)],text=True))
    k=kinetic_rates_per_s(np.array([r['temperature_K']]))[:,0]
    np.testing.assert_allclose(k,[r['k1_per_second'],r['k2_per_second']],rtol=2e-11)
    assert r['all_reactions_atom_balanced']

def test_openfoam_candidate_has_new_chemistry_and_pimple():
    s=(ROOT/'v11/openfoam/solver/reactiveLCOFoamV11.C').read_text()
    for token in ('pimple.loop()', 'pimple.correct()', 'pimple.correctNonOrthogonal()',
                  'LCO','Carbon','Li2O','CoO','Li2CO3','CO2','qEM','r0','Carbon[celli]'):
        assert token in s
    assert 'reactiveFixedBedFoamV7' not in s
    assert 'H2O' not in s
    assert (ROOT/'v11/openfoam/case/constant/kineticsV11').is_file()
    assert 'reactiveLCOFoamV11' in (ROOT/'v11/openfoam/case/system/controlDict').read_text()

def test_openfoam_runtime_not_falsely_verified():
    assert not (ROOT/'v11/results/openfoam_runtime/VERIFIED.txt').exists()
