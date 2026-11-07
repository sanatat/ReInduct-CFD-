"""Numerical verification only; NO claims of experimental reactor validation."""
from dataclasses import replace
from pathlib import Path
import shutil
import numpy as np
import pytest
from v12.compressible_transport import pressure_corrected_step, R
from v12.reactor import ATOMS,S,Config,run
from v12.openfoam.case_lint import lint

@pytest.fixture(scope='module')
def base():
    return run(Config(duration_s=2.,dt_s=.025,em_nr=16,em_nz=24))

def test_all_reaction_pathways_element_balanced():
    np.testing.assert_allclose(ATOMS.T@S,np.zeros((6,4)),atol=1e-13)

def test_radical_hydrogen_conservation_and_activation(base):
    assert base['H_radical_mean_mol_m3_pore']>0
    assert abs(base['element_residual_by_species_mol_m2']['H'])<1e-9

def test_total_atomic_balance(base):
    assert base['max_element_balance_mol_m2'] < 1e-9

def test_energy_including_chemical_bond_and_boundary_flux(base):
    assert base['relative_energy_ledger_error'] < 1e-10

def test_implicit_pressure_closure(base):
    assert base['max_implicit_pressure_closure_relative'] < 1e-7
    assert base['min_pressure_Pa']>1e5
    assert base['max_advective_CFL']<.7

def test_pressure_flux_sign_and_molar_conservation():
    n=12; eps=.45;dt=.01;T=np.full(n,1080.)
    y=np.array([0.,.12,.85,.03,0.]);p=101330.
    g=np.repeat((y*p/(R*1080))[:,None],n,axis=1)
    ng,J,F,d=pressure_corrected_step(g,T,np.zeros_like(g),dt,.06,eps,p,p-.5,inlet_fraction=y)
    assert J[0]>0 and J[-1]>0
    before=g.sum()*eps*.06/n
    after=ng.sum()*eps*.06/n
    np.testing.assert_allclose(after-before,dt*(J[0]-J[-1]),rtol=0,atol=1e-11)
    assert d['pressure_closure_relative']<1e-9

def test_closed_column_at_uniform_pressure_stationary():
    n=8;eps=.45;dt=.01;T=np.full(n,1050.)
    y=np.array([0.,.1,.8,.1,0.]);p=101325.
    g=np.repeat((y*p/(R*1050.))[:,None],n,axis=1)
    ng,J,F,d=pressure_corrected_step(g,T,np.zeros_like(g),dt,.04,eps,p,p,inlet_fraction=y)
    np.testing.assert_allclose(ng,g,rtol=0,atol=1e-9)
    assert abs(J).max()<1e-7

def test_no_radical_power_no_new_h_radicals():
    cfg=Config(duration_s=1.,enable_em=False,plasma_electron_power_W_m3=0.)
    r=run(cfg)
    assert r['H_radical_mean_mol_m3_pore'] == pytest.approx(0,abs=1e-12)

def test_no_reaction_no_lco_conversion():
    cfg=Config(duration_s=1.,enable_em=False,enable_radicals=False,enable_solid_reactions=False)
    r=run(cfg)
    assert abs(r['mean_LCO_conversion'])<1e-12
    assert r['max_element_balance_mol_m2']<1e-9

def test_em_feedback_updates(base):
    assert len(base['electromagnetic_updates'])>=2
    assert all(x['relative_power_error']<1e-9 for x in base['electromagnetic_updates'])
    ps=[x['annulus_EM_power_W'] for x in base['electromagnetic_updates']]
    assert ps[0]!=ps[-1]

def test_initial_thermo_config_invalid():
    with pytest.raises(ValueError):run(Config(pressure_out_Pa=150000.))
    with pytest.raises(ValueError):run(Config(inlet_H2_fraction=.95,inlet_CO2_fraction=.12))

@pytest.mark.parametrize('nz', [12,24])
def test_mesh_conservative(nz):
    r=run(Config(nz=nz,duration_s=1.,enable_em=False))
    assert r['max_element_balance_mol_m2']<1e-9
    assert r['relative_energy_ledger_error']<1e-9

@pytest.mark.parametrize('dt', [.05,.025])
def test_temporal_conservative(dt):
    r=run(Config(dt_s=dt,duration_s=1.,enable_em=False))
    assert r['max_element_balance_mol_m2']<1e-9
    assert r['relative_energy_ledger_error']<1e-9

def test_openfoam_static_preflight_does_not_claim_runtime():
    result=lint()
    assert result['static_case_syntax_check_passed']
    assert not result['openfoam_build_verified']

def test_lint_catches_missing_divergence_scheme(tmp_path):
    original=Path(__file__).resolve().parents[2]
    source=original/'v11/openfoam'
    target=tmp_path/'v11/openfoam'
    shutil.copytree(source,target,ignore=shutil.ignore_patterns('*.log','__pycache__'))
    schemes=target/'case/system/fvSchemes'
    schemes.write_text(schemes.read_text().replace('div(phi,CO2) Gauss upwind;', ''))
    r=lint(tmp_path)
    assert not r['static_case_syntax_check_passed']
    assert any('CO2' in s for s in r['errors'])
