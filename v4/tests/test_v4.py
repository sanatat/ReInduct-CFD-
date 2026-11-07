"""Tests for new v4 features; legacy v3 tests remain separately accessible."""
import json
import math
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from v4.battery_chemistry import batch, atomic_residual
from v4.literature_benchmark import fit_and_score
from v4.coupled_reactor import run

def test_lco_atomic_balance():
    assert all(v==0 for v in atomic_residual().values())

def test_lco_mass_balance():
    for c in [0,.1,.6,1.]:
        b=batch(conversion=c)
        assert b['absolute_mass_balance_error_g']<1e-9

def test_lco_limited_carbon():
    b=batch(feed_mol=10.,carbon_mol=1.,conversion=1.)
    assert abs(b['molar_amounts']['C'])<1e-12
    assert b['molar_amounts']['LiCoO2']>0

def test_cxx_rate_matches_paper_function(tmp_path):
    exe=tmp_path/'nio_kinetics'
    subprocess.run(['g++','-std=c++17','-Wall','-Wextra','-O2',str(ROOT/'v4/cpp/nio_kinetics.cpp'),'-o',str(exe)],check=True)
    for T,a in [(650,.1),(800,.2),(950,.7)]:
        actual=float(subprocess.check_output([str(exe),str(T),str(a)],text=True))
        expected=math.exp(19.5-90800/(8.314462618*T))/60*a**.63*(1-a)**1.39
        assert abs(actual-expected)<1e-12

def test_literature_contains_holdout(tmp_path):
    r=fit_and_score(str(tmp_path))
    assert len(r['heldout_temperatures_C'])==4
    assert r['heldout_high_T_MAE_um_s']>0
    assert math.isfinite(r['derived_Ea_J_mol'])

def test_mass_energy_divergence(tmp_path):
    r=run(nx=20,ny=8,dt=.002,end=.5,out=str(tmp_path/'case'))
    assert r['h2_balance_relative']<1e-8
    assert r['energy_balance_relative']<1e-8
    assert r['max_divergence_1_s']<1e-8
    assert r['mean_conversion']>0
    assert r['mean_velocity_m_s']>0

def test_reaction_off(tmp_path):
    r=run(nx=20,ny=8,dt=.002,end=.35,reaction=False,em=False,out=str(tmp_path/'case'))
    assert abs(r['mean_conversion'])<1e-15
    assert r['total_consumed_H2_mol_per_m']==0

def test_coupling_changes_solution(tmp_path):
    r1=run(nx=20,ny=8,dt=.002,end=.8,back_coupling=True,out=str(tmp_path/'coupled'),beta_drag=.04)
    r2=run(nx=20,ny=8,dt=.002,end=.8,back_coupling=False,out=str(tmp_path/'fixed'),beta_drag=.04)
    assert abs(r1['mean_velocity_m_s']-r2['mean_velocity_m_s'])>1e-7
    assert abs(r1['mean_conversion']-r2['mean_conversion'])>1e-10

def test_openfoam_candidate_case_structure():
    # Static-only packaging check: this is NOT a compile or runtime test.
    root=ROOT/'v4/openfoam_v2606'
    assert (root/'reactiveFixedBedFoam/reactiveFixedBedFoam.C').exists()
    for f in ['U','p','H2','H2O','Tg','Ts','X']:
        assert (root/'channelReactive2D/0'/f).exists()
    for f in ['blockMeshDict','fvSchemes','fvSolution','controlDict']:
        assert (root/'channelReactive2D/system'/f).exists()
    assert 'UNCOMPILED' in (root/'README.md').read_text()


def test_high_temperature_arrhenius_failure_is_visible(tmp_path):
    r=fit_and_score(tmp_path)
    assert r['predictions_um_s'][3]>r['observations_um_s'][3]*2
    assert r['heldout_high_T_MAE_um_s']>r['heldout_600C_absolute_error_um_s']
