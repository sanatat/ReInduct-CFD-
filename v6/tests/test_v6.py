import math
import subprocess
from pathlib import Path
import numpy as np
import pytest
from v6.axisymmetric_em import Grid,manufactured,run as em_run,sample_heating
from v6.axisymmetric_bed import run as bed_run,geometry,ergun_velocity
from v6.kinetics_cv import crossvalidate
from v5.thermochemistry import thermo

@pytest.mark.parametrize('n',[12,18,28])
def test_axisymmetric_mms_accuracy(n):
    m=manufactured(Grid(nr=n,nz=round(1.6*n)))
    assert m['relative_L2']<.007
    assert m['algebraic_relative_residual']<1e-11

def test_axisymmetric_mms_converges_quadratically():
    e1=manufactured(Grid(nr=14,nz=22))['relative_L2']
    e2=manufactured(Grid(nr=28,nz=44))['relative_L2']
    assert e2/e1<.30

def test_em_nonnegative_dissipation_and_solver_residual():
    s,f=em_run(Grid(nr=18,nz=28))
    assert s['joule_power_W']>0
    assert np.min(f['q'])>=0
    assert s['algebraic_relative_residual']<1e-10

def test_power_mesh_stability():
    a,_=em_run(Grid(nr=20,nz=32));b,_=em_run(Grid(nr=36,nz=58))
    assert abs(a['joule_power_W']/b['joule_power_W']-1)<.05

def test_em_zero_coil_zero_power():
    a,f=em_run(Grid(nr=12,nz=20),coil_current_A_m2=0.)
    assert a['joule_power_W']==0
    assert np.max(np.abs(f['q']))==0

def test_em_temperature_feedback_changes_power():
    g=Grid(nr=14,nz=23)
    a,_=em_run(g)
    b,_=em_run(g,bed_temperature_K=np.full((g.nz,g.nr),950.))
    assert not np.isclose(a['joule_power_W'],b['joule_power_W'])

def test_em_transfer_q_nonnegative():
    a,f=em_run(Grid(nr=16,nz=26))
    q=sample_heating(f,np.linspace(.04,.10,12),np.linspace(.001,.025,8))
    assert q.shape==(12,8) and np.min(q)>=0

@pytest.mark.parametrize('on',[False,True])
def test_axisym_reactor_conservative(on):
    s=bed_run(nz=16,nr=7,end_s=.6,dt_s=.006,em_nr=20,em_nz=32,
              induction=on,write=False)
    assert abs(s['hydrogen_inventory_balance_residual_mol'])<1e-12
    assert abs(s['water_plus_solid_inventory_balance_residual_mol'])<1e-12
    assert abs(s['sensible_energy_ledger_residual_J'])<1e-7
    assert 0 <= s['mean_NiO_conversion'] <= 1
    assert s['min_H2_mol_m3']>=0

def test_no_reaction_has_no_conversion():
    s=bed_run(nz=12,nr=6,end_s=.3,dt_s=.005,induction=False,reaction=False,write=False)
    assert s['mean_NiO_conversion']==0

def test_em_changes_reactor_heat():
    kw=dict(nz=12,nr=6,end_s=.6,dt_s=.005,em_nr=18,em_nz=30,write=False)
    a=bed_run(induction=False,**kw)
    b=bed_run(induction=True,**kw)
    assert b['max_solid_temperature_K']>a['max_solid_temperature_K']

def test_ergun_increased_pressure_increases_flow():
    g=geometry(8,5);t=np.full((8,5),850.);x=np.zeros((8,5))
    a=ergun_velocity(2.,t,x,g)[0];b=ergun_velocity(8.,t,x,g)[0]
    assert b>a>0

def test_cv_all_four_concentration_groups_heldout():
    r=crossvalidate(write=False)
    assert set(x['heldout_H2_fraction'] for x in r['single_arrhenius']['folds'])=={.26,.53,.75,1.}
    assert all(x['pooled_MAE_um_s']>=0 and np.isfinite(x['pooled_MAE_um_s']) for x in
              (r['single_arrhenius'],r['two_regime_first_order'],r['two_regime_free_order']))

def test_cpp_reaction_matches_python(build_kinetics):
    for T,c,w,X in [(600,8.,0.,0.),(750,9.,2.,.2),(890,5.,8.,.5),(990,12.,.01,.9)]:
        K=float(thermo(float(T))['Kp'])
        result=subprocess.run([str(build_kinetics),str(T),str(c),str(w),str(X),str(K)],
                              text=True,capture_output=True,check=True)
        expected=.45*750.*math.exp(-48000/(8.314462618*T))*max(c-w/K,0)*(1-X)
        assert math.isclose(float(result.stdout),expected,rel_tol=1e-12,abs_tol=1e-12)

@pytest.fixture(scope='module')
def build_kinetics(tmp_path_factory):
    dest=tmp_path_factory.mktemp('v6cpp')/'kinetics'
    subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',
                    'v6/cpp/kinetics.cpp','-o',str(dest)],check=True,capture_output=True)
    return dest


def test_cpp_invalid_input_rejected(build_kinetics):
    bad=subprocess.run([str(build_kinetics),'300','5','0','0','12'],capture_output=True)
    assert bad.returncode != 0

def test_lco_stoichiometric_scenarios_conserve_mass():
    from v6.battery_flowsheet import simulate
    for conv in [0,.25,.5,1.]:
        for carbon in [0,50,200]:
            r=simulate(conversion=conv,extra_carbon_g=carbon)
            assert abs(r['mass_balance_residual_g'])<1e-9
            assert all(v==0 for v in r['elemental_balance_reaction'].values())
            assert r['theoretical_Co_formed_g']>=0

def test_battery_recovery_is_explicitly_scenario_not_measured():
    from v6.battery_flowsheet import simulate
    r=simulate(conversion=.8,recovery_separation_fraction=.25)
    assert np.isclose(r['assumed_collected_Co_g'],.25*r['theoretical_Co_formed_g'])
    assert r['process_validation'].startswith('none')

def test_openfoam_static_preflight_not_compile_claim():
    from v6.openfoam.preflight import inspect
    r=inspect()
    assert r['static_preflight_passed']
    assert r['qEM_cell_count_verified']
    assert not r['openfoam_compiled'] and not r['openfoam_solver_run']
