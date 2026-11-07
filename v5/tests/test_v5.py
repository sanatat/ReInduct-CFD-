import json
import numpy as np
from v5.em2d import manufactured,run as induction
from v5.kinetics_validation import ROWS, fit, predict
from v5.thermochemistry import thermo, affinity_factor
from v4.battery_chemistry import atomic_residual,batch
from v4.coupled_reactor import run as reactor

def test_em_manufactured_second_order():
    a=manufactured(18,12)['relative_l2']
    b=manufactured(36,25)['relative_l2']
    c=manufactured(72,50)['relative_l2']
    assert a/b>3.5 and b/c>3.5 and c<1e-4

def test_em_sparse_linear_system_residual():
    s,_=induction(nx=41,ny=21,write=False)
    assert s['algebraic_relative_residual']<1e-9

def test_em_joule_power_is_nonnegative_and_mesh_stable():
    a,_=induction(nx=61,ny=30,write=False)
    b,f=induction(nx=81,ny=40,write=False)
    assert a['total_joule_power_W_per_m_out_of_plane']>0
    assert np.min(f['q'])>=0
    assert abs(a['total_joule_power_W_per_m_out_of_plane']/b['total_joule_power_W_per_m_out_of_plane']-1)<.015

def test_induction_changes_with_temperature_feedback():
    a,_=induction(nx=29,ny=15,write=False)
    b,_=induction(nx=29,ny=15,temperature=np.full((15,29),1050.),write=False)
    assert b['total_joule_power_W_per_m_out_of_plane'] != a['total_joule_power_W_per_m_out_of_plane']

def test_kinetics_exact_numerical_transcription():
    arr=np.asarray(ROWS)
    assert len(ROWS)==19
    assert np.isclose(arr[(arr[:,0]==800)&(arr[:,1]==1),2],2.80).any()
    assert np.isclose(arr[(arr[:,0]==1000)&(arr[:,1]==.26),2],1.34).any()

def test_kinetic_heldout_first_order_reasonable():
    T,p,y,sd=np.asarray(ROWS).T
    tr=p>=.53;te=p==.26
    coef=fit(T[tr],y[tr],sd[tr],p[tr],'segmented',1.)
    pred=predict(coef,T[te],p[te],'segmented',1.)
    assert np.mean(np.abs(pred-y[te]))<.12

def test_thermo_positive_K_and_reversible_affinity():
    for T in (550.,700.,900.):
        K=thermo(T)['Kp']
        assert K>1
        assert np.isclose(affinity_factor(T,1.,K),0,atol=1e-12)
        assert affinity_factor(T,1.,2*K)<0
        assert affinity_factor(T,1.,.01)>0

def test_battery_elements_and_mass_conserved():
    assert all(v==0 for v in atomic_residual().values())
    a=batch(feed_mol=10.,carbon_mol=7.5,conversion=.6)
    assert a['absolute_mass_balance_error_g']<1e-9

def test_em_reactor_two_way_and_conservation(tmp_path):
    r=reactor(nx=16,ny=8,dt=.002,end=.12,em_backend='2d',thermo_backend='tabulated',out=str(tmp_path))
    assert r['h2_balance_relative']<1e-9
    assert r['energy_balance_relative']<1e-9
    assert r['max_divergence_1_s']<1e-8
    assert r['last_joule_heat_source_mean_W_m3']>0
    assert r['max_h2_nonnegative']

def test_source_is_not_marked_experimentally_validated(tmp_path):
    r=reactor(nx=16,ny=8,dt=.002,end=.04,em_backend='2d',thermo_backend='tabulated',out=str(tmp_path))
    assert r['experimentally_validated'] is False
    assert r['openfoam_executed'] is False


def test_mhd_motional_force_dissipates_mechanical_energy():
    from v5.mhd2d import motional_lorentz
    u=np.array([[.02,.01],[.04,.015]])
    v=np.array([[.01,.03],[-.02,.03]])
    bx=np.full_like(u,.2+0.1j,dtype=complex)
    by=np.full_like(u,-.12+0.05j,dtype=complex)
    ax,ay,powvol=motional_lorentz(u,v,bx,by,3.,.1)
    assert (powvol<=1e-13).all()
    assert (powvol<0).any()


def test_em_plus_mhd_reactor_stable(tmp_path):
    r=reactor(nx=16,ny=8,dt=.002,end=.12,em_backend='2d',thermo_backend='tabulated',mhd_sigma_g=5.,out=str(tmp_path))
    assert r['cumulative_mhd_motional_mechanical_energy_J_per_m']<=0.
    assert r['max_divergence_1_s']<1e-8
    assert r['h2_balance_relative']<1e-9
    assert r['energy_balance_relative']<1e-8


def test_openfoam_heat_field_export_has_expected_cell_count():
    from pathlib import Path
    path=Path(__file__).resolve().parents[2]/'v4/openfoam_v2606/channelReactive2D/0/qEM'
    s=path.read_text()
    assert 'dimensions [1 -1 -3 0 0 0 0]' in s
    assert 'List<scalar>\n360\n(' in s
    source=(Path(__file__).resolve().parents[2]/'v4/openfoam_v2606/reactiveFixedBedFoam/reactiveFixedBedFoam.C').read_text()
    assert 'qEM/rhoCpS' in source
