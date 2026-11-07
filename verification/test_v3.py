"""Fast physics-oriented verification tests; all models use synthetic data."""
import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from em.induction1d import solve as em_legacy
from em.variable_induction import solve_variable,conductivity_from_temperature
from reactor2d.thermo import equilibrium_constant,reaction_rate
from reactor2d.reactor import simulate,ergun_pressure_drop
from cfd.mac_channel import simulate as channel


def test_em_variable_reduces_to_constant_sigma():
    reference=em_legacy(n=301)
    new=solve_variable(np.full(301,1e6))
    assert np.max(abs(reference['H']-new['H']))/10000<1e-9
    assert np.all(new['q']>=0)


def test_em_grid_refinement_to_analytic():
    def error(n):
        r=solve_variable(np.full(n,1e6))
        skin=np.sqrt(2/(2*np.pi*1e5*4*np.pi*1e-7*1e6))
        k=(1+1j)/skin
        exact=1e4*np.cosh(k*(0.02-r['y']))/np.cosh(k*.02)
        return np.max(abs(r['H']-exact))/1e4
    assert error(401)<error(101)
    assert error(801)<error(401)
    assert error(801)<2e-4


def test_em_power_quadratic_field():
    a=solve_variable(np.full(101,1e6),h0=4e3)
    b=solve_variable(np.full(101,1e6),h0=8e3)
    assert np.allclose(b['q'],4*a['q'],rtol=1e-10,atol=1e-9)
    sig0=conductivity_from_temperature(np.array([850,900,1000]))
    assert sig0[0]>sig0[1]>sig0[2]


def test_thermodynamic_limit_and_kinetics():
    T=np.array([850.,900.])
    K=equilibrium_constant(T)
    assert np.all(K>0)
    rates=reaction_rate(T,np.full(2,15.),np.zeros(2),np.zeros(2))
    assert rates[1]>rates[0]>0
    assert np.all(reaction_rate(T,np.ones(2),2*K,np.zeros(2))==0)
    assert np.all(reaction_rate(T,np.ones(2),np.zeros(2),np.ones(2))==0)
    assert ergun_pressure_drop(0.02)>ergun_pressure_drop(0.01)


def test_channel_poisseuille_and_mass_projection(tmp_path):
    s=channel(nx=24,ny=12,end=4.,out=tmp_path/'ns')
    assert s['max_divergence_1_s']<1e-10
    assert s['relative_l2_vs_exact_channel_profile']<.04
    assert s['max_transverse_velocity_m_s']<1e-8


def test_channel_hartmann_brakes_and_matches_analytic(tmp_path):
    s=channel(nx=24,ny=12,end=5.,B_T=.1,out=tmp_path/'mhd')
    assert s['max_divergence_1_s']<1e-10
    assert s['relative_l2_vs_exact_channel_profile']<.04
    assert s['Hartmann_number_half_channel']>1
    assert s['mean_U_m_s']<.016


def test_nonuniform_brinkman_creates_transverse_flow(tmp_path):
    s=channel(nx=24,ny=12,end=3.,porous_drag=.7,out=tmp_path/'porous')
    assert s['max_divergence_1_s']<1e-10
    assert s['max_transverse_velocity_m_s']>1e-6


def test_reactor_zero_sources_exact_thermal_state(tmp_path):
    s=simulate(nx=20,ny=8,end=2.,no_em=True,no_reaction=True,
               emissivity=0.,out=tmp_path/'zero')
    assert s['mean_conversion']==0
    assert abs(s['max_solid_T_K']-850.)<1e-8
    assert s['hydrogen_balance_rel']<1e-10
    assert s['energy_balance_abs_J_per_m']<1e-6


def test_reactor_coupled_conservation_and_feedback(tmp_path):
    s=simulate(nx=24,ny=8,end=2.,out=tmp_path/'react')
    assert 0<s['mean_conversion']<1
    assert s['max_solid_T_K']>850
    assert s['hydrogen_balance_rel']<1e-10
    assert s['energy_balance_rel']<1e-10
    assert s['em_field_recomputations']>=2


def test_cfd_velocity_to_reactor(tmp_path):
    channel(nx=24,ny=8,end=2.5,porous_drag=.7,out=tmp_path/'flow')
    s=simulate(nx=24,ny=8,end=2.,velocity_file=tmp_path/'flow/velocity.npz',out=tmp_path/'coupled')
    assert s['velocity_model']=='computed_MAC_velocity'
    assert s['max_transverse_velocity_m_s']>0
    assert s['hydrogen_balance_rel']<1e-10
    assert s['energy_balance_rel']<1e-10


def test_induction_strength_and_conversion(tmp_path):
    a=simulate(nx=18,ny=6,end=2.,h0=3000.,out=tmp_path/'weak')
    b=simulate(nx=18,ny=6,end=2.,h0=10000.,out=tmp_path/'strong')
    assert b['max_solid_T_K']>a['max_solid_T_K']
    assert b['mean_conversion']>a['mean_conversion']
    assert b['cumulative_EM_heat_J_per_m']>a['cumulative_EM_heat_J_per_m']


def test_radical_plasma_surrogate_conserves_H_and_solid_O():
    from plasma.radical0d import solve
    r,s=solve()
    assert s['H_atoms_balance_abs_mol_m3']<1e-7
    assert s['solid_O_to_H2O_balance_abs_mol_m3']<1e-7
    assert s['radical_max_mol_m3']>0
    assert 0<s['NiO_conversion_final']<1
    _,without_diss=solve(k_diss=0.)
    assert without_diss['NiO_conversion_final']==0


def test_separate_region_conjugate_heat_conservation():
    from multiregion.conjugate1d import solve
    for Ng,Ns in ((8,10),(23,35)):
        x,t,s=solve(Ng=Ng,Ns=Ns)
        assert s['relative_flux_error']<1e-11
        assert s['max_face_heat_flux_discrepancy_W_m2']<1e-7
        assert 800<t.min()<t.max()<900
    _,_,contact=solve(Rcontact=.01)
    _,_,plain=solve()
    assert contact['analytic_heat_flux_W_m2']<plain['analytic_heat_flux_W_m2']
