"""Numerical regression tests; no experimental/OpenFOAM-runtime claims."""
from dataclasses import replace
import numpy as np
import pytest
from v13.gas_properties import cp_molar,enthalpy_molar,gas_energy,invert_energy
from v13.thermo_pressure import advance
from v13.em_feedback import solve_em_thermal_feedback
from v13.compressible_transport import R
from v13.reactor import Config,run


def test_shomate_sensible_heat_matches_derivative():
    for t in (900.,1020.,1080.,1150.,1270.,1350.):
        numeric=(enthalpy_molar(t+.01)-enthalpy_molar(t-.01))/.02
        np.testing.assert_allclose(numeric,cp_molar(t),atol=2e-3,rtol=1e-5)


def test_heat_capacity_positive_and_multicomponent():
    cp=cp_molar(np.linspace(950,1373,41))
    assert cp.shape==(5,41)
    assert np.all(cp>20)
    assert np.ptp(cp[:,20])>20 # CO2 vs hydrogen heat capacity differs greatly


def test_gas_energy_temperature_inversion():
    c=np.array([[2.,1.],[1.,3.],[10.,10.],[1.,2.],[.02,.05]])
    temp=np.array([1020.,1150.]); eps=.45
    e=gas_energy(c,temp,eps)
    assert np.all(e>0)
    np.testing.assert_allclose(invert_energy(e,c,eps,T_guess=np.array([1000.,1000.])),temp,atol=1e-8)


def test_shomate_continuity_at_piecewise_thresholds():
    for t in (500.,1000.,1200.,1300.):
        low=enthalpy_molar(t-.0001)
        hi=enthalpy_molar(t+.0001)
        assert np.max(np.abs(hi-low))<1. # J/mol


def test_pressure_and_temperature_are_consistent_same_time_step():
    n=9;eps=.45;temp=np.full(n,1080.)
    y=np.array([.0,.1,.8,.1,0.])
    g=np.repeat((101330/(R*1080)*y)[:,None],n,axis=1)
    E=gas_energy(g,temp,eps)
    s=np.zeros_like(g)
    eother=np.full(n,2000.)
    new,j,sj,info,en,tnew,hflux,diag=advance(g,temp,E,s,eother,dt=.025,
           length=.06,porosity=eps,pin=101330,pout=101329.5,inlet_fraction=y,T_in=1080.)
    assert diag['iterations']>1
    assert diag['pressure_EOS_residual_Pa']<1.1e-4
    np.testing.assert_allclose(info['pressure_Pa'],R*tnew*new.sum(axis=0),atol=1.1e-4)
    assert np.max(np.abs(en-gas_energy(new,tnew,eps)))<1e-4


def test_fixed_point_em_known_solution():
    e=np.full(4,10e6); cps=1e6;t0=e/cps;dt=1.
    # q(T) = 100 + 2 (T-10); solid RHS = 0.
    fn=lambda T:(100.+2*(T-10),{'annulus_EM_power_W':float(np.mean(100.+2*(T-10)))})
    q,report=solve_em_thermal_feedback(e,np.zeros(4),dt,cps,t0,fn)
    expected=100./(1.-2e-6)
    np.testing.assert_allclose(q,expected,rtol=1e-7)
    assert report['feedback_iterations']<10


def test_em_feedback_rejects_invalid_params():
    with pytest.raises(ValueError):
        solve_em_thermal_feedback(np.ones(2),np.ones(2),1.,1.,np.ones(2),lambda t:(t,{}),relaxation=0)


@pytest.fixture(scope='module')
def solution():
    return run(Config(duration_s=2.,em_nr=16,em_nz=24))


def test_full_reactor_conservation_and_picard(solution):
    assert solution['max_element_balance_mol_m2']<1e-9
    assert solution['relative_energy_ledger_error']<1e-10
    assert solution['max_end_step_EOS_residual_Pa']<1.1e-4
    assert solution['max_pressure_energy_Picard_iterations']<60
    assert solution['max_EM_feedback_relative_residual']<1e-7


def test_cold_or_invalid_temperature_rejected():
    with pytest.raises(ValueError):run(Config(T_initial_K=950.))


def test_zero_induction_no_em_iterations():
    r=run(Config(duration_s=.2,enable_em=False,enable_radicals=False))
    assert r['max_EM_feedback_iterations_per_refresh']==0
    assert r['max_element_balance_mol_m2']<1e-9


def test_gas_properties_change_with_temperature():
    assert abs(cp_molar(1000.)[1]-cp_molar(1350.)[1])>3.


def test_no_plasma_no_hydrogen_atom_creation():
    r=run(Config(duration_s=.2,enable_radicals=False,enable_em=False))
    assert r['H_radical_mean_mol_m3_pore']==0.
