"""V11 numerical checks. Tests do not constitute experimental reactor validation."""
import numpy as np
from dataclasses import replace
import pytest
from v11.integrated_reactor import (Config,ATOMS,STOICH,kinetic_rates_per_s,
 calculate_em_heat,simulate,ergun_superficial)

@pytest.fixture(scope='module')
def base():
    return simulate(Config(duration_s=80.,dt_s=.1,nz=20,em_grid_nr=16,em_grid_nz=24,em_updates=3))

def test_reaction_element_stoichiometry():
    np.testing.assert_allclose(ATOMS.T@STOICH,0.,atol=1e-14)

def test_stage_rates_positive_increase_with_temperature():
    rates=kinetic_rates_per_s(np.array([1000.,1100.,1200.]))
    assert rates.shape==(2,3) and (rates>0).all()
    assert (np.diff(rates,axis=1)>0).all()

def test_no_extrapolation():
    with pytest.raises(ValueError,match='outside'):
        kinetic_rates_per_s(np.array([972.]))

def test_em_remap_conserves_power():
    cfg=Config(nz=18,em_grid_nr=16,em_grid_nz=24)
    _,r=calculate_em_heat(np.full(cfg.nz,1080.),cfg)
    assert r['relative_power_error']<1e-12

def test_em_temperature_feedback():
    cfg=Config(nz=16,em_grid_nr=16,em_grid_nz=24)
    _,a=calculate_em_heat(np.full(cfg.nz,1080.),cfg)
    _,b=calculate_em_heat(np.full(cfg.nz,1200.),cfg)
    assert abs(a['annulus_EM_power_W']-b['annulus_EM_power_W'])/a['annulus_EM_power_W']>.03

def test_conservative_elements_and_energy(base):
    assert base['max_abs_element_residual']<1e-8
    assert base['relative_energy_residual']<1e-10

def test_nonzero_battery_reaction_and_products(base):
    assert 0<base['mean_conversion_LCO']<1
    assert base['outlet_CO_mol_m3_pore']>=0
    assert base['mean_lithium_carbonate_mol_m3_bulk']>0

def test_feedback_em_updates(base):
    r=base['em_updates']
    assert len(r)==3
    assert r[0]['annulus_EM_power_W']!=r[-1]['annulus_EM_power_W']
    assert all(x['relative_power_error']<1e-12 for x in r)

def test_off_reactions_zero_conversion():
    cfg=Config(duration_s=30.,dt_s=.1,nz=16,enable_em=False,enable_chemistry=False)
    r=simulate(cfg)
    assert abs(r['mean_conversion_LCO'])<1e-13
    assert r['max_abs_element_residual']<1e-9

def test_off_heat_and_chemistry_energy_conservation():
    cfg=Config(duration_s=30.,dt_s=.1,nz=16,enable_em=False,
               enable_chemistry=False,enable_radiation=False)
    r=simulate(cfg)
    assert r['relative_energy_residual']<1e-10

def test_carbon_inventory_limit():
    cfg=Config(duration_s=40.,dt_s=.1,nz=16,enable_em=False,
               carbon_LCO_ratio=.00001)
    r=simulate(cfg)
    assert r['mean_conversion_LCO']>=0
    assert r['max_abs_element_residual']<1e-8

def test_ergun_u_decreases_as_pressure_drop_decreases():
    c=Config()
    gas=np.broadcast_to(np.array([0.,1.,9.])[:,None],(3,12))
    u=ergun_superficial(c,np.full(12,1080.),gas)
    u2=ergun_superficial(replace(c,pressure_gradient_Pa_m=4.),np.full(12,1080.),gas)
    assert 0<u2<u

def test_invalid_input_rejected():
    with pytest.raises(ValueError):simulate(Config(porosity=1.2))
