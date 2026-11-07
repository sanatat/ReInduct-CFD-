"""v8 model-specific unit, dimensional, and numerical regression checks."""
from pathlib import Path
import numpy as np
import pytest
from v8.physics import (layer_speed_to_bulk_rate,molar_density_ratio,h2_mole_fraction,
                        ideal_gas_mixture,volumetric_radiation_loss,ergun_pressure_gradient,
                        RHO_NI,MW_NI)
from v7.validation.interfacial_law import interfacial_to_bulk_rate
from v8.lco_kinetics import lco_conversion,ginstling_brounshtein,gb_inverse
from v8.lco_holdout import estimate
from v8.thermal_em import solve as em_thermal_solve

ROOT=Path(__file__).resolve().parents[2]


def test_ni_product_layer_not_oxide_molar_density():
    assert molar_density_ratio()==pytest.approx(1.731,rel=.002)
    expected=1e-6*(RHO_NI/MW_NI)*2
    assert layer_speed_to_bulk_rate(1.,2.)==pytest.approx(expected,rel=1e-12)


def test_product_layer_porosity_reduces_rate():
    dense=layer_speed_to_bulk_rate(1.,2.)
    porous=layer_speed_to_bulk_rate(1.,2.,ni_layer_porosity=.25)
    assert porous==pytest.approx(.75*dense)
    with pytest.raises(ValueError):layer_speed_to_bulk_rate(1.,2.,ni_layer_porosity=1)


def test_tgas_not_tsolid_in_h2_fraction_and_kinetic_source():
    c=6.
    r1=interfacial_to_bulk_rate(1100.,c,.1,Tgas_K=800.)
    r2=interfacial_to_bulk_rate(1100.,c,.1,Tgas_K=1000.)
    assert r2/r1==pytest.approx(1.25)
    assert h2_mole_fraction(c,800.,101325.)<h2_mole_fraction(c,1000.,101325.)


def test_kinetics_mass_mapping_matches_analytical():
    a=interfacial_to_bulk_rate(1100.,7.,.2,area_m2_m3=3.,ni_product_porosity=.15)
    b=interfacial_to_bulk_rate(1100.,7.,.2,area_m2_m3=3.,ni_product_porosity=0)
    assert a==pytest.approx(.85*b)


def test_ideal_gas_mixture_closes_fractions_and_dilutes_density():
    a=ideal_gas_mixture(.6,.2,850,101325)
    b=ideal_gas_mixture(.6,.2,1000,101325)
    assert a['yH2']+a['yH2O']+a['yN2']==pytest.approx(1.)
    assert a['rho_kg_m3']>b['rho_kg_m3']>0
    assert a['cp_mass_J_kg_K']>0
    assert a['volumetric_heat_capacity_J_m3_K']>0


def test_invalid_gas_mixture_rejected():
    with pytest.raises(ValueError):ideal_gas_mixture(.9,.3,850,101325)
    with pytest.raises(ValueError):h2_mole_fraction(100,1000,101325)


def test_radiation_sign_and_t4_law():
    q=volumetric_radiation_loss(950,700,.7,80)
    assert q>0
    assert volumetric_radiation_loss(700,950,.7,80)==pytest.approx(-q)
    assert volumetric_radiation_loss(850,850,.7,80)==0


def test_ergun_pressure_drop_increases_with_velocity():
    a=ergun_pressure_gradient(.01,.45,.003,3e-5,.2)
    b=ergun_pressure_gradient(.02,.45,.003,3e-5,.2)
    assert b>a>0


def test_bed_conductivity_em_feedback_converges_and_changes_power():
    result,hist=em_thermal_solve(nr=18,nz=30)
    assert result['converged'] and hist[-1]['max_temperature_change_K']<2e-5
    assert .1<abs(result['power_feedback_relative_change_pct'])<30
    assert 700<result['peak_temperature_K']<1100


def test_radiation_affects_temperature():
    hot,_=em_thermal_solve(wall_emissivity=0)
    cool,_=em_thermal_solve(wall_emissivity=.7)
    assert cool['peak_temperature_K']<hot['peak_temperature_K']


def test_gb_kinetic_integral_inverse():
    for alpha in (.1,.3,.6,.95):
        assert gb_inverse(ginstling_brounshtein(alpha))==pytest.approx(alpha,abs=1e-9)


def test_lco_model_monotonic_time_and_temperature():
    a=lco_conversion(3000,973.15)
    b=lco_conversion(6000,973.15)
    c=lco_conversion(6000,1073.15)
    assert 0<=a<=b<=c<=1
    with pytest.raises(ValueError):lco_conversion(6000,1500)


def test_lco_holdout_never_trains_on_validation_temperatures():
    r=estimate()
    t={x['T_C'] for x in r['train']['rows']}
    v={x['T_C'] for x in r['held_out']['rows']}
    assert t.isdisjoint(v)
    assert np.isfinite(r['held_out']['MAE_conversion_fraction'])
    assert r['caveat']


def test_openfoam_candidate_has_corrected_physics():
    cpp=(ROOT/'v7/openfoam/solver/reactiveFixedBedFoamV7.C').read_text()
    assert 'rhoNi_mol_m3' in cpp
    assert 'rhoNiO_mol_m3' not in cpp
    assert 'max(Tg[i],1.0)' in cpp
    assert 'qRad' in cpp
    assert 'phiHbyA' in cpp
    assert 'maxDivergencePhi=' in cpp


def test_openfoam_kinetic_table_molar_density():
    tbl=(ROOT/'v7/openfoam/case/constant/kineticTable').read_text()
    assert 'NiProductLayerMolarDensity' in tbl and 'NiProductLayerPorosity' in tbl
    assert 'NiOMolarDensity' not in tbl
