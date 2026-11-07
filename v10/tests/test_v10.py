"""Version 10 regression evidence: no fabricated reactor measurements."""
from dataclasses import replace
from pathlib import Path
import json
import numpy as np
import pytest
from v10.plasma_transport_1d import Settings,run,electric_field,initialize,integrals
from v10.mhd_hartmann_heat import solve
from v10.lco_published_rates import EARLY,LATE,arrhenius_fit,held_out_rate_test,build
from v10.external_phases import audit
from v10.openfoam.postcheck import check


def test_plasma_atom_charge_and_electrostatic_energy_conservation():
    r,s,_=run()
    assert r['absolute_H_atom_inventory_error']<1e-10
    assert r['absolute_integrated_charge_error']<1e-10
    assert r['absolute_total_energy_including_Poisson_field_error']<1e-8
    assert all(np.all(np.isfinite(s[k])) for k in ('h2','h','ne','ni','ee','eg'))


def test_plasma_reaction_not_just_frozen_fields():
    cfg=Settings()
    _,s,_=run(cfg)
    start=initialize(cfg)
    assert s['h'].mean()>start['h'].mean()
    assert s['ni'].mean()>start['ni'].mean()


def test_plasma_explicit_stability_guard():
    with pytest.raises(ValueError,match='stability'):
        run(Settings(dt=.2,duration=.8))


def test_plasma_periodic_poisson_analytic_sine():
    c=Settings(n=128)
    x=(np.arange(c.n)+.5)*c.length/c.n
    ne=np.ones(c.n)*.01
    amplitude=.001
    ni=ne+amplitude*np.sin(2*np.pi*x/c.length)
    e,phi=electric_field(ni,ne,c)
    expected=c.e_external-c.poisson_coupling*amplitude/(2*np.pi/c.length)*np.cos(2*np.pi*x/c.length)
    assert np.max(np.abs(e-expected)) < 2e-7


def test_plasma_poisson_rejects_unbalanced_periodic_charge():
    c=Settings(n=32)
    with pytest.raises(ValueError,match='zero net charge'):
        electric_field(np.ones(c.n)*.02,np.ones(c.n)*.01,c)


def test_plasma_halved_timestep_convergence():
    r1, s1, _=run(Settings(dt=.0008,duration=.8))
    r2, s2, _=run(Settings(dt=.0004,duration=.8))
    r3, s3, _=run(Settings(dt=.0002,duration=.8))
    e12=np.linalg.norm(s1['h']-s2['h'])
    e23=np.linalg.norm(s2['h']-s3['h'])
    assert e23<.6*e12


def test_mhd_analytical_profile_and_power():
    r=solve(n=97)
    assert r['relative_L2_flow_error']<1e-4
    assert abs(r['mechanical_budget_relative_error'])<1e-10
    assert r['Joule_heating_W_m2']>0
    assert r['viscous_dissipation_W_m2']>0
    assert r['max_temperature_K']>300


def test_mhd_grid_convergence():
    errors=[solve(n=n)['relative_L2_flow_error'] for n in (25,49,97)]
    assert errors[0]>3*errors[1]>9*errors[2]


def test_mhd_zero_magnetic_field_poiseuille():
    r=solve(B=0,n=63)
    assert r['relative_L2_flow_error']<1e-10
    assert r['Joule_heating_W_m2']==0


def test_mhd_magnetic_braking_reduces_flow():
    r0=solve(B=0)
    r1=solve(B=3)
    assert r0['pressure_input_W_m2']>r1['pressure_input_W_m2']


def test_lco_published_data_stages_exactly_5_each():
    assert len(EARLY)==len(LATE)==5
    assert EARLY[1][2]==.00115 and LATE[3][2]==.002515


def test_lco_arrhenius_activation_energies_reproduce_published():
    assert abs(arrhenius_fit(EARLY)['E_a_J_per_mol']-121e3)<2e3
    assert abs(arrhenius_fit(LATE)['E_a_J_per_mol']-95e3)<2e3


def test_lco_rate_holdout_does_not_use_held_out_point():
    folds=held_out_rate_test(EARLY)['rows']
    assert len(folds)==5
    assert all(np.isfinite(f['log_error']) for f in folds)
    assert abs(folds[2]['relative_error'])>.1


def test_independent_study_is_not_mislabeled_as_numeric_conversion_validation(tmp_path):
    result=audit(tmp_path/'external.json')
    assert result['can_quantitatively_validate_conversion'] is False
    assert result['metallic_cobalt_vs_cobalt_oxide_mismatch'] is True
    assert result['lithium_carbonate_phase_compatible'] is True


def test_published_rate_provenance_serialized(tmp_path):
    result=build(tmp_path)
    table=(tmp_path/'published_rate_table.csv').read_text()
    assert '10.1007/s11663-022-02712-1' in table
    assert table.count('model-fit constant')==10
    assert result['external_reactor_conversion_validated'] is False


def test_openfoam_actual_runtime_fails_on_missing_log(tmp_path):
    with pytest.raises(ValueError,match='No solver log'):
        check(tmp_path)


def test_openfoam_actual_runtime_fails_on_fake_incomplete_log(tmp_path):
    (tmp_path/'log.solver').write_text('Time = 0.001\nmeanX=0.4\n')
    with pytest.raises(ValueError,match='End'):
        check(tmp_path)


def test_openfoam_pipeline_includes_compile_mesh_mapping_and_solver():
    script=Path('v10/openfoam/Allrun').read_text()
    for action in ('wmake','blockMesh','checkMesh','exportCellGeometry','v7.mapping','reactiveFixedBedFoamV7','v10.openfoam.postcheck'):
        assert action in script
    assert 'cp -a' in script
    assert 'WM_PROJECT_VERSION' in script


def test_external_reactor_validation_still_fails_closed():
    from v9.external_validation import evaluate
    with pytest.raises(ValueError,match='real measured points'):
        evaluate('v9/data/independent_reactor_measurements_TEMPLATE.csv')
