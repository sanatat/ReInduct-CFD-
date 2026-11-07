"""New physically meaningful regression checks for v9, distinct from 76 older tests."""
from dataclasses import replace
from pathlib import Path
import json
import pytest
from v9.plasma_energy import PlasmaCase,run,DISSOCIATION_J_PER_MOL_H2
from v9.lco_process import scenario
from v9.validation_protocol import leave_one_temperature_out
from v9.openfoam.preflight import report


def test_radical_atoms_and_energy_balances():
    r,_=run()
    assert r['atomic_H_balance_abs_max_mol_m3']<1e-11
    assert r['total_energy_balance_abs_max_J_m3']<1e-7
    assert r['final_H_mol_m3']>0


def test_no_plasma_input_no_dissociation_without_energy():
    r,_=run(PlasmaCase(electron_energy_J_m3=0.,electron_heating_W_m3=0.,duration_s=1.))
    assert r['final_H_mol_m3']==0.
    assert r['final_H2_mol_m3']==1.


def test_recombination_gas_energy_conserved():
    case=PlasmaCase(h2_mol_m3=0.,h_mol_m3=1.,electron_energy_J_m3=0.,
                    electron_heating_W_m3=0.,duration_s=2.)
    r,_=run(case)
    assert r['final_H_mol_m3']<1
    assert r['final_neutral_energy_J_m3']>0
    assert r['atomic_H_balance_abs_max_mol_m3']<1e-11
    assert r['total_energy_balance_abs_max_J_m3']<1e-7


def test_halving_dt_changes_plasma_solution_small():
    a,_=run(PlasmaCase(dt_s=.01))
    b,_=run(PlasmaCase(dt_s=.005))
    assert abs(a['final_H_mol_m3']-b['final_H_mol_m3'])<5e-4


def test_li_mass_balance_low_temp():
    s=scenario(800,240,10,7.5)
    assert abs(s['mass_residual_g'])<1e-8
    assert all(v==0 for v in s['elemental_balance'].values())
    assert s['inventories_mol']['Co']>0
    assert s['inventories_mol']['Li2CO3']>0


def test_li_mass_balance_high_temp_and_carbon_limit():
    s=scenario(1050,240,10,2.0)
    assert abs(s['mass_residual_g'])<1e-8
    assert all(v==0 for v in s['elemental_balance'].values())
    assert s['inventories_mol']['C']>=-1e-12
    assert s['carbon_limited_conversion_LCO'] <= (2./15.) + 1e-12
    assert 'Li2O' in s['inventories_mol']


def test_li_time_monotonic_until_carbon_limit():
    y=[scenario(800,t,10,7.5)['carbon_limited_conversion_LCO'] for t in (0,30,60,120,240)]
    assert y[0]==0
    assert all(b>=a for a,b in zip(y,y[1:]))


def test_phase_gap_not_interpolated():
    with pytest.raises(ValueError):scenario(950,240)


def test_holdout_report_requires_independent_validation_flag_false():
    r=leave_one_temperature_out()
    assert len(r['cases'])==4
    assert all(t not in row['training_temperature_C'] for row in r['cases'] for t in [row['held_out_temperature_C']])
    assert r['independent_reactor_validation'] is False
    assert r['summary_MAE_conversion']>=0


def test_openfoam_preflight_cannot_mark_run_passed(tmp_path):
    p=tmp_path/'preflight.json';r=report(p)
    assert p.is_file() and not r['solver_compiled_and_executed']


def test_existing_openfoam_runner_fail_closed():
    text=Path('v9/openfoam/Allrun').read_text()
    assert 'wmake' in text and 'blockMesh' in text and 'checkMesh' in text
    assert 'v9/results/openfoam_runtime' in text
    assert 'VERIFIED.txt' in text


def test_external_validation_rejects_empty_template():
    from v9.external_validation import evaluate
    with pytest.raises(ValueError,match='real measured points'):
        evaluate('v9/data/independent_reactor_measurements_TEMPLATE.csv')


def test_external_validation_rejects_same_source_or_mismatched_chemistry(tmp_path):
    from v9.external_validation import evaluate
    p=tmp_path/'observations.csv'
    import csv
    cols=['source_doi','experiment_id','reactor_id','chemistry','temperature_C','time_min','carbon_mol_per_lco_mol','observed_conversion','measurement_uncertainty_abs']
    rows=[]
    for i in range(3):
        rows.append(['10.1007/s11663-022-02712-1',str(i),'r','LiCoO2+C','800','120','0.75','0.5','0.05'])
    with p.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(cols);w.writerows(rows)
    with pytest.raises(ValueError,match='not independent'):
        evaluate(p)
    for r in rows:r[0]='10.1234/OTHER';r[3]='NMC+H2'
    with p.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(cols);w.writerows(rows)
    with pytest.raises(ValueError,match='Mismatched chemistry'):
        evaluate(p)
