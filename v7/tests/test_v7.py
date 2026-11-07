import json
from pathlib import Path
import numpy as np
import pytest
from v7 import mapping
from v7.validation.nio_holdout import read_data,run,fit_global_arrhenius
from v6.axisymmetric_em import Grid,run as em_run

def test_published_data_19_values():
    t=read_data();assert t.shape==(19,4)
    assert float(t[(t[:,0]==800)&(t[:,1]==1)][0,2])==2.8
    assert float(t[(t[:,0]==675)&(t[:,1]==1)][0,2])==4.04
    assert np.all(t[:,2]>0) and np.all(t[:,3]>0)

def test_kinetic_predictions_nonnegative():
    d=read_data();m=fit_global_arrhenius(d);assert np.all(m(d[:,0],d[:,1])>0)

def test_paper_holdout_split_is_disjoint(tmp_path):
    report=run(tmp_path)
    assert report['holdout_800C']['global_arrhenius']['n']==4
    assert report['leave_H2_26pct_out']['global_arrhenius']['n']==4
    assert len(report['leave_one_temperature_out_interior'])==4

def test_annular_sector_volume_matches_analytical():
    c,V=mapping.synthetic_annular_wedge()
    analytical=np.deg2rad(4)*.5*(.025**2-.005**2)*(.100-.040)
    assert sum(V)==pytest.approx(analytical,rel=1e-12)
    assert c.shape==(576,3)

def test_map_conserves_power_and_source_units(tmp_path):
    _,f=em_run(Grid(nr=26,nz=42))
    c,V=mapping.synthetic_annular_wedge()
    pref=mapping.em_ref_power(f,.005,.025,.04,.10,nr=80,nz=100)
    q,rep=mapping.map_conservatively(f,c,V,pref)
    assert np.all(q>=0) and np.isfinite(q).all()
    assert rep['post_mapping_relative_power_error']==pytest.approx(0.,abs=1e-12)
    p=tmp_path/'qEM';mapping.write_foam_scalar(p,q)
    txt=p.read_text();assert '[1 -1 -3 0 0 0 0]' in txt
    assert '576' in txt and 'type wedge' in txt

def test_mesh_load_rejects_bad_csv(tmp_path):
    f=tmp_path/'bad.csv';f.write_text('x_m,y_m,z_m,volume_m3\n0,0,0,-1\n')
    with pytest.raises(ValueError):mapping.load_actual_cells(f)

def test_source_uses_molar_units_not_dimensionless():
    case=Path(__file__).resolve().parents[1]/'openfoam/case/0'
    for k in ['H2','H2O']:
        assert 'dimensions [0 -3 0 0 1 0 0];' in (case/k).read_text()
    assert 'dimensions [1 -1 -3 0 0 0 0];' in (case/'qEM.template').read_text()

from v7.validation.closed_reactor import solve as closed_solve

def test_reactor_closed_conservation():
    result,_=closed_solve()
    assert abs(result['H2_inventory_error_mol_m3_bulk'])<1e-10
    assert abs(result['H2O_inventory_error_mol_m3_bulk'])<1e-10
    assert abs(result['energy_error_J_m3'])<5e-3

def test_energy_coupling_changes_temperature():
    a,_=closed_solve(q_ref_W_m3=18000)
    b,_=closed_solve(q_ref_W_m3=0)
    assert a['max_Ts_K']>b['max_Ts_K']
    assert a['final_X']>=b['final_X']

def test_time_step_consistency_in_conversion():
    a,_=closed_solve(dt=0.004);b,_=closed_solve(dt=0.002)
    assert abs(a['final_X']-b['final_X'])/max(b['final_X'],1e-20)<.01

from v7.validation.interfacial_law import (
    rate_table_from_all_data,interfacial_to_bulk_rate,export_foam_table
)
def test_interface_rate_units_and_conversion():
    k=rate_table_from_all_data()
    assert k.shape==(6,2)
    # For equal T and c, larger reactive area doubles volumetric rate.
    a=interfacial_to_bulk_rate(1100.,7.,.25,area_m2_m3=1)
    b=interfacial_to_bulk_rate(1100.,7.,.25,area_m2_m3=2)
    assert a>0 and b==pytest.approx(2*a)
    assert interfacial_to_bulk_rate(1100.,0.,.25)==0
    assert interfacial_to_bulk_rate(1100.,7.,1.)==0

def test_no_out_of_range_extrapolation():
    with pytest.raises(ValueError):interfacial_to_bulk_rate(750.,10.,0)

def test_exported_foam_table_has_all_temperature_nodes(tmp_path):
    path=tmp_path/'kineticTable';out=tmp_path/'law.json'
    record=export_foam_table(path,out)
    assert len(record['table_T_C'])==6
    assert 'rateUmPerSAtYH2One' in path.read_text()

from v7.validation.battery_phase_envelope import LOW,HIGH,scenario
from v4.battery_chemistry import atomic_residual

def test_battery_reactions_atom_and_mass_balance():
    for reaction in (LOW,HIGH):
        assert all(x==0 for x in atomic_residual(reaction).values())
    for t in (700,850,900,1000,1050,1100):
        s=scenario(t)
        assert abs(s['mass_balance_error_g'])<1e-9

def test_battery_phase_envelope_does_not_extrapolate():
    with pytest.raises(ValueError):scenario(950.)
    with pytest.raises(ValueError):scenario(1200.)
