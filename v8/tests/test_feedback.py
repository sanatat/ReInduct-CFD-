import csv
from pathlib import Path
import numpy as np
import pytest
from v7.mapping import synthetic_annular_wedge
from v8.openfoam.feedback import read_foam_internal_scalar,remap_from_temperature

def test_read_uniform_foam_field(tmp_path):
    f=tmp_path/'Ts';f.write_text('FoamFile { format ascii; }\ninternalField uniform 850;\nboundaryField {}')
    assert np.array_equal(read_foam_internal_scalar(f,3),np.array([850.,850.,850.]))

def test_read_nonuniform_foam_field(tmp_path):
    f=tmp_path/'Ts';f.write_text('FoamFile { format ascii; }\ninternalField nonuniform List<scalar>\n3\n(\n850\n860\n870\n);\nboundaryField {}')
    assert list(read_foam_internal_scalar(f,3))==[850,860,870]
    with pytest.raises(ValueError):read_foam_internal_scalar(f,4)

def test_feedback_from_synthetic_openfoam_like_fields(tmp_path):
    C,V=synthetic_annular_wedge(nr=8,nz=10)
    csvpath=tmp_path/'cells.csv'
    with csvpath.open('w',newline='') as fp:
        w=csv.writer(fp);w.writerow(['x_m','y_m','z_m','volume_m3'])
        for c,v in zip(C,V):w.writerow([*c,v])
    f1=tmp_path/'Ts850';f1.write_text('FoamFile { format ascii; }\ninternalField uniform 850;')
    f2=tmp_path/'Ts950';f2.write_text('FoamFile { format ascii; }\ninternalField uniform 950;')
    a=remap_from_temperature(csvpath,f1,tmp_path/'a',nr=14,nz=22)
    b=remap_from_temperature(csvpath,f2,tmp_path/'b',nr=14,nz=22)
    assert abs(a['post_mapping_relative_power_error'])<1e-12
    assert abs(b['post_mapping_relative_power_error'])<1e-12
    assert abs(a['final_mapped_power_W']/b['final_mapped_power_W']-1)>.001
    assert (tmp_path/'b/qEM').exists()
