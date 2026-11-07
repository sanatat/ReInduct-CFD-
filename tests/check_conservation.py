#!/usr/bin/env python3
"""Run substantive invariants and limiting cases with the compiled C++ benchmark."""
import subprocess
import tempfile
from pathlib import Path
import numpy as np

exe=Path('build/reactor1d').resolve()
assert exe.is_file(), 'Build first: g++ -std=c++17 -O2 src/reactor1d.cpp -o build/reactor1d'
with tempfile.TemporaryDirectory() as d:
    root=Path(d)
    def sim(tag, *flags):
        dest=root/tag
        cmd=[str(exe), '--out',str(dest),*flags]
        subprocess.run(cmd,stdout=subprocess.DEVNULL,check=True)
        f=np.genfromtxt(dest/'fields.csv',delimiter=',',names=True)
        h=np.genfromtxt(dest/'history.csv',delimiter=',',names=True)
        assert np.all(f['solid_conversion']>=-1e-12)
        assert np.all(f['solid_conversion']<=1+1e-12)
        assert np.all(f['H2_mol_m3gas']>=-1e-12)
        assert np.all(f['H2O_mol_m3gas']>=-1e-12)
        assert np.max(h['H2_balance_relative'])<1e-9, h['H2_balance_relative'].max()
        return f,h
    base,hbase=sim('base')
    no_reaction,_=sim('no_reaction','--A','0')
    assert np.max(no_reaction['solid_conversion'])==0
    assert abs(np.max(no_reaction['H2O_mol_m3gas']))<1e-12
    cold,_=sim('no_heat','--heat','0','--A','0')
    assert np.max(np.abs(cold['Tsolid_K']-850.0))<1e-9
    assert np.max(np.abs(cold['Tgas_K']-850.0))<1e-9
    heated,_=sim('more_heat','--heat','900000')
    assert np.max(heated['Tsolid_K'])>np.max(base['Tsolid_K'])
    assert np.mean(heated['solid_conversion'])>np.mean(base['solid_conversion'])
print('PASS: H2 conservation, nonnegativity, conversion bounds, zero reaction,')
print('      no-heating equilibrium, heat sensitivity.')
