import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'em'))
from induction1d import solve
from coupled_thermal import thermal
import numpy as np

def test_coupling():
    em=solve(n=101)
    T0=thermal(em,duration=0)
    T=thermal(em,duration=10)
    assert np.allclose(T0,850.)
    assert T.max()>850.
    em0=dict(em);em0['q']=np.zeros_like(em['q'])
    assert np.max(abs(thermal(em0,duration=10)-850))<1e-8
    # Induction heating scales with H0 squared for constant material parameters.
    stronger=solve(n=101,h0=20_000.)
    assert np.allclose(stronger['q'],4*em['q'])
    assert thermal(stronger,duration=10).max()>T.max()
    return float(T.max())
if __name__=='__main__':print('EM->heat passed; Tmax K:',test_coupling())
