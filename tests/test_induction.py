import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'em'))
from induction1d import solve
import numpy as np

def test_exact_complex_slab():
    r=solve(n=801)
    x=r['x']; L=x[-1]; d=r['delta']; k=(1+1j)/d
    exact=r['params']['H0_A_m']*np.cosh(k*(L-x))/np.cosh(k*L)
    err=np.max(np.abs(r['H']-exact))/r['params']['H0_A_m']
    assert err<0.002,(err,r['delta'])
    assert np.all(r['q']>=0)
    assert r['power_per_area']>0
    return err

def test_refinement():
    def error(n):
        r=solve(n=n)
        k=(1+1j)/r['delta'];L=r['x'][-1]
        exact=r['params']['H0_A_m']*np.cosh(k*(L-r['x']))/np.cosh(k*L)
        return np.max(np.abs(r['H']-exact))/r['params']['H0_A_m']
    assert error(401)<error(101)
    assert error(801)<error(401)
    return [error(n) for n in (101,401,801)]

if __name__=='__main__':
    print('EM exact solution max relative error',test_exact_complex_slab())
    print('EM grid refinement errors',test_refinement())
