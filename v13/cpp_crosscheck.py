"""Check compiled C++17 Shomate output against the independent Python code."""
import subprocess
import numpy as np
from v13.gas_properties import cp_molar,enthalpy_molar

def run(binary='v13/results/gas_shomate_check'):
    text=subprocess.check_output([binary],text=True)
    x=np.array([[float(v) for v in row.split(',')] for row in text.strip().splitlines()])
    np.testing.assert_array_equal(x[:,0],np.arange(5))
    np.testing.assert_allclose(x[:,1],cp_molar(1080.),rtol=1e-12,atol=1e-10)
    np.testing.assert_allclose(x[:,2],enthalpy_molar(1080.),rtol=1e-12,atol=1e-8)
    print('C++17/Python Shomate cross-check: PASS')
if __name__=='__main__':run()
