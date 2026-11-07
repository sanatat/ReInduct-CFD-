"""Export 2D quasi-static induction qEM as OpenFOAM volScalarField.

Mapping convention: structured blockMesh (Nx Ny 1), x index fastest; two
empty front/back patches, values on a cell-centered mesh. This is a
one-way prescribed heat field, not an OpenFOAM-native EM solver.
"""
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from v5.em2d import run as induction

def export(dest='v4/openfoam_v2606/channelReactive2D/0/qEM',nx=30,ny=12):
    _,f=induction(nx=61,ny=31,write=False)
    rg=RegularGridInterpolator((f['y'],f['x']),f['q'],bounds_error=False,fill_value=0.)
    xs=(np.arange(nx)+.5)*.12/nx
    ys=(np.arange(ny)+.5)*.02/ny
    xg,yg=np.meshgrid(xs,ys)
    values=rg(np.column_stack((yg.reshape(-1),xg.reshape(-1))))
    path=Path(dest);path.parent.mkdir(parents=True,exist_ok=True)
    lines=['FoamFile','{','    version 2.0;','    format ascii;','    class volScalarField;','    object qEM;','}',
      'dimensions [1 -1 -3 0 0 0 0];','internalField nonuniform List<scalar>',str(len(values)),'(']
    lines += [f'    {v:.12g}' for v in values]
    lines += [');','boundaryField','{',
        '    inlet { type zeroGradient; }','    outlet { type zeroGradient; }',
        '    walls { type zeroGradient; }','    frontAndBack { type empty; }','}']
    path.write_text('\n'.join(lines)+'\n')
    return float(np.mean(values)),float(np.max(values)),len(values)
if __name__=='__main__':print('qEM mean, maximum, cells:',export())
