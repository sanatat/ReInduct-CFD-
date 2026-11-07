#!/usr/bin/env python3
"""Plot the illustrative 1-D reactor benchmark (numpy/matplotlib)."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser()
p.add_argument('directory',type=Path,nargs='?',default=Path('results/base'))
a=p.parse_args()
f=np.genfromtxt(a.directory/'fields.csv',delimiter=',',names=True)
h=np.genfromtxt(a.directory/'history.csv',delimiter=',',names=True)
fig,ax=plt.subplots(figsize=(7,4))
ax.plot(f['z_m'],f['Tsolid_K'],label='Solid')
ax.plot(f['z_m'],f['Tgas_K'],label='Gas')
ax.set(xlabel='Axial coordinate (m)',ylabel='Temperature (K)',title='Reduced-order bed model: final temperature')
ax.legend();fig.tight_layout();fig.savefig(a.directory/'temperature.png',dpi=160);plt.close(fig)
fig,ax=plt.subplots(figsize=(7,4))
ax.plot(f['z_m'],f['solid_conversion'])
ax.set(xlabel='Axial coordinate (m)',ylabel='Solid conversion (-)',ylim=(0,1),title='Representative NiO reduction')
fig.tight_layout();fig.savefig(a.directory/'conversion.png',dpi=160);plt.close(fig)
fig,ax=plt.subplots(figsize=(7,4))
ax.plot(h['time_s'],h['conversion_mean'],label='Mean solid conversion')
ax.set(xlabel='Time (s)',ylabel='Mean conversion (-)',ylim=(0,1))
fig.tight_layout();fig.savefig(a.directory/'conversion_history.png',dpi=160);plt.close(fig)
print('Saved plots to:',a.directory)
