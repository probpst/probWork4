#!/usr/bin/env python3
"""Orszag-Tang vortex at t = 0.5: density / pressure / |B| maps from the mid-plane slice,
energy conservation and div B from the monitor."""
import sys, os, glob
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
from readp3d import read_slice, read_monitor

f = sorted(glob.glob('Output_data/OrszagTang_slice_*.dat'))[-1]
s = read_slice(f)
ni, nj = int(s['i'].max()), int(s['j'].max())
def fld(k):
    a = np.zeros((ni, nj)); a[s['i'].astype(int)-1, s['j'].astype(int)-1] = s[k]; return a.T
rho, p = fld('rho'), fld('p')
Bm = np.sqrt(fld('Bx')**2 + fld('By')**2)
m = read_monitor('Output_data/OrszagTang_monitor.dat')
fig, ax = plt.subplots(2, 2, figsize=(10, 9))
for a, (d, name) in zip(ax.ravel()[:3], [(rho, r'$\rho$'), (p, 'p'), (Bm, '|B|')]):
    im = a.imshow(d, origin='lower', extent=[0, 1, 0, 1], cmap='viridis'); a.set_title(name + ' at t = 0.5'); plt.colorbar(im, ax=a, shrink=0.8)
a = ax[1, 1]
a.plot(m['t'], (m['etot'] - m['etot'][0])/m['etot'][0], label='(E_tot - E_0)/E_0')
a.plot(m['t'], m['divb'], label='max |div B| h/|B|')
a.set_yscale('symlog', linthresh=1e-15); a.legend(); a.grid(alpha=.3); a.set_xlabel('t'); a.set_title('conservation')
fig.suptitle(f'Orszag-Tang vortex, {ni} x {nj}, ENO2 + HLLD + CT (GS) + SSP-RK3')
fig.tight_layout(); fig.savefig('orszag_tang_result.png', dpi=130)
with open('orszag_tang_result.txt', 'w') as g:
    g.write(f'grid {ni} x {nj}, t_end = {m["t"][-1]:.4f}\n')
    g.write(f'max |dE_tot/E_0| over the run = {np.abs((m["etot"] - m["etot"][0])/m["etot"][0]).max():.3e}\n')
    g.write(f'max |div B| (normalised) over the run = {m["divb"].max():.3e}\n')
    g.write(f'rho range at t_end: {rho.min():.4f} .. {rho.max():.4f};  p range: {p.min():.4f} .. {p.max():.4f}\n')
print(open('orszag_tang_result.txt').read())
