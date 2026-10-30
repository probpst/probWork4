#!/usr/bin/env python3
"""L1 error of B_z after one period against the initial condition; convergence order."""
import sys, os, glob
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
from readp3d import read_slice

res = {}
for scheme in (4, 5):
    Ns, errs = [], []
    for N in (16, 32, 64, 128):
        d = f'N{N}_s{scheme}'
        files = sorted(glob.glob(os.path.join(d, 'Output_data', '*_slice_*.dat')))
        if len(files) < 2:
            continue
        s0 = read_slice(files[0]); s1 = read_slice(files[-1])
        e = np.mean(np.abs(s1['Bz'] - s0['Bz']))
        Ns.append(N); errs.append(e)
    res[scheme] = (np.array(Ns), np.array(errs))
fig, ax = plt.subplots(figsize=(6, 4.5))
lines = []
for scheme, (Ns, errs) in res.items():
    if len(Ns) == 0:
        continue
    lab = {4: 'ENO2 (Cnvt_scheme 4)', 5: 'ENO3 (Cnvt_scheme 5)'}[scheme]
    ax.loglog(Ns, errs, 'o-', label=lab)
    orders = -np.diff(np.log(errs))/np.diff(np.log(Ns))
    lines.append(f'{lab}: N = {Ns.tolist()}  L1(Bz) = {[f"{e:.3e}" for e in errs]}  orders = {[f"{o:.2f}" for o in orders]}')
if len(res[4][0]):
    N0, e0 = res[4][0][0], res[4][1][0]
    ax.loglog(res[4][0], e0*(res[4][0]/N0)**-2.0, 'k--', lw=0.8, label='2nd order')
ax.set_xlabel('N (cells along x)'); ax.set_ylabel('L1 error of $B_z$ after one period')
ax.set_title('Circularly polarised Alfven wave, oblique (tan a = 2)')
ax.grid(alpha=.3, which='both'); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig('alfven_convergence.png', dpi=130)
open('alfven_convergence.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
