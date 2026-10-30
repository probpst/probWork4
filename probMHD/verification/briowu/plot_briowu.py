#!/usr/bin/env python3
"""Brio-Wu shock tube: probMHD (800 cells, ENO2 / ENO3 / 1st-order upwind) against the
3200-cell reference run of the same code at t = 0.1."""
import sys, os, glob
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
from readp3d import read_slice

def last_slice(d):
    f = sorted(glob.glob(os.path.join(d, 'Output_data', '*_slice_*.dat')))[-1]
    s = read_slice(f)
    n = int(s['i'].max())
    x = (np.arange(n) + 0.5)/n
    return x, s

runs = [('ref', '3200 cells, ENO2 (reference)', 'k-'), ('.', '800 cells, ENO2 (Cnvt_scheme 4)', 'r.'),
        ('eno3', '800 cells, ENO3 (Cnvt_scheme 5)', 'b.'), ('upw1', '800 cells, 1st-order upwind', 'g.')]
fig, ax = plt.subplots(2, 3, figsize=(14, 8))
ax = ax.ravel()
vars_ = [('rho', r'$\rho$'), ('p', 'p'), ('u', r'$v_x$'), ('v', r'$v_y$'), ('By', r'$B_y$'), ('T', r'$p/\rho$')]
errs = {}
for d, lab, st in runs:
    try:
        x, s = last_slice(d)
    except IndexError:
        continue
    for a, (v, name) in zip(ax, vars_):
        a.plot(x, s[v], st, ms=2, lw=1, label=lab)
        a.set_title(name)
    if d == 'ref':
        xr, sr = x, s
    else:
        # L1 error of rho against the reference (averaged onto the coarse cells)
        rr = sr['rho'].reshape(len(x), -1).mean(axis=1)
        errs[lab] = np.abs(s['rho'] - rr).mean()
for a in ax:
    a.set_xlim(0, 1); a.grid(alpha=.3)
ax[0].legend(fontsize=8)
fig.suptitle('Brio & Wu shock tube, gamma = 2, t = 0.1  (probMHD, HLLD + CT, SSP-RK3)')
fig.tight_layout()
fig.savefig('briowu_result.png', dpi=130)
with open('briowu_errors.txt', 'w') as f:
    for k, v in errs.items():
        f.write(f'L1(rho) vs 3200-cell reference  {k:40s} {v:.4e}\n')
print(open('briowu_errors.txt').read())
