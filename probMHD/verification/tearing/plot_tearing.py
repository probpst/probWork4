#!/usr/bin/env python3
"""Growth rates of the m=2/n=1 perturbation from the kinetic-energy history, gamma = d ln(sqrt E_kin)/dt,
fitted over the exponential phase (E_kin between 10x its early minimum and 1/20 of its maximum)."""
import sys, os, glob
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
from readp3d import read_monitor

cases = ['ideal_28', 'ideal_56', 'S100_28', 'S100_56', 'S300_56', 'S300_56_eno3']
fig, ax = plt.subplots(figsize=(8, 5))
lines = []
res = {}
for c in cases:
    f = f'{c}/Output_data/tearing_{c}_monitor.dat'
    if not os.path.exists(f):
        continue
    m = read_monitor(f)
    t, ek = m['tA'], m['ekin']
    ax.semilogy(t, ek, label=c)
    lnE = np.log(np.maximum(ek, 1e-300))
    i0 = np.argmin(ek[: max(len(ek)//3, 2)])
    sel = (t > t[i0]) & (ek > 10*ek[i0]) & (ek < 0.05*ek.max())
    if sel.sum() > 5:
        p = np.polyfit(t[sel], lnE[sel], 1)
        res[c] = 0.5*p[0]
        ax.semilogy(t[sel], np.exp(np.polyval(p, t[sel])), 'k--', lw=0.8)
        lines.append(f'{c:14s} gamma t_A = {0.5*p[0]:.4f}   (fit window t = {t[sel][0]:.1f} .. {t[sel][-1]:.1f})')
    else:
        lines.append(f'{c:14s} no exponential phase found (E_kin max/min = {ek.max()/max(ek.min(),1e-300):.2e})')
ax.set_xlabel('t / t_A'); ax.set_ylabel('E_kin (code units)'); ax.grid(alpha=.3); ax.legend(fontsize=8)
ax.set_title('m = 2 / n = 1 tearing perturbation, periodic cylinder (ENO2 unless noted)')
fig.tight_layout(); fig.savefig('tearing_growth.png', dpi=130)
if 'S100_56' in res and 'S300_56' in res:
    lines.append(f'exponent from S = 100 -> 300 on the 56 x 64 grid: {np.log(res["S300_56"]/res["S100_56"])/np.log(3.0):.3f}  (FKR: -0.600)')
open('tearing_growth.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
