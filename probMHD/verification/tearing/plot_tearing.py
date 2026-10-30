#!/usr/bin/env python3
"""Growth rate of the m=2/n=1 tearing mode from the kinetic-energy history (gamma = d ln sqrt(E_kin) / dt)
for two Lundquist numbers; the FKR scaling gives gamma ~ S^(-3/5)."""
import sys, os, glob
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
from readp3d import read_monitor

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
res = {}
for S in (300, 1000):
    f = f'S{S}/Output_data/tearing_S{S}_monitor.dat'
    if not os.path.exists(f):
        continue
    m = read_monitor(f)
    t, ek = m['tA'], m['ekin']
    ax[0].semilogy(t, ek, label=f'S = {S}')
    # linear phase: fit ln(E_kin) where E_kin is between 10x the minimum (after the initial transient) and 1/10 of max
    lnE = np.log(np.maximum(ek, 1e-300))
    i0 = np.argmin(ek[: len(ek)//2])
    sel = (t > t[i0]) & (ek > 10*ek[i0]) & (ek < 0.1*ek.max())
    if sel.sum() > 5:
        p = np.polyfit(t[sel], lnE[sel], 1)
        gam = 0.5*p[0]
        res[S] = gam
        ax[0].semilogy(t[sel], np.exp(np.polyval(p, t[sel])), 'k--', lw=0.8)
ax[0].set_xlabel('t / t_A'); ax[0].set_ylabel('E_kin (code units)'); ax[0].set_title('m=2/n=1 tearing mode, periodic cylinder')
ax[0].legend(); ax[0].grid(alpha=.3)
lines = []
if len(res) == 2:
    Ss = np.array(sorted(res)); gs = np.array([res[s] for s in Ss])
    slope = np.log(gs[1]/gs[0])/np.log(Ss[1]/Ss[0])
    ax[1].loglog(Ss, gs, 'o-', label='probMHD')
    ax[1].loglog(Ss, gs[0]*(Ss/Ss[0])**(-0.6), 'k--', label=r'$S^{-3/5}$ (FKR)')
    ax[1].set_xlabel('S'); ax[1].set_ylabel(r'$\gamma\, t_A$'); ax[1].legend(); ax[1].grid(alpha=.3, which='both')
    lines.append(f'gamma(S=300) = {gs[0]:.4e}   gamma(S=1000) = {gs[1]:.4e}   measured exponent = {slope:.3f}  (FKR: -0.600)')
for S, g in res.items():
    lines.append(f'S = {S}: gamma t_A = {g:.4e}')
fig.tight_layout(); fig.savefig('tearing_growth.png', dpi=130)
open('tearing_growth.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
