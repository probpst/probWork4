#!/usr/bin/env python3
"""Kinetic-energy histories of the three turbulence settings."""
import sys, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
from readp3d import read_monitor
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for t, lab in [(0, 'SW_Turb_MHD = 0 (none)'), (1, '1 (MHD Smagorinsky)'), (3, '3 (anomalous transport)')]:
    f = f'turb{t}/Output_data/KSTAR_turb{t}_monitor.dat'
    if not os.path.exists(f):
        continue
    m = read_monitor(f)
    ax[0].semilogy(m['tA'], m['ekin'], label=lab)
    ax[1].plot(m['tA'], m['eth']/m['eth'][0], label=lab)
ax[0].set_xlabel('t / t_A'); ax[0].set_ylabel('E_kin (code units)'); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[1].set_xlabel('t / t_A'); ax[1].set_ylabel('E_th / E_th(0)'); ax[1].grid(alpha=.3)
fig.suptitle('KSTAR axisymmetric equilibrium, m = 3 perturbation: turbulence switch')
fig.tight_layout(); fig.savefig('turbulence_switch.png', dpi=130); print('turbulence_switch.png')
