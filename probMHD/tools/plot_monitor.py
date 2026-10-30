#!/usr/bin/env python3
"""Plot the global history of a run: energies, max |v|, div B, perturbation energies."""
import sys, os, glob
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from readp3d import read_monitor

f = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob('Output_data/*_monitor.dat'))[0]
m = read_monitor(f)
fig, ax = plt.subplots(2, 2, figsize=(11, 7))
ax[0, 0].plot(m['tA'], m['ekin'], label='E_kin'); ax[0, 0].plot(m['tA'], m['emagpert'], label='E_mag,pert')
ax[0, 0].set_yscale('log'); ax[0, 0].set_title('perturbation energies (code units)'); ax[0, 0].legend()
ax[0, 1].plot(m['tA'], (m['etot'] - m['etot'][0])/abs(m['etot'][0])); ax[0, 1].set_title('(E_tot - E_tot0)/E_tot0')
ax[1, 0].plot(m['tA'], m['vmax']); ax[1, 0].set_yscale('log'); ax[1, 0].set_title('max |v| / v_A')
ax[1, 1].plot(m['tA'], m['divb']); ax[1, 1].set_yscale('log'); ax[1, 1].set_title('max |div B| h / |B|')
for a in ax.ravel():
    a.set_xlabel('t / t_A'); a.grid(alpha=.3)
fig.tight_layout(); out = os.path.splitext(os.path.basename(f))[0] + '.png'
fig.savefig(out, dpi=130); print(out)
