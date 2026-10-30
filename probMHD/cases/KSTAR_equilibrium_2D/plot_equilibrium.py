#!/usr/bin/env python3
"""KSTAR axisymmetric equilibrium hold: max |v|/v_A and the thermal / magnetic energies against time
for the coarse grid without and with equilibrium subtraction and for the fine grid; the pressure
and |B| maps of the coarse run on the (R, Z) plane."""
import sys, os, glob
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
from readp3d import read_monitor, read_grid, read_function, read_names

runs = [('nosub', 'KSTAR_eq2D', '33 x 60, free'), ('sub', 'KSTAR_eq2D_sub', '33 x 60, SW_Eq_Subtract = 1'),
        ('fine', 'KSTAR_eq2D_fine', '99 x 150, free')]
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
lines = []
for d, case, lab in runs:
    f = f'{d}/Output_data/{case}_monitor.dat'
    if not os.path.exists(f):
        continue
    m = read_monitor(f)
    ax[0].semilogy(m['tA'], np.maximum(m['vmax'], 1e-20), label=lab)
    ax[1].plot(m['tA'], m['eth']/m['eth'][0], label=lab)
    lines.append(f'{lab:32s} t_end = {m["tA"][-1]:6.2f} t_A  max|v|/v_A: end {m["vmax"][-1]:.3e}  peak {m["vmax"].max():.3e}  '
                 f'E_th/E_th0 end {m["eth"][-1]/m["eth"][0]:.4f}  E_mag/E_mag0 end {m["emag"][-1]/m["emag"][0]:.6f}  max divB {m["divb"].max():.2e}')
ax[0].set_xlabel('t / t_A'); ax[0].set_ylabel('max |v| / v_A'); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[0].set_title('equilibrium hold: residual flow')
ax[1].set_xlabel('t / t_A'); ax[1].set_ylabel('E_th / E_th(0)'); ax[1].grid(alpha=.3); ax[1].set_title('thermal energy')
# pressure map of the coarse run (last volume file)
try:
    g = sorted(glob.glob('nosub/Volume_solutions_grids/*_cells.grd'))[0]
    X, Y, Z = read_grid(g)
    p3 = sorted(glob.glob('nosub/Volume_solutions_grids/KSTAR_eq2D_*.p3d'))[-1]
    v = read_function(p3); names = read_names(sorted(glob.glob('nosub/Volume_solutions_grids/*.nam'))[0])
    R = np.hypot(X, Y)[:, :, 0]; Zc = Z[:, :, 0]
    ip = names.index('p [Pa]')
    tc = ax[2].tricontourf(R.ravel(), Zc.ravel(), v[ip][:, :, 0].ravel()/1e3, 30, cmap='magma')
    plt.colorbar(tc, ax=ax[2], label='p [kPa]')
    ax[2].plot(R, Zc, 'w-', lw=0.2, alpha=0.5); ax[2].plot(R.T, Zc.T, 'w-', lw=0.2, alpha=0.5)
    ax[2].set_aspect('equal'); ax[2].set_xlabel('R [m]'); ax[2].set_ylabel('Z [m]'); ax[2].set_title('p at t_end (33 x 60 cells)')
except (IndexError, ValueError) as e:
    ax[2].set_title(f'(no volume file: {e})')
fig.tight_layout(); fig.savefig('equilibrium_hold.png', dpi=130)
open('equilibrium_hold.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
