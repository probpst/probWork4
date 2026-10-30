#!/usr/bin/env python3
"""Periodic cylinder grid for probMHD (tearing-mode tests): r in [r_in, r_out],
theta periodic (no closing node), z periodic over L_z (closing plane repeated: Periodic_dup = 0, 0, 1).
Written as PLOT3D Fortran unformatted float64, single block (Grid_format = 2).

    python3 make_cylinder_grid.py out.xyz --nr 48 --ntheta 64 --nz 16 --rin 0.05 --rout 1.0 --Lz 18.85
"""
import argparse, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('out')
ap.add_argument('--nr', type=int, default=48, help='radial cells')
ap.add_argument('--ntheta', type=int, default=64, help='poloidal cells (periodic, no duplicate node)')
ap.add_argument('--nz', type=int, default=16, help='axial cells (periodic, no duplicate node)')
ap.add_argument('--rin', type=float, default=0.05)
ap.add_argument('--rout', type=float, default=1.0)
ap.add_argument('--Lz', type=float, default=2*np.pi*3.0)
ap.add_argument('--pack', type=float, default=0.0, help='radial clustering towards r_s (0 = uniform)')
ap.add_argument('--rs', type=float, default=0.5, help='clustering centre')
a = ap.parse_args()

u = np.linspace(0, 1, a.nr + 1)
if a.pack > 0:
    # tanh clustering around rs
    s = (a.rs - a.rin)/(a.rout - a.rin)
    r = a.rin + (a.rout - a.rin)*(u + a.pack*np.sin(np.pi*u)*(s - u))
else:
    r = a.rin + (a.rout - a.rin)*u
th = 2*np.pi*np.arange(a.ntheta)/a.ntheta           # periodic, no duplicate
z = a.Lz*np.arange(a.nz + 1)/a.nz                    # periodic WITH the closing plane (a translation: Periodic_dup = 0, 0, 1)
R, TH, Z = np.meshgrid(r, th, z, indexing='ij')
X = R*np.cos(TH); Y = R*np.sin(TH)
# right-handed (r, theta, z) -> (i, j, k) with J > 0
ni, nj, nk = X.shape
with open(a.out, 'wb') as f:
    def rec(b):
        f.write(np.array([len(b)], np.int32).tobytes()); f.write(b); f.write(np.array([len(b)], np.int32).tobytes())
    rec(np.array([ni, nj, nk], np.int32).tobytes())
    rec(np.concatenate([X.ravel(order='F'), Y.ravel(order='F'), Z.ravel(order='F')]).astype(np.float64).tobytes())
print(f'{a.out}: {ni} x {nj} x {nk} nodes, cells {ni-1} x {nj} x {nk-1}, r in [{a.rin}, {a.rout}], Lz = {a.Lz:.4f}')
print('name list: Periodic = 0, 1, 1   Periodic_dup = 0, 0, 1   (theta closed without a duplicate node, z closed by the repeated plane)')
