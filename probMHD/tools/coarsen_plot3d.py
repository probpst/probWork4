#!/usr/bin/env python3
"""Sub-sample a single-block PLOT3D grid (Fortran unformatted float64, e.g. the KSTAR
flux-aligned grid kstar_flux_aligned_3M.xyz) by integer factors, keeping the periodic
directions consistent (the factor must divide the periodic point count).

    python3 coarsen_plot3d.py kstar_flux_aligned_3M.xyz kstar_coarse.xyz --fi 4 --fj 5 --fk 5 --periodic 0 1 1
    -> 25 x 60 x 20 nodes
Also: --nphi N keeps only the first N toroidal planes of the sub-sampled grid (a toroidal wedge of
N*fk*3.6 degrees, periodic when the wedge divides 360).
"""
import argparse, numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from readp3d import read_grid

ap = argparse.ArgumentParser()
ap.add_argument('inp'); ap.add_argument('out')
ap.add_argument('--fi', type=int, default=1); ap.add_argument('--fj', type=int, default=1); ap.add_argument('--fk', type=int, default=1)
ap.add_argument('--periodic', type=int, nargs=3, default=[0, 1, 1])
ap.add_argument('--nphi', type=int, default=0, help='keep only the first nphi toroidal planes after sub-sampling (wedge)')
ap.add_argument('--close', action='store_true', help='with --nphi: append the closing plane (nphi+1 planes; Periodic_dup = 1, Periodic_rot = wedge angle)')
ap.add_argument('--wedge-deg', type=float, default=0.0, help='with --nphi 1 --close: the closing plane is plane 1 rotated by this angle (any wedge angle, e.g. 0.5 for an axisymmetric run)')
ap.add_argument('--multiblock', action='store_true', help='input has a leading nblocks record')
a = ap.parse_args()

X, Y, Z = read_grid(a.inp)
n = X.shape
for c, f in enumerate((a.fi, a.fj, a.fk)):
    if a.periodic[c] and n[c] % f != 0:
        raise SystemExit(f'factor {f} does not divide the periodic point count {n[c]} in direction {c+1}')
    if not a.periodic[c] and (n[c] - 1) % f != 0:
        raise SystemExit(f'factor {f} does not divide the cell count {n[c]-1} in direction {c+1}')
sl = (slice(None, None, a.fi), slice(None, None, a.fj), slice(None, None, a.fk))
X, Y, Z = X[sl], Y[sl], Z[sl]
if a.nphi > 0:
    nkeep = a.nphi + (1 if a.close else 0)
    if nkeep > X.shape[2]:
        raise SystemExit('not enough toroidal planes for the wedge plus its closing plane')
    X, Y, Z = X[:, :, :nkeep], Y[:, :, :nkeep], Z[:, :, :nkeep]
if a.wedge_deg > 0:
    if a.nphi != 1 or not a.close:
        raise SystemExit('--wedge-deg needs --nphi 1 --close')
    ang = np.radians(a.wedge_deg)
    X0, Y0, Z0 = X[:, :, :1], Y[:, :, :1], Z[:, :, :1]
    X = np.concatenate([X0, X0*np.cos(ang) - Y0*np.sin(ang)], axis=2)
    Y = np.concatenate([Y0, X0*np.sin(ang) + Y0*np.cos(ang)], axis=2)
    Z = np.concatenate([Z0, Z0], axis=2)
ni, nj, nk = X.shape
with open(a.out, 'wb') as f:
    def rec(b):
        f.write(np.array([len(b)], np.int32).tobytes()); f.write(b); f.write(np.array([len(b)], np.int32).tobytes())
    rec(np.array([ni, nj, nk], np.int32).tobytes())
    rec(np.concatenate([X.ravel(order='F'), Y.ravel(order='F'), Z.ravel(order='F')]).astype(np.float64).tobytes())
R = np.hypot(X, Y)
print(f'{a.out}: {ni} x {nj} x {nk} nodes; R in [{R.min():.4f}, {R.max():.4f}] m, Z in [{Z.min():.4f}, {Z.max():.4f}] m')
if a.nphi > 0:
    phi = np.degrees(np.arctan2(Y[0, 0, :], X[0, 0, :]))
    ext = (phi[1]-phi[0])*(nk - (1 if a.close else 0)) if a.wedge_deg == 0 else a.wedge_deg
    print(f'toroidal planes at phi = {phi.min():.2f} .. {phi.max():.2f} deg; wedge extent {ext:.2f} deg'
          + (f' -> Periodic_rot = 0, 0, {ext:.4f} with Periodic_dup = 0, 0, 1' if a.close else ''))
