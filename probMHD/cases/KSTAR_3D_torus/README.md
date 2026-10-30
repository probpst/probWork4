# KSTAR_3D_torus — full torus, 3-D, m = 3 / n = 1 perturbation

Grid `kstar_33x60x20_torus.xyz`: the delivered KSTAR grid sub-sampled 3 x 5 x 5 (33 x 60 x 20 cells,
20 toroidal planes of 18 deg, periodic in j and k without duplicate nodes exactly as the original
file).  Solov'ev equilibrium, internal-energy equation, equilibrium subtraction, resistivity
S = 1e6, velocity perturbation 1e-3 v_A with cos(3 theta - phi) on the q = 3 surface region
(the equilibrium has q from 2.63 to 5.36, so the 3/1, 4/1, 5/1 rational surfaces exist; there is
no q = 1 or q = 2 surface, see the grid's README).

    cp ../../probMHD . && mpirun -np 4 ./probMHD
    python3 ../../tools/plot_monitor.py Output_data/KSTAR_3D_monitor.dat

The 18 deg toroidal spacing of the coarse grid gives a chord error of order (Delta phi)^2 / 8 ~ 1e-2 in
the hoop force and a spurious B_n / |B| ~ 3e-3 on the flux-surface boundaries; the full-resolution
grid (3.6 deg) has 1.2e-4 / 1e-4.  The equilibrium subtraction removes the residual of the initial
state in either case.  With the full 3M-point grid the Alfven CFL gives dt ~ 1e-10 s (see the
main README on the explicit time step).
