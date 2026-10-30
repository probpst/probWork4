# KSTAR_equilibrium_2D — axisymmetric hold of the Solov'ev equilibrium on the KSTAR grid

The delivered flux-aligned grid (`kstar_flux_aligned_3M.xyz`, 100 x 300 x 100 points) sub-sampled to
33 x 60 cells in (i, j) and reduced to ONE toroidal cell: a 0.5 deg wedge closed by rotational
periodicity (`tools/coarsen_plot3d.py ... --nphi 1 --close --wedge-deg 0.5`, `Periodic_rot = 0, 0, 0.5`).
The initial state is the Cerfon-Freidberg Solov'ev equilibrium of the grid generator (R0 1.8 m,
a 0.5 m, kappa 1.8, delta 0.6, B0 3.5 T, A = -0.155, psi0 = 8.3188 Wb/rad: p_axis 162 kPa,
R_axis 1.8586 m - the code prints both and `tools/check_solovev.py` compares the seven coefficients
with the Python generator), n_axis = 5e19 m-3.  Code units: L0 = 1 m, B0 = 3.5 T,
rho0 = 1.672e-7 kg/m3 -> v_A0 = 7.6e6 m/s, t_A = 1.3e-7 s.

| folder | grid | `SW_Eq_Subtract` | purpose |
|---|---|---|---|
| `nosub` | 33 x 60 | 0 | the free discrete equilibrium: the residual of J x B - grad p of the coarse grid drives a flow |
| `sub`   | 33 x 60 | 1 | the same with the initial residual subtracted: the equilibrium is an exact fixed point |
| `fine`  | 99 x 150 | 0 | the free equilibrium on a 3 x finer grid: the residual flow drops with resolution |

`SW_Energy = 2` (internal-energy equation) in all three: with beta ~ 1e-3 the total-energy form loses the
pressure (see the main README).  Run each folder with `mpirun -np 1 ./probMHD` (the fine grid with
`-np 2`), then `python3 plot_equilibrium.py` (figure `equilibrium_hold.png`, numbers in
`equilibrium_hold.txt`).
