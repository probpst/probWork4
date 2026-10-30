# probMHD verification

Every case is a folder with its `probMHD.nml`, a run script where several runs are involved, the
plot script and the resulting figure / numbers.  Run from the folder with the solver copied in
(`cp ../../probMHD .`).

| folder | test | what is checked | result |
|---|---|---|---|
| `briowu/` | Brio & Wu (1988) shock tube, gamma = 2, 800 cells, t = 0.1 | the wave structure (fast rarefaction, compound wave, contact, slow shock, fast rarefaction) against an 8000-cell run of the same code; ENO2 (probCFD), ENO3 and 1st-order upwind compared; energy conservation | `briowu_result.png`, `briowu_errors.txt` |
| `orszag_tang/` | Orszag-Tang vortex, 256 x 256, periodic, t = 0.5 | total energy conserved to round-off, div B at round-off (constrained transport), the known density / pressure structure | `orszag_tang_result.png`, `orszag_tang_result.txt` |
| `alfven/` | circularly polarised Alfven wave, oblique (tan a = 2), one period, N = 16 ... 128 | order of accuracy of the ENO2 and ENO3 reconstructions on a smooth solution | `alfven_convergence.png`, `alfven_convergence.txt` |
| `tearing/` | m = 2 / n = 1 tearing mode in a periodic cylinder (q0 = 1.2, q = 2 at r = 0.69), S = 300 and 1000 | exponential growth of the kinetic energy; the growth-rate ratio against the FKR scaling S^(-3/5); the axis-hole boundary and the resistive EMF of the constrained transport | `tearing_growth.png`, `tearing_growth.txt` |
| `../cases/KSTAR_equilibrium_2D/` | Solov'ev equilibrium of the KSTAR grid, axisymmetric, 20 Alfven times | the discrete residual flow of the free equilibrium (33 x 60 and 99 x 150) and the exact hold with `SW_Eq_Subtract = 1`; the Solov'ev coefficients against the Python generator (`tools/check_solovev.py`) | `equilibrium_hold.png`, `equilibrium_hold.txt` |

## Results

### Standard MHD tests

| test | result |
|---|---|
| Brio-Wu, 800 cells (`briowu_result.png`) | all five waves at the right positions against the 3200-cell reference; L1(rho) 2.0e-3 (ENO2), 1.3e-3 (ENO3), 7.1e-3 (1st-order upwind); total energy conserved to 5e-15 |
| Orszag-Tang, 256 x 256 (`orszag_tang_result.png`) | dE_tot/E_0 = 0.0 over the run (round-off), max normalised div B 3.5e-12; rho in 0.091 .. 0.495 and p in 0.031 .. 0.509 at t = 0.5 (the usual structure) |
| CP Alfven wave, one period (`alfven_convergence.png`) | ENO3: L1(Bz) 1.16e-2, 1.63e-3, 2.08e-4, 2.64e-5 for N = 16, 32, 64, 128 -> order 2.8 .. 3.0; ENO2 (minmod): 3.6e-2, 1.3e-2, 5.4e-3, 1.7e-3 -> order 1.3 .. 1.7 (the minmod limiter clips the smooth extrema, as it does in probCFD) |

### Sheared equilibria: the numerical resistivity sets the growth

The tearing study (`tearing/`) did NOT give the FKR scaling at the resolution that this environment allows
(one effective CPU core):

* with `eta = 0` the seeded m = 2 / n = 1 perturbation still grows: gamma t_A = 0.29 on 28 x 32 x 8 cells and 0.65
  on 56 x 64 x 8 cells, and the growing structure has m ~ 8 (4 cells per wavelength), not the seeded m = 2;
* the same growth appears on a Cartesian force-free sheared slab (`IC_type = 7`) with the Rusanov flux and is absent
  with the cell-centred (`CT_scheme = 0`) update and with uniform fields.

So on these grids the perturbation dynamics of a current-carrying equilibrium are set by the numerical resistivity of
the constrained-transport scheme (of order c Delta x), which exceeds the physical eta = 1/300 .. 1/1000 by far:
the physical S-scaling cannot be measured until the resistive layer is resolved (several cells inside
delta ~ S^-2/5 a), i.e. a few hundred radial cells for S = 1e3 .. 1e4.  The scripts (`make_cases.py`,
`run_all.sh`, `plot_tearing.py`) are set up for that study: raise the grid with `tools/make_cylinder_grid.py` and run
`S100_*`, `S300_*` on a machine with more cores.

### KSTAR equilibrium (`../cases/KSTAR_equilibrium_2D/equilibrium_hold.png`)

| run | behaviour |
|---|---|
| 33 x 60, free (`SW_Eq_Subtract = 0`) | the discrete residual of J x B - grad p drives a flow that reaches 0.06 v_A within 0.3 t_A and 0.3 v_A by 3 t_A; the first ring of cells at the inner hole loses its field and density (the run was stopped at 3 t_A) |
| 99 x 150, free | the same breakdown, already at 0.8 t_A (the smaller cells at the hole set the pace) |
| 33 x 60, `SW_Eq_Subtract = 1` | exact fixed point: max v / v_A = 1.9e-9 after 8 t_A, E_th and E_mag constant to 1e-6, div B 1e-11 |

The Solov'ev coefficients of the Fortran set-up equal those of the Python grid generator to 1e-14
(`tools/check_solovev.py`), and the printed p_axis = 162.18 kPa and R_axis = 1.8586 m equal the grid metadata.
Tokamak runs must therefore use `SW_Eq_Subtract = 1` together with `SW_Energy = 2`; the perturbed runs in
`../cases/KSTAR_turbulence_switch` show how far a 1e-3 v_A perturbation then develops with and without the
turbulence models.

