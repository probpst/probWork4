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

RESULTS_PLACEHOLDER
