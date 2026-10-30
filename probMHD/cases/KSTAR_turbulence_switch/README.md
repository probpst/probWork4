# KSTAR_turbulence_switch — the turbulence-model switch on the axisymmetric KSTAR equilibrium

Three runs of the same perturbed equilibrium (33 x 60 flux-aligned cells, 0.5 deg wedge,
m = 3 / n = 0 velocity perturbation of 1e-3 v_A, resistive S = 1e5, well-balanced, internal-energy
equation, 15 Alfven times) that differ only in `$Turbulence`:

| folder | `SW_Turb_MHD` | model |
|---|---|---|
| `turb0` | 0 | none: the ENO scheme alone (implicit LES) |
| `turb1` | 1 | MHD Smagorinsky: `nu_sgs = (Cs Delta)^2 |S|`, `eta_sgs = (C_eta Delta)^2 |J|`, `chi_sgs = nu_sgs/Pr_sgs` |
| `turb3` | 3 | anomalous perpendicular transport: `chi_perp += chi_anom`, `D_perp += D_anom`, `nu += nu_anom` (constants) |

Run each folder with `mpirun -np 1 ./probMHD`, then `python3 ../../tools/plot_monitor.py Output_data/<case>_monitor.dat`
or `python3 compare.py` for the three kinetic-energy histories in one figure.

## Result of the delivered runs (15 Alfven times, one effective CPU core)

| run | E_kin at the end | max v / v_A | E_th / E_th(0) | min rho / rho(0) |
|---|---|---|---|---|
| turb0 (none) | 1.6e-5 | 0.36 | 2.5 | 0.075 |
| turb1 (Smagorinsky, Cs = C_eta = 0.1) | 1.6e-5 | 0.37 | 2.4 | 0.10 |
| turb3 (anomalous chi = 1, D = 0.3, nu = 0.3 m2/s) | 8.6e-6 | 0.33 | 1.8 | 0.08 |

The kinetic energy stays small (the 1e-3 v_A perturbation grows by a factor ~ 4 in amplitude over the run), but a
fast flow (~ 0.3 v_A) appears within 1.5 t_A in the first ring of cells at the inner hole (rho_pol = 0.05) and the
density of that ring drains while the thermal energy rises - the same inner-boundary breakdown as the free
equilibrium runs, only slower.  The anomalous-transport model (turb3) damps it most; the sub-grid model (turb1)
hardly changes it because its coefficients vanish with the velocity gradients.  The figure `turbulence_switch.png`
shows the three histories.  Treat the inner-hole boundary as the current weak point of the tokamak set-up
(see the main README).
