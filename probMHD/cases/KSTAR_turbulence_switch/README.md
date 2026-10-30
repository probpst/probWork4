# KSTAR_turbulence_switch — the turbulence-model switch on the axisymmetric KSTAR equilibrium

Three runs of the same perturbed equilibrium (33 x 60 flux-aligned cells, 0.5 deg wedge,
m = 3 / n = 0 velocity perturbation of 1e-3 v_A, resistive S = 1e5, well-balanced, internal-energy
equation) that differ only in `$Turbulence`:

| folder | `SW_Turb_MHD` | model |
|---|---|---|
| `turb0` | 0 | none: the ENO scheme alone (implicit LES) |
| `turb1` | 1 | MHD Smagorinsky: `nu_sgs = (Cs Delta)^2 |S|`, `eta_sgs = (C_eta Delta)^2 |J|`, `chi_sgs = nu_sgs/Pr_sgs` |
| `turb3` | 3 | anomalous perpendicular transport: `chi_perp += chi_anom`, `D_perp += D_anom`, `nu += nu_anom` (constants) |

Run each folder with `mpirun -np 1 ./probMHD`, then `python3 ../../tools/plot_monitor.py Output_data/<case>_monitor.dat`
or `python3 compare.py` for the three kinetic-energy histories in one figure.
