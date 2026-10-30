#!/usr/bin/env python3
"""Circularly polarised Alfven wave (Toth 2000; Gardiner & Stone 2005): oblique
propagation at tan(alpha) = 2 on [0, sqrt(5)] x [0, sqrt(5)/2], N x N/2 cells, one period (t = 1).
Writes one case folder per resolution and scheme."""
import os, sys, math
L = math.sqrt(5.0)
for scheme in (4, 5):
    for N in (16, 32, 64, 128):
        d = f'N{N}_s{scheme}'
        os.makedirs(d, exist_ok=True)
        nml = f"""$Case_title
  Case_name            = "Alfven_N{N}_s{scheme}"
/end
$Control
  Units                = 1
  Del_T                = 0.0
  End_T                = 1.0
  CFL                  = 0.4
  Gamma                = 1.6666666666666667
/end
$Grid
  Grid_format          = 3
  Ncell                = {N}, {N//2}, 1
  Box_min              = 0.0, 0.0, 0.0
  Box_max              = {L:.15f}, {L/2:.15f}, {L/N:.15f}
  Periodic             = 1, 1, 1
/end
$Initial_Condition
  IC_type              = 3
  IC_param             = 2.0, 0.1
/end
$Numerical_Scheme
  Cnvt_scheme          = {scheme}
  Riemann_solver       = 3
  Time_scheme          = 3
  CT_scheme            = 1
/end
$Export_Interval
  Output_Intvl         = 0
  Log_Intvl            = 1000
  Monitor_Intvl        = 10
  Output_slice         = 1
/end
"""
        open(os.path.join(d, 'probMHD.nml'), 'w').write(nml)
print('cases written')
