#!/usr/bin/env python3
"""Cylindrical m = 2 / n = 1 tearing mode: q(r) = q0 (1 + (r/r0)^4)^(1/2), q0 = 1.2, r0 = 0.6,
B_z(0) = 1, R0 = 3 (L_z = 2 pi R0), p = 0.1 (uniform), resistive, inviscid.
Lundquist numbers S = 1/eta (L = 1, v_A(0) = 1).  Grid 28 x 32 x 8 cells, r in [0.15, 1]."""
import os
for S in (300, 1000):
    d = f'S{S}'
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, 'probMHD.nml'), 'w').write(f"""$Case_title
  Case_name            = "tearing_S{S}"
/end
$Control
  Units                = 1
  Del_T                = 0.0
  End_T                = {300.0 if S == 1000 else 200.0}
  CFL                  = 0.5
  Gamma                = 1.6666666666666667
/end
$Grid
  Grid_file            = "../cyl_28x32x8.xyz"
  Grid_format          = 2
  Periodic             = 0, 1, 1
  Periodic_dup         = 0, 0, 1
  BC_imin              = 4
  BC_imax              = 1
/end
$Initial_Condition
  IC_type              = 4
  IC_param             = 1.2, 0.6, 2.0, 18.849555921538759, 0.1
  p_init               = 0.1
  Perturb_amp          = 1.0e-4
  Perturb_m            = 2
  Perturb_n            = 1
/end
$Numerical_Scheme
  Cnvt_scheme          = 4
  Riemann_solver       = 3
  Time_scheme          = 2
  CT_scheme            = 1
  SW_Eq_Subtract       = 1
/end
$Physics
  eta                  = {1.0/S}
  nu                   = 0.0
/end
$Export_Interval
  Output_Intvl         = 0
  Log_Intvl            = 2000
  Monitor_Intvl        = 20
  Output_slice         = 1
/end
""")
print('cases written')
