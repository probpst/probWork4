#!/usr/bin/env python3
"""Cylindrical m = 2 / n = 1 tearing mode: q(r) = q0 (1 + (r/r0)^4)^(1/2), q0 = 1.2, r0 = 0.6, B_z(0) = 1,
R0 = 3 (L_z = 2 pi R0), p = 0.1 uniform, inviscid.  Cases: (grid, Lundquist number S = 1/eta, scheme, End_T).
'ideal' (eta = 0) runs measure the growth driven by the NUMERICAL resistivity of the scheme on each grid."""
import os
CASES = [  # name, grid, eta, Cnvt_scheme, End_T
    ('ideal_28',   'cyl_28x32x8.xyz', 0.0,     4, 40.0),
    ('ideal_56',   'cyl_56x64x8.xyz', 0.0,     4, 40.0),
    ('S100_28',    'cyl_28x32x8.xyz', 1.0/100, 4, 120.0),
    ('S100_56',    'cyl_56x64x8.xyz', 1.0/100, 4, 120.0),
    ('S300_56',    'cyl_56x64x8.xyz', 1.0/300, 4, 160.0),
    ('S300_56_eno3', 'cyl_56x64x8.xyz', 1.0/300, 5, 160.0),
]
for name, grid, eta, scheme, tend in CASES:
    os.makedirs(name, exist_ok=True)
    open(os.path.join(name, 'probMHD.nml'), 'w').write(f"""$Case_title
  Case_name            = "tearing_{name}"
/end
$Control
  Units                = 1
  Del_T                = 0.0
  End_T                = {tend}
  CFL                  = 0.5
  Gamma                = 1.6666666666666667
/end
$Grid
  Grid_file            = "../{grid}"
  Grid_format          = 2
  Periodic             = 0, 1, 1
  Periodic_dup         = 0, 0, 1
  BC_imin              = 1
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
  Cnvt_scheme          = {scheme}
  Riemann_solver       = 3
  Time_scheme          = 2
  CT_scheme            = 1
  SW_Eq_Subtract       = 1
/end
$Physics
  eta                  = {eta}
  nu                   = 0.0
/end
$Export_Interval
  Output_Intvl         = 0
  Log_Intvl            = 2000
  Monitor_Intvl        = 20
  Output_slice         = 1
/end
""")
print('cases written:', [c[0] for c in CASES])
