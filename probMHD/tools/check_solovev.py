#!/usr/bin/env python3
"""Cross-check of the Fortran Solov'ev coefficients (printed in the probMHD log) against the
Python generator of the KSTAR grid (generate_kstar_grid.py, Cerfon & Freidberg construction).

    python3 check_solovev.py log.txt
"""
import sys, re, numpy as np
import sympy as sp

def python_coeffs(eps, kap, dlt, A):
    x, y = sp.symbols('x y', positive=True)
    psi_p = x**4/8 + A*(sp.Rational(1, 2)*x**2*sp.log(x) - x**4/8)
    basis = [sp.Integer(1), x**2, y**2 - x**2*sp.log(x), x**4 - 4*x**2*y**2,
             2*y**4 - 9*y**2*x**2 + 3*x**4*sp.log(x) - 12*x**2*y**2*sp.log(x),
             x**6 - 12*x**4*y**2 + 8*x**2*y**4,
             8*y**6 - 140*y**4*x**2 + 75*y**2*x**4 - 15*x**6*sp.log(x) + 180*x**4*y**2*sp.log(x) - 120*x**2*y**4*sp.log(x)]
    alpha = np.arcsin(dlt)
    N1 = -(1 + alpha)**2/(eps*kap**2); N2 = (1 - alpha)**2/(eps*kap**2); N3 = -kap/(eps*np.cos(alpha)**2)
    P_out, P_in, P_top = (1+eps, 0.0), (1-eps, 0.0), (1-dlt*eps, kap*eps)
    def constraints(f):
        fl = sp.lambdify((x, y), f, 'numpy'); fx = sp.lambdify((x, y), sp.diff(f, x), 'numpy')
        fy = sp.lambdify((x, y), sp.diff(f, y), 'numpy'); fxx = sp.lambdify((x, y), sp.diff(f, x, 2), 'numpy')
        fyy = sp.lambdify((x, y), sp.diff(f, y, 2), 'numpy')
        return [fl(*P_out), fl(*P_in), fl(*P_top), fx(*P_top), fyy(*P_out) + N1*fx(*P_out),
                fyy(*P_in) + N2*fx(*P_in), fxx(*P_top) + N3*fy(*P_top)]
    M = np.array([constraints(b) for b in basis]).T
    return np.linalg.solve(M, -np.array(constraints(psi_p)))

log = open(sys.argv[1]).read()
m = re.search(r'coefficients\s+(.*)', log)
fort = np.array([float(v) for v in m.group(1).split()])
py = python_coeffs(0.5/1.8, 1.8, 0.6, -0.155)
print('Fortran :', fort)
print('Python  :', py)
print('max |diff| / max |c| = %.3e' % (np.abs(fort - py).max()/np.abs(py).max()))
