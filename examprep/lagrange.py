#!/usr/bin/env python3
# ============================================================================
# File:     Skeleton for Python files
# Author:   Erik Johannes Husom
# Created:  2019-06-19
# ----------------------------------------------------------------------------
# Description:
#
# ============================================================================
import sympy as sym

def Lagrange_polynomial(x, i, points):
    p = 1
    for k in range(len(points)):
        if k != i:
            p *= (x - points[k])/(points[i] - points[k])
    return p



def Lagrange_polynomials_01(x, N):
    if isinstance(x, sym.Symbol):
        h = sym.Rational(1, N-1)
    else:
        h = 1.0/(N-1)
    points = [i*h for i in range(N)]
    psi = [Lagrange_polynomial(x, i, points) for i in range(N)]
    return psi, points


x = sym.Symbol('x')
psi, points = Lagrange_polynomials_01(x, N=3)
print(psi)
print(points)
