#!/usr/bin/env python3
# ============================================================================
# File:     Skeleton for Python files
# Author:   Erik Johannes Husom
# Created:  2019-06-19
# ----------------------------------------------------------------------------
# Description:
#
# ============================================================================
import numpy as np
import sympy as sym
sym.init_printing(use_unicode=True)


x0 = 0
x1 = 0.5
x2 = 1

def phi0(x):
    return np.piecewise(x, [x < 0, 
                            np.logical_and(x >= 0, x1 >= x), 
                            x > x1],
                        [0,
                         lambda x: -2*x + 1,
                         0])


def phi1(x):
    return np.piecewise(x, [x < 0, 
                            np.logical_and(x >= 0, x1 >= x),
                            np.logical_and(x > x1, x < x2), 
                            x > x2],
                        [0,
                         lambda x: 2*x,
                         lambda x: -2*x + 2,
                         0])

def phi2(x):
    return np.piecewise(x, [x < x1, 
                            np.logical_and(x >= x1, x2 >= x), 
                            x > x2],
                        [0,
                         lambda x: 2*x - 1,
                         0])

x = sym.Symbol('x')
n = 3
# phi = [(1 - 2*x)*(1 - x), 2*x*(2 - 2*x), x*(2*x - 1)]
phi = [phi0(x), phi1(x), phi2(x)]

f = 1 + 2*x-x**2


A = sym.zeros(n, n)
b = sym.zeros(n)

for i in range(n):
    for j in range(n):
        A[i, j] = sym.integrate(phi[i]*phi[j], (x, 0, 1))
    b[i] = sym.integrate(phi[i]*f, (x, 0, 1))
print('A:')
A


