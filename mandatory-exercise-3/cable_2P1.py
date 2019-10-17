#!/usr/bin/env python3
# ============================================================================
# File:     cable_2P1.py
# Author:   Erik Johannes Husom
# Created:  2019-10-10
# ----------------------------------------------------------------------------
# Description:
# Third mandatory exercise in IN5270, University of Oslo, fall 2019.
# Details of the problem is found in cable_2P1.pdf
# ============================================================================
import matplotlib.pyplot as plt
import numpy as np
import sympy as sym

def u_e(x):
    return 0.5 * x**2 - x


def phi0(x):
    pass

def phi1(x):
    pass

def phi2(x):
    pass

def phi3(x):
    pass


# x: variable, f: right hand side, constant
x, f = sym.symbols('x f')
# B: boundary function
B = 0
dBdx = sym.diff(B, x)


N = 3
psi = {0: [x**(i+1)*(1-x) for i in range(N+1)]}
psi[1] = [sym.diff(psi_i, x) for psi_i in psi[0]]

def integrand_lhs(psi, i, j):
    return psi[1][i]*psi[1][j]

def integrand_rhs(psi, i):
    return f*psi[0][i] - dBdx*psi[1][i]

Omega = [0, 1]

def Lagrange_polynomial(x, i, points):
    p = 1
    for k in range(len(points)):
        if k != i:
            p *= (x - points[k])/(points[i] - points[k])
    return p

nodes = [0, 0.5, 1]
elements = [[0, 1], [1, 2]]



def plot():
    
    pass

if __name__ == '__main__':

    x = 0.5
    i = 1
    points = [0, 0.5, 1]
    print(Lagrange_polynomial(x, i, points))
