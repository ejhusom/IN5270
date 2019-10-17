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



def phi_r(r, X, d):
    if isinstance(X, sym.Symbol):
        h = sym.Rational(1, d)  # node spacing
        nodes = [2*i*h - 1 for i in range(d+1)]
    else:
        # assume X is numeric: use floats for nodes
        nodes = np.linspace(-1, 1, d+1)
    return Lagrange_polynomial(X, r, nodes)

def Lagrange_polynomial(x, i, points):
    p = 1
    for k in range(len(points)):
        if k != i:
            p *= (x - points[k])/(points[i] - points[k])
    return p

def basis(d=1):
    """Return the complete basis."""
    X = sym.Symbol('X')
    phi = [phi_r(r, X, d) for r in range(d+1)]
    return phi

def element_matrix(phi, Omega_e, symbolic=True):
    n = len(phi)
    A_e = sym.zeros(n, n)
    X = sym.Symbol('X')
    if symbolic:
        h = sym.Symbol('h')
    else:
        h = Omega_e[1] - Omega_e[0]
    detJ = h/2  # dx/dX

    for p in range(len(phi)):
        phi[p] = sym.diff(phi[p], X)
    
    for r in range(n):
        for s in range(r, n):
            A_e[r,s] = sym.integrate(phi[r]*phi[s]*detJ, (X, -1, 1))
            A_e[s,r] = A_e[r,s]
    return A_e

def element_vector(f, phi, Omega_e, symbolic=True):
    n = len(phi)
    b_e = sym.zeros(n, 1)
    # Make f a function of X
    X = sym.Symbol('X')
    if symbolic:
        h = sym.Symbol('h')
    else:
        h = Omega_e[1] - Omega_e[0]
    x = (Omega_e[0] + Omega_e[1])/2 + h/2*X  # mapping
    #f = f.subs('x', x)  # substitute mapping formula for x
    detJ = h/2  # dx/d
    for r in range(n):
        b_e[r] = sym.integrate(f*phi[r]*detJ, (X, -1, 1))
    return b_e


def assemble(nodes, elements, phi, f, symbolic=True):
    N_n, N_e = len(nodes), len(elements)
    if symbolic:
        A = sym.zeros(N_n, N_n)
        b = sym.zeros(N_n, 1)
    else:
        A = np.zeros((N_n, N_n))
        b = np.zeros((N_n, 1))



    for e in range(N_e):
        Omega_e = [nodes[elements[e][0]], nodes[elements[e][-1]]]

        A_e = element_matrix(phi, Omega_e, symbolic)
        b_e = element_vector(f, phi, Omega_e, symbolic)

        for r in range(len(elements[e])):
            for s in range(len(elements[e])):
                A[elements[e][r],elements[e][s]] += A_e[r,s]
            b[elements[e][r]] += b_e[r]
    return A, b

def solve(A, b, symbolic=True):
    
    if symbolic:
        c = A.LUsolve(b)           # sympy arrays, symbolic Gaussian elim.
    else:
        c = np.linalg.solve(A, b)  # numpy arrays, numerical solve

    return c

if __name__ == '__main__':


    nodes = [0, 0.5, 1]
    elements = [[0, 1], [1, 2]]
    phi = basis(d=1)
    f = 1
    A, b = assemble(nodes, elements, phi, f, symbolic=True)

    print(A)
    print(b)

    c = solve(A, b, symbolic=True)
    x_arr = np.linspace(0, 1, 100)

#    plt.figure()
#    for i in range(len(c) - 1):
#        plt.plot(x[i:i+1], c[i:i+1], label='numeric')
#    plt.plot(x_arr, u_e(x_arr), label='u_e')
#    plt.show()

