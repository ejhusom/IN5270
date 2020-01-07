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
import sympy as sp
from scipy.integrate import quad
import matplotlib.pyplot as plt

# def p1_elements():

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

f = lambda x: 1 + 2*x -x**2

A1 = np.zeros((2, 2))
A2 = np.zeros((2, 2))
A = np.zeros((3, 3))
b = np.zeros(3)

all_phi = [phi0, phi1, phi2]
phi_e1 = [phi0, phi1]
phi_e2 = [phi1, phi2]


for i in range(2):
    for j in range(2):
        func1 = lambda x: phi_e1[i](x) * phi_e1[j](x)
        func2 = lambda x: phi_e2[i](x) * phi_e2[j](x)
        A1[i, j] = quad(func1, 0, 0.5)[0]
        A2[i, j] = quad(func2, 0.5, 1)[0]

xs = [0, 0.5, 1]


for i in range(3):
    func = lambda x: all_phi[i](x) * f(x)
    b[i] = quad(func, 0, 1)[0]


A[0:2,0:2] = A1
A[1:3, 1:3] += A2


A = np.array([[2/15, 1/15, -1/30], [1/15, 8/15, 1/15], [-1/30, 1/15, 2/15]])

b = np.array([11/60, 17/15, 7/20])
b = np.array([-0.483333333333, 1.13333333, 7/20])

def phi0(x):
    return (1-2*x)*(1-x)

def phi0(x):
    return 2*x*(2-2*x)

def phi2(x):
    return x*(2*x-1)

# c = [1, 7/4, 2]
c = np.linalg.solve(A, b)


print('A:', A)
print('b:', b)
print('c:', c)

x = np.linspace(0, 1, 100)

u = np.zeros(len(x))


u = c[0]*phi0(x) + c[1]*phi1(x) + c[2]*phi2(x)


plt.plot(x, u, label='approx')
plt.plot(x, f(x), label='exact')
plt.legend()
plt.show()
