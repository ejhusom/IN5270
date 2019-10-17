#!/usr/bin/env python3
# ============================================================================
# File:     cable_2P1.py
# Author:   Erik Johannes Husom
# Created:  2019-10-17
# ----------------------------------------------------------------------------
# Description:
# Third mandatory exercise in IN5270, University of Oslo, fall 2019.
# Details of the problem is found in cable_2P1.pdf
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt

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
def u_e(x):
    return 0.5 * x**2 - x

def plot(u, x):

    plt.figure()
    plt.plot(x, u, label='u')
    plt.plot(x, u_e(x), label='u_e')
    plt.legend()
    plt.savefig('cable_2P1.png')
    plt.show()

if __name__ == '__main__':


    # Method 1: Excluding the unknown at x = 0
    A = np.array([[4, -2], [-2, 2]])
    b = np.array([-0.5, -0.25])
    c = np.linalg.solve(A, b)
    x = np.linspace(0, 1, 100)
    u = 0*phi0(x) + c[0]*phi1(x) + c[1]*phi2(x)

    plot(u, x)

    # Method 2: Modifying the linear system:
    A = np.array([[1, 0, 0], [-2, 4, -2], [0, -2, 2]])
    b = np.array([0, -0.5, -0.25])
    c = np.linalg.solve(A, b)
    x = np.linspace(0, 1, 100)
    u = c[0]*phi0(x) + c[1]*phi1(x) + c[2]*phi2(x)

    plot(u, x)

