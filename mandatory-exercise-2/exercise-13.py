#!/usr/bin/env python3
# ==================================================================
# File:     exercise-13.py
# Author:   Erik Johannes Husom
# Created:  2019-02-09
# ------------------------------------------------------------------
# Description:
# Finite difference method for wave motion.
# ==================================================================
import numpy as np
import matplotlib.pyplot as plt

def solver(L, T, Nx, Nt):
    '''Numerically solving second-order ODE.'''
    dx = L/N
    x = np.linspace(0, L, Nx+1)
    dt = T/N
    t = np.linspace(0, T, Nt+1)
    u = np.zeros(Nx+1)      # u[n-1]
    u_1[i] = np.zeros(Nx+1) # u[n]
    u_2[i] = np.zeros(Nx+1) # u[n+1]
    
    for i in range(0, Nx+1):
        u_1[i] = 0


    for n in range(1, N):
        u[n] =
        
    return u, t

# Example parameters and initial conditions
Nx = 10
L = 1
Nt = 10
T = 1
# Quadratic source term
f_lin = lambda t: w**2*(V*t + I)
f_quad = lambda t: 2*K + w**2*(K*t**2 + V*t + I)

u_lin, t_lin = solver(I, V, w, f_lin, T, N)
u_quad, t_quad = solver(I, V, w, f_quad, T, N)

plt.plot(t_lin, u_lin, '-', label='numerical linear')
plt.plot(t_quad, u_quad, '-', label='numerical quadratic')
plt.legend()
plt.show()


