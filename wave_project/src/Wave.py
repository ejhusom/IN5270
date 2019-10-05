#!/usr/bin/env python3
# ============================================================================
# File:     wave.py
# Author:   Erik Johannes Husom
# Created:  2019-10-04
# ----------------------------------------------------------------------------
# Description:
# Solving 2D linear wave equation with damping and variable coefficients.
# ============================================================================
import numpy as np



class Wave():




    def __init__(self, I, V, q, f, b, Lx, dx, Ly, dy, T, dt, version='scalar'):


        self.I = I
        self.V = V
        self.q = q
        self.f = f
        self.b = b
        self.Lx, self.dx, self.Ly, self.dy = Lx, dx, Ly, dy
        self.dt, self.T = dt, T
        self.version = version

        self.Nx = int(round(Lx/float(dx)))
        self.Ny = int(round(Ly/float(dy)))
        self.Nt = int(round(T/float(dt)))

        self.t = np.linspace(0, self.Nt*self.dt, self.Nt+1)

        self.u = np.zeros((self.Nx, self.Ny, self.Nt))
        self.x = np.linspace(0, self.Lx, self.Nx+1)
        self.y = np.linspace(0, self.Ly, self.Ny+1)

        self.dt2 = self.dt**2
        self.dtdx2 = self.dt2/(2*self.dx**2)
        self.dtdy2 = self.dt2/(2*self.dy**2)

        # Set initial conditions and advance
        if version=='vectorized':
            self.u[:,:,0] = self.I(self.x, self.y)
            self.u[:,:,-1] = self.u[:,:,0] - self.dt*self.V(self.x, self.y)
        else:
            for i in range(0, self.Nx):
                for j in range(0, self.Ny):
                    self.u[i,j,0] = self.I(self.x[i],self.y[j])

            for i in range(0, self.Nx):
                for j in range(0, self.Ny):
                    self.u[i,j,-1] = self.u[i,j,0] - self.dt*self.V(self.x[i],self.y[j])
            
            self.advance_scalar()
        



    def advance_scalar(self):
        

        for n in range(0, self.Nt-1):
            for i in range(1, self.Nx-1):
                for j in range(1, self.Ny-1):
                    self.scheme_scalar(n, i, i+1, i-1, j, j+1, j-1)

            for i in range(1, self.Nx-1):
                j = 0
                self.scheme_scalar(n, i, i+1, i-1, j, j+1, j+1)

                j = self.Ny-1
                self.scheme_scalar(n, i, i+1, i-1, j, j-1, j-1)

            for j in range(1, self.Ny-1):
                i = 0
                self.scheme_scalar(n, i, i+1, i+1, j, j+1, j-1)

                i = self.Nx-1
                self.scheme_scalar(n, i, i-1, i-1, j, j+1, j-1)

            i = 0
            j = 0
            self.scheme_scalar(n, i, i+1, i+1, j, j+1, j+1)
            i = 0
            j = self.Ny-1
            self.scheme_scalar(n, i, i+1, i+1, j, j-1, j-1)
            i = self.Nx-1
            j = 0
            self.scheme_scalar(n, i, i-1, i-1, j, j+1, j+1)
            i = self.Nx-1
            j = self.Ny-1
            self.scheme_scalar(n, i, i-1, i-1, j, j-1, j-1)



    def advance_vec(self):
        pass

    def scheme_scalar(self, n, i, ip, im, j, jp, jm):

        u = self.u
        q = self.q
        f = self.f
        x = self.x
        y = self.y

        u[i,j,n+1] = 2*u[i,j,n] + (0.5*self.b*self.dt - 1)*u[i,j,n-1] + \
            self.dtdx2*((q(x[i],y[j]) + q(x[ip],y[j]))*(u[ip,j,n] - u[i,j,n])-\
            (q(x[i],y[j]) + q(x[im],y[j]))*(u[i,j,n]-u[im,j,n])) + \
            self.dtdy2*((q(x[i],y[jp]) + q(x[i],y[j]))*(u[i,jp,n] - u[i,j,n])- \
            (q(x[i],y[j]) + q(x[i],y[jm]))*(u[i,j,n] - u[i,jm,n])) + \
            self.dt2*f(x[i],y[j],self.dt*n)
        
        u[i,j,n+1] /= 1 + 0.5*self.b*self.dt


    def scheme_vec(self, n, istart, istop, ipstart, ipstop, imstart, imstop,\
                            jstart, jstop, jpstart, jpstop, jmstart, jmstop):


        u[istart:istop,jstart:jstop,n+1] = \
        2*u[istart:istop,jstart:jstop,n] - (1 - 0.5*b*dt)*u[istart:istop,jstart:jstop,n-1] + \
        dtdx2*((q(x[i2start:i2stop],y[jstart:jstop]) + q(x[istart:istop],y[jstart:jstop]))\
               *(u[i2start:i2stop,jstart:jstop,n] - u[istart:istop,jstart:jstop,n]) \
               - (q(x[istart:istop],y[jstart:jstop]) + q(x[i3start:i3stop],y[jstart:jstop]))\
               *(u[istart:istop,jstart:jstop,n] - u[i3start:i3stop,jstart:jstop,n])) + \
        dtdy2*((q(x[istart:istop],y[j2start:j2stop]) + q(x[istart:istop],y[jstart:jstop]))\
               *(u[istart:istop,j2start:j2stop,n] - u[istart:istop,jstart:jstop,n]) \
               - (q(x[istart:istop],y[jstart:jstop]) + q(x[istart:istop],y[j3start:j3stop]))\
               *(u[istart:istop,jstart:jstop,n] -u[istart:istop,j3start:j3stop,n])) + \
        dt2*f(x[istart:istop],y[jstart:jstop], dt*n)
        u[istart:istop,jstart:jstop,n+1] /=  1 + 0.5*b*dt
        


if __name__ == '__main__':
    pass
