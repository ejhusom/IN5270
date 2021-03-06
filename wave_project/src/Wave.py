#!/usr/bin/env python3
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
        self.x = np.linspace(0, self.Lx, self.Nx)
        self.y = np.linspace(0, self.Ly, self.Ny)

        self.dt2 = self.dt**2
        self.dtdx2 = self.dt2/(2*self.dx**2)
        self.dtdy2 = self.dt2/(2*self.dy**2)

        # Set initial conditions and advance
        if version=='vectorized':
            self.u[:,:,0] = self.I(self.x, self.y)
            self.u[:,:,-1] = self.u[:,:,0] - self.dt*self.V(self.x, self.y)

            self.advance_vec()

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

        for n in range(0, self.Nt-1):
            self.scheme_vec(n)

            #boundary conditions
            #j = 0;
            self.scheme_vec2(n, jstart = 0, jstop = 1, jpstart = 1, jpstop = 2, jmstart = 1, jmstop = 2)
            
            #j = Ny-1
            self.scheme_vec2(n, jstart = -1, jstop = None, jpstart = -2, jpstop = -1, jmstart = -2, jmstop = -1)

            
            #i = 0
            self.scheme_vec(n, istart = 0, istop = 1, ipstart = 1, ipstop = 2, imstart = 1, imstop = 2)
                   
            #i = Nx-1
            self.scheme_vec(n, istart = -1, istop = None, ipstart = -2, ipstop = -1, imstart = -2, imstop = -1)
            
            #Create corners
            #i = 0; j = 0
            self.scheme_vec(n, istart = 0, istop = 1, ipstart = 1, ipstop = 2, imstart = 1, imstop = 2,\
                          jstart = 0, jstop = 1, jpstart = 1, jpstop = 2, jmstart = 1, jmstop = 2)
                   
            #i = 0; j = Ny-1
            self.scheme_vec(n,
                       istart = 0, istop = 1, ipstart = 1, ipstop = 2, imstart = 1, imstop = 2,\
                       jstart = -1, jstop = None, jpstart = -2, jpstop = -1, jmstart = -2, jmstop = -1)
                           
            #i = Nx-1; j = 0
            self.scheme_vec(n, 
                       istart = -1, istop = None, ipstart = -2, ipstop = -1, imstart = -2, imstop = -1,\
                       jstart = 0, jstop = 1, jpstart = 1, jpstop = 2, jmstart = 1, jmstop = 2)
                           
                   
            #i = Nx-1; j = Ny-1
            self.scheme_vec(n,
                       istart = -1, istop = None, ipstart = -2, ipstop = -1, imstart = -2, imstop = -1,\
                       jstart = -1, jstop = None, jpstart = -2, jpstop = -1, jmstart = -2, jmstop = -1)


    def scheme_scalar(self, n, i, ip, im, j, jp, jm):

        u = self.u
        q = self.q
        x = self.x
        y = self.y

        u[i,j,n+1] = 2*u[i,j,n] + (0.5*self.b*self.dt - 1)*u[i,j,n-1] + \
            self.dtdx2*((q(x[i],y[j]) + q(x[ip],y[j]))*(u[ip,j,n] - u[i,j,n])-\
            (q(x[i],y[j]) + q(x[im],y[j]))*(u[i,j,n]-u[im,j,n])) + \
            self.dtdy2*((q(x[i],y[jp]) + q(x[i],y[j]))*(u[i,jp,n] - u[i,j,n])- \
            (q(x[i],y[j]) + q(x[i],y[jm]))*(u[i,j,n] - u[i,jm,n])) + \
            self.dt2*self.f(x[i],y[j],self.dt*n)
        
        u[i,j,n+1] /= 1 + 0.5*self.b*self.dt


    def scheme_vec(self, n,
            istart=1, istop=-1, ipstart=2, ipstop=None, imstart=0, imstop=-2,\
            jstart=1, jstop=-1, jpstart=2, jpstop=None, jmstart=0, jmstop=-2):

        
        u = self.u
        q = self.q
        x = self.x
        y = self.y

        u[istart:istop,jstart:jstop,n+1] = \
        2*u[istart:istop,jstart:jstop,n] - (1 - 0.5*self.b*self.dt)*u[istart:istop,jstart:jstop,n-1] + \
        self.dtdx2*((q(x[ipstart:ipstop],y[jstart:jstop]) + q(x[istart:istop],y[jstart:jstop]))\
               *(u[ipstart:ipstop,jstart:jstop,n] - u[istart:istop,jstart:jstop,n]) \
               - (q(x[istart:istop],y[jstart:jstop]) + q(x[imstart:imstop],y[jstart:jstop]))\
               *(u[istart:istop,jstart:jstop,n] - u[imstart:imstop,jstart:jstop,n])) + \
        self.dtdy2*((q(x[istart:istop],y[jpstart:jpstop]) + q(x[istart:istop],y[jstart:jstop]))\
               *(u[istart:istop,jpstart:jpstop,n] - u[istart:istop,jstart:jstop,n]) \
               - (q(x[istart:istop],y[jstart:jstop]) + q(x[istart:istop],y[jmstart:jmstop]))\
               *(u[istart:istop,jstart:jstop,n] -u[istart:istop,jmstart:jmstop,n])) + \
        self.dt2*self.f(x[istart:istop],y[jstart:jstop], self.dt*n)
        
        u[istart:istop,jstart:jstop,n+1] /=  1 + 0.5*self.b*self.dt
        
        
    def scheme_vec2(self, n,\
                   istart = 1, istop = -1, ipstart = 2, ipstop = None, imstart = 0, imstop = -2,\
                   jstart = 1, jstop = -1, jpstart = 2, jpstop = None, jmstart = 0, jmstop = -2):

             
        u = self.u
        q = self.q
        x = self.x
        y = self.y
                   
        u[istart:istop,jstart:jstop,n+1] = \
        2*u[istart:istop,jstart:jstop,n] - (1 - 0.5*self.b*self.dt)*u[istart:istop,jstart:jstop,n-1] + \
        self.dtdx2*((q(x[ipstart:ipstop],y[jstart:jstop]) + q(x[istart:istop],y[jstart:jstop])).reshape(-1,1)\
               *(u[ipstart:ipstop,jstart:jstop,n] - u[istart:istop,jstart:jstop,n]) \
               - (q(x[istart:istop],y[jstart:jstop]) + q(x[imstart:imstop],y[jstart:jstop])).reshape(-1,1)\
               *(u[istart:istop,jstart:jstop,n] - u[imstart:imstop,jstart:jstop,n])) + \
        self.dtdy2*((q(x[istart:istop],y[jpstart:jpstop]) + q(x[istart:istop],y[jstart:jstop])).reshape(-1,1)\
               *(u[istart:istop,jpstart:jpstop,n] - u[istart:istop,jstart:jstop,n]) \
                - (q(x[istart:istop],y[jstart:jstop]) + q(x[istart:istop],y[jmstart:jmstop])).reshape(-1,1)\
                *(u[istart:istop,jstart:jstop,n] -u[istart:istop,jmstart:jmstop,n])) + \
                self.dt2*(self.f(x[istart:istop],y[jstart:jstop], self.dt*n)).reshape(-1,1)
        u[istart:istop,jstart:jstop,n+1] /=  1 + 0.5*self.b*self.dt

if __name__ == '__main__':
    pass
