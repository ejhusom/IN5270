#!/usr/bin/env python3
from Wave import *
import nose.tools as nt
from sympy import *

def test_constant_solution():
    b = 0
    Lx, Ly, T = 1, 1, 1
    dx, dy, dt = 0.1, 0.1, 0.1

    V = lambda x, y: 0
    I = lambda x, y: 8
    q = lambda x, y: 4
    f = lambda x, y, n: 0
    exact_solution = 8

    wave_constant= Wave(I, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="scalar")

    diff = abs(exact_solution - wave_constant.u).max()
    print('Test of constant solution (scalar version):')
    print("Max difference: ", diff)
    nt.assert_almost_equal(diff, 0, places=15)


def test_constant_solution_vec():
    b = 0
    Lx, Ly = 10, 10

    dx, dy = 0.1, 0.1
    dt = 0.01
    T = 1.

    V = lambda x, y: 0
    I = lambda x, y: 8

    def q(x, y):
        return np.array(4)

    def f(x, y, n):
        return np.array(0)


    exact_solution = 8

    wave_constant= Wave(I, V, q, f, b, Lx, dx, Ly, dy, T, dt,
            version="vectorized")

    diff = abs(exact_solution - wave_constant.u).max()
    print('Test of constant solution (vectorized version):')
    print("Max difference: ", diff)
    nt.assert_almost_equal(diff, 0, places=15)



def test_plug():
    b = 0
    Lx, Ly = 1, 1
    dx, dy = 0.05, 0.05
    dt = 0.1
    T = 4.

    V = lambda x, y: 0

    def q(x, y):
        return np.array(4)

    def f(x, y, n):
        return np.array(0)

    def Ix(x, y):
        result = np.zeros((np.size(x),np.size(y)))
        result[:,:] = 1
        result[np.where(abs(x-Lx/2.0) > 0.1),:] = 0
        return result

    def Isx(x,y):
        if abs(x-Lx/2.0) > 0.1:
            return 0
        else:
            return 1


    def Iy(x, y):
        result = np.zeros((np.size(x),np.size(y)))
        result[:,:] = 1
        result[:,np.where(abs(y-Ly/2.0) > 0.1)] = 0
        return result

    def Isy(x,y):
        if abs(y-Ly/2.0) > 0.1:
            return 0
        else:
            return 1
        
 
    wave_x_scalar = Wave(Isx, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="scalar")
    wave_y_scalar = Wave(Isy, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="scalar")
    wave_x_vec = Wave(Ix, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="vectorized")
    wave_y_vec = Wave(Iy, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="vectorized")
    diff = abs(wave_x_scalar.u - wave_x_vec.u).max()
    print('Test of plug (x-direction):')
    print("Max difference: ", diff)
    nt.assert_almost_equal(diff, 0, places=13)
    diff = abs(wave_y_scalar.u - wave_y_vec.u).max()
    print('Test of plug (y-direction):')
    print("Max difference: ", diff)
    nt.assert_almost_equal(diff, 0, places=13)



def test_undampened():
    b = 0
    Lx, Ly, T = 10, 10, 10

    V = lambda x, y: 0
    I = lambda x, y: u_e(x, y, 0)
    q = lambda x, y: 10
    f = lambda x, y, n: 0

    def u_e(x, y, t):
        A = 4
        mx, my = 0.1, 0.1
        kx = mx*np.pi/Lx
        ky = my*np.pi/Ly
        omega = 0.1
        x, y = np.meshgrid(x,y)
        return A*np.cos(kx*x)*np.cos(ky*y)*np.cos(omega*t)
            
        
        
    rate_theoretical = 2
    n = 10
    E = []
    c = 0.1
    h_list = np.linspace(5, 0.5, n)

    for h in h_list: 

        Nx = int(round(Lx/float(h)))
        Ny = int(round(Ly/float(h)))
        x = np.linspace(0, Lx, Nx)
        y = np.linspace(0, Ly, Ny)
        
        wave = Wave(I, V, q, f, b, Lx, h, Ly, h, T, c*h, version="scalar")
        v_e = u_e(x, y, T)
        
        E.append(abs(v_e - wave.u[:,:,-1]).max())

    E = np.array(E)
    rate = np.zeros(n-1)

    for i in range(1, n):
        rate[i-1] = np.log(E[i-1]/E[i])/np.log(h_list[i-1]/h_list[i])

    #print(E/h_list**2)
    print(rate)
    diff = abs(rate_theoretical - rate[-1])
    print('Test of undampened solution:')
    print(f'Difference: {diff}')
    nt.assert_almost_equal(diff, 0, places=1)
 

def manufactured():
    f = Symbol('f')
    x = Symbol('x')
    y = Symbol('y')
    t = Symbol('t')
    A = Symbol('A')
    omega = Symbol('omega')
    b = Symbol('b')
    kx = Symbol('kx')
    ky = Symbol('ky')
    u = A*cos(kx*x)*cos(ky*y)*cos(omega*t)
    u_t = diff(u, t)
    u_tt = diff(u_t, t)
    u_x = diff(u, x)
    u_xx = diff(u_x, x)
    u_y = diff(u, y)
    u_yy = diff(u_y, y)
    f = u_tt - u_xx - u_yy
    print(f)

def test_mms():
    b = 0
    Lx, Ly, T = 10, 10, 10
    
    def V(x, y):
        A = 1
        B = 0
        mx = 1
        my = 1
        kx = mx*p.pi/Lx
        ky = my*p.pi/Ly
        c = b/2. 
        omega = np.sqrt(kx**2*q(x,y) + ky**2*q(x,y) - c**2)
        return (omega*B - c*A)*p.cos(kx*x)*p.cos(ky*y)

    def I(x, y):
        return u_e(x, y, 0)
    
    def q(x, y):
        return x

    def f(x, y, t):
        A = 1
        B = 0
        mx = 1
        my = 1
        kx = mx*p.pi/Lx
        ky = my*p.pi/Ly
        c = b/2. 
        omega = np.sqrt(kx**2*q(x,y) + ky**2*q(x,y) - c**2)

        #The source term, using q = x, and dq/dx = 1, dq/dy = 0
        return (c**2-b*c-omega**2 + x*(kx**2+ky**2)*np.cos(kx*x)*p.cos(ky*y) +\
                kx*np.sin(kx*x)*np.cos(ky*y))*(A*np.cos(omega*t) + B*np.sin(omega*t)) +\
                (b*omega - 2*c*omega)*(-A*np.sin(omega*t) + B*np.cos(omega*t))

  

    def u_e(x,y,t):
        A = 1
        B = 0
        mx = 1
        my = 1
        kx = mx*np.pi/Lx
        ky = my*np.pi/Ly
        c = b/2. 
        omega = np.sqrt(kx**2*q(x,y) + ky**2*q(x,y) - c**2)
        x,y = np.meshgrid(x,y)
        return (A*np.cos(omega*t) +
                B*np.sin(omega*t))*np.cos(kx*x)*np.cos(ky*y)*np.exp(-c*t)

        
    expected_rate = 2
    n = 10
    E = []
    c = 0.1
    h_list = np.linspace(5, 0.05, n)
    #h_list = [0.01]
    for h in h_list: 
        Nx = int(round(Lx/float(h)))
        Ny = int(round(Ly/float(h)))

        x = np.linspace(0, Lx, Nx)
        y = np.linspace(0, Ly, Ny)
        
        u = solver(I, V, q, f, b, Lx, h, Ly, h, T, c*h, version="vectorized")
        v_e = u_e(x,y,T)
        
        E.append(abs(v_e - u[:,:,-1]).max())

    E = np.array(E)
    print(E)
    rate = np.zeros(n-1)
    for i in range(1, n):
        rate[i-1] = np.log(E[i-1]/E[i])/np.log(h_list[i-1]/h_list[i])

    print(rate)
    diff = abs(expected_rate - rate[-1])
    nt.assert_almost_equal(diff, 0 ,places=1)


def physical(h, bottom, I0):
    b = 0.2
    Lx, Ly = 2, 2
    c = 0.1
    dx = h
    dy = h
    dt = c*h
    T = 2

    Nx = int(round(Lx/float(dx)))
    Ny = int(round(Ly/float(dy)))
    Nt = int(round(T/float(dt)))

    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    t = np.linspace(0, T, Nt)
    
    V = lambda x, y: 0

    def f(x, y, n):
        return np.array(0)

    #Initial conditions
    def I(x, y):
        #I0 = 2
        Ia = 2
        Im = 0
        Is = 0.5
        return I0 + Ia*np.exp(-((x - Im)/Is)**2) 


    def I2(x, y):
        #I0 = 2
        Ia = 1
        Im = 0
        Is = 0.2
        return I0 + Ia*np.exp(-((x - Im)/Is)**2-((y.reshape(-1,1) - Im)/Is)**2) 

    #3 different kinds of bottom shapes
    def B1(x,y):
        B0 = 0
        Ba = 2.5
        Bmy = 1
        Bmx = 1
        Bs = 0.4
        b = 1
        return B0 + Ba*np.exp(-((x - Bmx)/Bs)**2-((y - Bmy)/(b*Bs))**2) 


    def B2(x,y):
        B0 = 0
        Ba = 2.5
        Bmy = 1
        Bmx = 1
        Bs = 0.4
        b = 1
        
        index = 0 <  np.sqrt(x**2 + y**2)
        index2 = np.sqrt(x**2 + y**2) <= Bs
        index = index+index2
        
        results =  B0 + Ba*np.cos(np.pi*((x - Bmx)/(2*Bs)))*np.cos(np.pi*(y - Bmy)/(2*Bs)) 

        #Temporary solution to make sure that the boundaries for the box are working correctly. 
        #Fix this once it have been tested
        for i in range(len(x)):
            for j in range(len(y)):
                if (0 > np.sqrt((x[i])**2 + (y[j])**2) >= Bs):
                    results[i,j] = B0
        
        return results
    
    


    
    def H(x,y):
        #I0 = 2
        if (bottom == 1):
            return I0 - B1(x,y)
        elif (bottom == 2):
            return I0 - B2(x,y)

            
    def q(x, y):
        return 9.81*H(x,y)


    wave = Wave(I, V, q, f, b, Lx, dx, Ly, dy, T, dt, version="vectorized")
    

    #Save the arrays
    np.save("u", wave.u)
    if (bottom == 1):
        np.save("h",B1(x,y.reshape(-1,1)))
    elif (bottom == 2):
        np.save("h",B2(x,y.reshape(-1,1)))
        
    np.save("x",x)
    np.save("y",y)
    
    
    


    
if __name__ == '__main__':
    #test_constant_solution()
    #test_constant_solution_vec()
    #test_plug()
    #test_undampened()
    #test_mms()
    physical(0.01,2,4)
    #manufactured()
