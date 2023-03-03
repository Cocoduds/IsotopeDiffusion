# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 11:29:59 2023

@author: Colton and Nikhil
"""
# =============================================================================
# X direction is along the grain boundry (oxygen to aluminum or vice versa)
# Y direction is across the grain boundry/bulk
# =============================================================================

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#%% DIFFERENCE SCHEMES

def diffusion(x, h, D): #diffusion in x direction
    uiPlus1 = np.copy(x)
    uiMinus1 = np.copy(x)
    for y in range(len(x[:,0])):
        uiPlus1[y,:] = np.append(x[y,1:],x[y,len(x)-1])
        uiMinus1[y,:] = np.append(100,x[y,:-1])
    return (np.multiply(D,(uiPlus1 + -2*x + uiMinus1) / h**2))

def periodicdiffusion(x,h,D): #diffusion in y direction
    uiPlus1 = np.copy(x)
    uiMinus1 = np.copy(x)
    for y in range(len(x[0,:])):
        uiPlus1[:,y] = np.append(x[1:,y],x[:1,y])
        uiMinus1[:,y] = np.append(x[-1:,y],x[:-1,y])
    return (np.multiply(D,(uiPlus1 + -2*x + uiMinus1) / h**2))


#%% RUNGE KUTTA METHODS

def RK4(x, dt, dx, D):
    fa = diffusion(x, dx, D)
    fb = diffusion(x + fa*dt/2, dx, D)
    fc = diffusion(x + fb*dt/2, dx, D)
    fd = diffusion(x + fc*dt, dx, D)
    return x + 1/6 *(fa + 2*fb + 2*fc + fd) * dt

def periodicRK4(x,dt,dx,D):
    fa = periodicdiffusion(x, dx, D)
    fb = periodicdiffusion(x + fa*dt/2, dx, D)
    fc = periodicdiffusion(x + fb*dt/2, dx, D)
    fd = periodicdiffusion(x + fc*dt, dx, D)
    return x + 1/6 *(fa + 2*fb + 2*fc + fd) * dt


#%% Initializing

def initial_left(x): #unsaturated oxygen
    s = (int((xmax-xmin)/h),int((xmax-xmin)/h))
    y = np.zeros(s)
    y[:,5] = 100
    return(y)

#Setting step size

h=0.1
xmin=0
xmax=1
dt=0.001
s=(int((xmax-xmin)/h),int((xmax-xmin)/h))
x = np.linspace(int(xmin),int(xmax),int((xmax-xmin)/h)) #this just defines the axis for plots
y = np.linspace(int(xmin),int(xmax),int((xmax-xmin)/h)) #this is only used in 3d plotting

# MATRIX OF DIFFUSION COEFICIENTS
D = np.ones(s)
for i in range (1,int((xmax-xmin)/h)):
    D[i,:] = 0.01

#Setting the starting conditions
u = initial_left(x)
#v = initial_right(x)


print('x=',x)
plt.figure(1)
plt.plot(x,u)
print('x=', x)
print('y=', u)
plt.show()


#%% THIS LOOP IS FOR 1D only
# for i in range (0, 20000000):
#     u = RK4(u, dt, h, D, 1)
#     # v = RK4(v, dt, h, D/2, 2)
#     maxes.append(np.max(u))
#     times.append(i*dt)
#     if i%500 == 0:
#         # plt.plot(x,v, label = 't = '+str(i*0.001))
#         plt.plot(x,u, label = 't = '+str(i*0.001))
#         plt.xlim(-1,20)
#         plt.ylim(0,100)
#         plt.axvline(x = 0, color = 'b')
#         plt.axvline(x = xmax, color = 'b')
#         plt.show()
#     if u[len(u)-1] > 30:
#         u = np.append(u,0)
#         # v = np.append(0,v)
#         xmax = xmax + h
#         x = np.append(x,xmax+h)
#     # if v[0] > 25:
#     #     u = np.append(u,0)
#     #     v = np.append(0,v)
#     #     xmax = xmax + h
#     #     x = np.append(x,xmax+h)
    

#%%This loop is for 3D

for i in range (0,100000):
    u = RK4(u, dt, h, D)
    u = periodicRK4(u, dt, h, D)
    u[:,0] = 100
    
    # if i%50 == 0: THIS IS FOR A 3D PLOT
    #     fig = plt.figure()
    #     ax = plt.axes(projection='3d')
    #     ax.contour3D(x, y, u, 50, cmap='coolwarm')
    #     ax.set_xlabel('x')
    #     ax.set_ylabel('y')
    #     ax.set_zlabel('z')
    #     ax.set_title('t='+ str(i))
    #     plt.show()
    
    if i%5000 == 0: #2D projections
        plt.plot(x, u[0,:])
        plt.show()
        
    
    
    if np.sum(u[:,len(u[0,:])-1]) > 30: #Moving boundry
        u = np.c_[u,np.zeros(len(u[:,0]))]
        xmax = xmax + h
        x = np.append(x,xmax+h)
        D = np.c_[D,[1,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]]

        