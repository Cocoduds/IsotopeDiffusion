# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 11:29:59 2023

@author: dudle
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def dispersion(x, h, D):
    uiPlus1 = np.copy(x)
    uiMinus1 = np.copy(x)
    for y in range(len(x[:,0])):
        uiPlus1[y,:] = np.append(x[y,1:],x[y,len(x)-1])
        uiMinus1[y,:] = np.append(100,x[y,:-1])
    return (np.multiply(D,(uiPlus1 + -2*x + uiMinus1) / h**2))

def periodicdispersion(x,h,D):
    uiPlus1 = np.copy(x)
    uiMinus1 = np.copy(x)
    for y in range(len(x[:,0])):
        uiPlus1[:,y] = np.append(x[1:,y],x[:1,y])
        uiMinus1[:,y] = np.append(x[-1:,y],x[:-1,y])
    return (np.multiply(D,(uiPlus1 + -2*x + uiMinus1) / h**2))

def RK4(x, dt, dx, D):
    fa = dispersion(x, dx, D)
    fb = dispersion(x + fa*dt/2, dx, D)
    fc = dispersion(x + fb*dt/2, dx, D)
    fd = dispersion(x + fc*dt, dx, D)
    return x + 1/6 *(fa + 2*fb + 2*fc + fd) * dt

def periodicRK4(x,dt,dx,D):
    fa = periodicdispersion(x, dx, D)
    fb = periodicdispersion(x + fa*dt/2, dx, D)
    fc = periodicdispersion(x + fb*dt/2, dx, D)
    fd = periodicdispersion(x + fc*dt, dx, D)
    return x + 1/6 *(fa + 2*fb + 2*fc + fd) * dt

def initial_left(x): 
    s = (int((xmax-xmin)/h),int((xmax-xmin)/h))
    y = np.zeros(s)
    y[:,0] = 100
    return(y)


h=0.1
xmin=0
xmax=1
dt=0.001
s=(int((xmax-xmin)/h),int((xmax-xmin)/h))
D = np.ones(s)
for i in range (1,int((xmax-xmin)/h)):
    D[i,:] = 0.01

x = np.linspace(int(xmin),int(xmax),int((xmax-xmin)/h))
xinit = np.linspace(int(xmin),int(xmax),int((xmax-xmin)/h))
print('x=',x)
plt.figure(1)
u = initial_left(x)
#v = initial_right(x)
plt.plot(x,u)
print('x=', x)
print('y=', u)
plt.show()
maxes = []
times=[]

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
    
for i in range (0,100000):
    u = RK4(u, dt, h, D)
    u = periodicRK4(u, dt, h, D)
    maxes.append(np.max(u))
    times.append(i*dt)
    u[:,0] = 100
    if i%50 == 0:
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.contour3D(x, xinit, u, 50, cmap='coolwarm')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title('t='+ str(i))
        plt.show()
    if np.sum(u[:,len(u[:,0])-1]) > 5:
        u = np.c_[u,np.zeros(len(u[:,0]))]
        xmax = xmax + h
        x = np.append(x,xmax+h)
        D = np.c_[D,[1,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]]

        