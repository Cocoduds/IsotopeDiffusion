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
    for i in range(len(x[:,0])):
        uiPlus1[i,:] = np.append(x[i,1:],x[i,len(x[0,:])-1])
        uiMinus1[i,:] = np.append(100,x[i,:-1])
    return (np.multiply(D,(uiPlus1 + -2*x + uiMinus1) / h**2))

def periodicdiffusion(x, h, D): #diffusion in y direction
    uiPlus1 = np.copy(x)
    uiMinus1 = np.copy(x)
    for i in range(len(x[0,:])):
        uiPlus1[:,i] = np.append(x[1:,i],x[:1,i])
        uiMinus1[:,i] = np.append(x[-1:,i],x[:-1,i])
    return (np.multiply(D,(uiPlus1 + -2*x + uiMinus1) / h**2))

def mirrordiffusion(x, h, D): #for y direction half boundry
    uiPlus1 = np.copy(x)
    uiMinus1 = np.copy(x)
    for i in range(len(x[0,:])):
        uiPlus1[:,i] = np.append(x[1:,i],x[len(x[:,0])-1,i])
        uiMinus1[:,i] = np.append(x[0,i],x[:-1,i])
    return (np.multiply(D,(uiPlus1 + -2*x + uiMinus1) / h**2))



#%% RUNGE KUTTA METHODS

def RK4(x, dt, dx, D):
    fa = diffusion(x, dx, D)
    fb = diffusion(x + fa*dt/2, dx, D)
    fc = diffusion(x + fb*dt/2, dx, D)
    fd = diffusion(x + fc*dt, dx, D)
    return x + 1/6 *(fa + 2*fb + 2*fc + fd) * dt

def periodicRK4(x, dt, dx, D):
    fa = periodicdiffusion(x, dx, D)
    fb = periodicdiffusion(x + fa*dt/2, dx, D)
    fc = periodicdiffusion(x + fb*dt/2, dx, D)
    fd = periodicdiffusion(x + fc*dt, dx, D)
    return x + 1/6 *(fa + 2*fb + 2*fc + fd) * dt

def mirrorRK4(x, dt, dx, D):
    fa = mirrordiffusion(x, dx, D)
    fb = mirrordiffusion(x + fa*dt/2, dx, D)
    fc = mirrordiffusion(x + fb*dt/2, dx, D)
    fd = mirrordiffusion(x + fc*dt, dx, D)
    return x + 1/6 *(fa + 2*fb + 2*fc + fd) * dt

#EULER METHODS

def euler(x, dt, dx, D):
    return((dt * diffusion(x, dx, D)))

def mirroreuler(x, dt, dx, D):
    return((dt * mirrordiffusion(x, dx, D)))

def periodiceuler(yn, stepsize, t, h):
    return(yn + (stepsize * periodicdiffusion(t, yn, h)))


#%% Initializing

def initial_left(x): #unsaturated oxygen
    stemp = (int((ymax-ymin)/h),int((xmax-xmin)/h))
    ytemp = np.zeros(stemp)
    return(ytemp)

#Setting step size

h=0.1 #SPATIAL STEP SIZE
xmin=0 
xmax=1 #INITIAL RANGE BETWEEN OXYGEN AND ALLUMINIUM
ymin=0
ymax=0.5 #DISTANCE IN BULK DIRECTION
dt=0.001 #TEMPORAL STEP SIZE

s=(int((ymax-ymin)/h),int((xmax-xmin)/h)) #number of nodes (used in defining D matrix)
x = np.linspace(int(xmin),int(xmax),int((xmax-xmin)/h)) #this just defines the axis for plots
y = np.linspace(int(ymin),int(ymax),int((ymax-ymin)/h)) #this is only used in 3d plotting

# MATRIX OF DIFFUSION COEFICIENTS
D = np.ones(s)
for i in range (1,int((ymax-ymin)/h)):
    D[i,:] = 0.01

#Setting the starting conditions
u = initial_left(x)
#v = initial_right(x)

tracer=[]

print('x=',x)
plt.figure(1)
plt.plot(x,u[0,:])
print('y=', u[0,:])
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
    u[:,0] = 100
    u = u + RK4(u, dt, h, D) + mirrorRK4(u, dt, h, D)
    
    # if i%50 == 0: THIS IS FOR A 3D PLOT
    #     fig = plt.figure()
    #     ax = plt.axes(projection='3d')
    #     ax.contour3D(x, y, u, 50, cmap='coolwarm')
    #     ax.set_xlabel('x')
    #     ax.set_ylabel('y')
    #     ax.set_zlabel('z')
    #     ax.set_title('t='+ str(i))
    #     plt.show()
    
    if i%500 == 0: #2D projections
        plt.plot(x, u[0,:], label = 'Grain boundry')
        plt.plot(x, u[int((ymax-ymin)/h)-1,:], label = 'Middle of bulk')
        plt.legend()
        plt.xlabel('Oxygen ----> Aluminum')
        plt.ylabel('02 Concentration')
        plt.title(str(i)+' time steps')
        plt.ylim(0,100)
        plt.xlim(0,1.5)
        plt.show()
        # plt.plot([1,2,3,4,5],u[:,1])
        # plt.plot([1,2,3,4,5],u[:,2])
        # plt.plot([1,2,3,4,5],u[:,3])
        # plt.plot([1,2,3,4,5],u[:,4])
        # plt.plot([1,2,3,4,5],u[:,0])
        # plt.show()
        

    # if i%200000 == 0: #2D projections
        # plt.close()
        # plt.plot(x, oTotal[0,:], label = 'Total O gb')
        # print(oTotal[0,:])
        # plt.plot(x,o18[0,:], label = 'O18% gb')
        # plt.plot(x,o18[2,:], label = 'O18% bulk')
        # plt.xlabel('concentration')
        # plt.axvline(x = 0, color = 'b')
        # plt.axvline(x = xmax, color = 'b')
        # # plt.ylim(0,1.05)
        # plt.legend()
        # plt.show()
    
    if np.sum(oTotal[:,len(oTotal[0,:])-1]) > 3-(3*10e-2): #Moving boundry
        # plt.close()
        # plt.plot(x, oTotal[0,:], label = 'Total O gb')
        # plt.plot(x,o18[0,:], label = 'O18% gb')
        # plt.plot(x,o18[2,:], label = 'O18% bulk')
        # plt.xlabel('concentration')
        # plt.axvline(x = 0, color = 'b')
        # plt.axvline(x = xmax, color = 'b')
        # plt.ylim(0,1.05)
        # plt.legend()
        # plt.show()
    
        oTotal = np.c_[oTotal,(100-10e-4)*np.zeros(len(oTotal[:,0]))]
        o18 = np.c_[o18,o18[:,(len(o18[:,0]))]]
        xmax = xmax + h
        x = np.append(x,xmax)
        D = np.c_[D,[1,0.001,0.001]]
        
        print('O18 Concentration is: ', o18[0, (len(o18)-1,)], ' at GB and' , o18[2, (len(o18)-1,)], ' in bulk')
        tracer.append(o18[0, (len(o18)-1,)])
        
        
#%%
plt.close()
plt.plot(np.linspace(0,len(tracer), len(tracer)), tracer)
plt.ylabel('concentration of O18')
plt.xlabel('node #')
plt.title('concentration of tracer in new nodes')
plt.savefig('2d.png')
