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
<<<<<<< Updated upstream
        uiPlus1[i,:] = np.append(x[i,1:],x[i,len(x[0,:])-1])
=======
        uiPlus1[i,:] = np.append(x[i,1:],x[i,len(x[0,:])-1]) # THIS IS FOR NO DIFFUSION INTO METAL
        #uiPlus1[i,:] = np.append(x[i,1:],metalBoundryTotal/len(x[:,0])) #this boundry condiditon is eq. 23 (unweighted average atm)
>>>>>>> Stashed changes
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


<<<<<<< Updated upstream
=======
def isotopeDiffusion(x, co, h, Do, L, fl, Omega, cop):
    uiPlus1 = np.copy(x)
    uiMinus1 = np.copy(x)
    coPlus1 = np.copy(co)
    coMinus1 = np.copy(co)
    metalBoundryTotal = np.sum(u[:,len(u[0,:])-1]) #addition of concentrations at metal boundry
    for i in range(len(x[:,0])):
        #uiPlus1[i,:] = np.append(x[i,1:],x[i,len(x[0,:])-1]) # THIS IS FOR NO DIFFUSION INTO METAL
        uiPlus1[i,:] = np.append(x[i,1:],metalBoundryTotal/len(x[:,0])) #this boundry condiditon is eq. 23 (unweighted average atm)
        uiMinus1[i,:] = np.append(100,x[i,:-1])
        coPlus1[i,:] = np.append(co[i,1:],metalBoundryTotal/len(co[:,0])) 
        coMinus1[i,:] = np.append(100,co[i,:-1])
        
    return (np.multiply(Do,(uiPlus1 + -2*x + uiMinus1) / h**2)) 
    - ((1/cop)*np.multiply(co,(np.multiply(Do,(uiPlus1 + -2*x + uiMinus1) / h**2))))
    - ((1/cop))*(np.multiply(Do,(uiPlus1 - uiMinus1) / 2*h))
    - ((3*fl)/(Omega*L))*np.devide((np.multiply(Do,(uiPlus1 - uiMinus1) / 2*h)),co)
    
def mirrorIsoptopeDiffusion(x, co, h, Do, L, fl, Omega, cop):
    uiPlus1 = np.copy(x)
    uiMinus1 = np.copy(x)
    coPlus1 = np.copy(co)
    coMinus1 = np.copy(co)
    metalBoundryTotal = np.sum(u[:,len(u[0,:])-1]) #addition of concentrations at metal boundry
    for i in range(len(x[:,0])):
        #uiPlus1[i,:] = np.append(x[i,1:],x[i,len(x[0,:])-1]) # THIS IS FOR NO DIFFUSION INTO METAL
        uiPlus1[i,:] = np.append(x[i,1:],metalBoundryTotal/len(x[:,0])) #this boundry condiditon is eq. 23 (unweighted average atm)
        uiMinus1[i,:] = np.append(100,x[i,:-1])
        coPlus1[i,:] = np.append(co[i,1:],metalBoundryTotal/len(co[:,0])) 
        coMinus1[i,:] = np.append(100,co[i,:-1])
        
    return (np.multiply(Do,(uiPlus1 + -2*x + uiMinus1) / h**2)) 
    - ((1/cop)*np.multiply(co,(np.multiply(Do,(uiPlus1 + -2*x + uiMinus1) / h**2))))
    - ((1/cop))*(np.multiply(Do,(uiPlus1 - uiMinus1) / 2*h))
    - ((3*fl)/(Omega*L))*np.devide((np.multiply(Do,(uiPlus1 - uiMinus1) / 2*h)),co)
>>>>>>> Stashed changes

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



def isotopeRK4(x, co, h, Do, L, fl, Omega, cop, dt):
    fa = isotopeDiffusion(x, co, h, Do, L, fl, Omega, cop)
    fb = isotopeDiffusion(x + fa*dt/2, co, h, Do, L, fl, Omega, cop)
    fc = isotopeDiffusion(x + fb*dt/2, co, h, Do, L, fl, Omega, cop)
    fd = isotopeDiffusion(x + fc*dt, co, h, Do, L, fl, Omega, cop)
    return 1/6 *(fa + 2*fb + 2*fc + fd) * dt

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

def initalOxygen(x):
    stemp = (int((ymax-ymin)/dx)+1,int((xmax-xmin)/dx)+1)
    y = np.ones(stemp)
    return(y*99)

#Setting step size

h=0.1 #SPATIAL STEP SIZE
xmin=0 
xmax=1 #INITIAL RANGE BETWEEN OXYGEN AND ALLUMINIUM
ymin=0
ymax=0.5 #DISTANCE IN BULK DIRECTION
dt=0.001 #TEMPORAL STEP SIZE

<<<<<<< Updated upstream
s=(int((ymax-ymin)/h),int((xmax-xmin)/h)) #number of nodes (used in defining D matrix)
x = np.linspace(int(xmin),int(xmax),int((xmax-xmin)/h)) #this just defines the axis for plots
y = np.linspace(int(ymin),int(ymax),int((ymax-ymin)/h)) #this is only used in 3d plotting
=======
cop=1
L=1
fl=1
Omega=1

s=(int((ymax-ymin)/dx)+1,int((xmax-xmin)/dx)+1) #number of nodes (used in defining D matrix)
x = np.linspace(int(xmin),int(xmax),int((xmax-xmin)/dx)+1) #this just defines the axis for plots
y = np.linspace(int(xmin),int(xmax),int((xmax-xmin)/dx)+1) #this is only used in 3d plotting
>>>>>>> Stashed changes

# MATRIX OF DIFFUSION COEFICIENTS
D = np.ones(s)
for i in range (1,int((ymax-ymin)/h)):
    D[i,:] = 0.01

#Setting the starting conditions
u = initial_left(x)
o18 = initial_left(x)
oTotal = initalOxygen(x)
#v = initial_right(x)


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
<<<<<<< Updated upstream
for i in range (0,100000):
    u[:,0] = 100
    u = u + RK4(u, dt, h, D) + mirrorRK4(u, dt, h, D)
=======
for i in range (0,200):
    oTotal = oTotal + RK4(oTotal, dt, dx, D)
    o18 = o18 + isotopeRK4(o18, oTotal, dx, D, L, fl, Omega, cop, dt)
    oTotal[:,0] = 100
    o18[:,0] = 100
>>>>>>> Stashed changes
    
    # if i%50 == 0: THIS IS FOR A 3D PLOT
    #     fig = plt.figure()
    #     ax = plt.axes(projection='3d')
    #     ax.contour3D(x, y, u, 50, cmap='coolwarm')
    #     ax.set_xlabel('x')
    #     ax.set_ylabel('y')
    #     ax.set_zlabel('z')
    #     ax.set_title('t='+ str(i))
    #     plt.show()
<<<<<<< Updated upstream
    
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
=======
        
    if i%1 == 0: #2D projections
        plt.plot(x, oTotal[0,:], label = 'total O')
        plt.plot(x,o18[0,:], label = 'O18%')
        plt.xlabel('concentration')
        plt.axvline(x = 0, color = 'b')
        plt.axvline(x = xmax, color = 'b')
        plt.legend()
        plt.show()
    
    if np.sum(u[:,len(oTotal[0,:])-1]) > 598: #Moving boundry
        oTotal = np.c_[oTotal,np.zeros(len(oTotal[:,0]))]
        xmax = xmax + dx
        x = np.append(x,xmax)
        D = np.c_[D,[1,0.01,0.01,0.01,0.01,0.01]]
>>>>>>> Stashed changes
        
    
    if np.sum(u[:,len(u[0,:])-1]) > 400: #MOVING BOUNDRY
        u = np.c_[u,np.zeros(len(u[:,0]))]
        xmax = xmax + h #the added step could vary graphically from initial step
        x = np.append(x,xmax)
        D = np.c_[D,[1,0.01,0.01,0.01,0.01]] #THIS IS STILL HARDCODED
        print(u)

        