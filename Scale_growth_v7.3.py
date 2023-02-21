"""
Created on Tue Jan 31 00:27:22 2023

@author: coltondudley
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from datetime import datetime

version  = 7.3
print('SCALE GROWTH MODEL Version',version)
#os.chdir("/Users/finnis/Documents/Projects/Notes/Scale Growth")
working_directory = os.getcwd()
print('Current working directory =', working_directory)
today = datetime.today()
now = datetime.now()
today = today.strftime('%Y.%m.%d')
now = now.strftime('%H-%M-%S')
output_folder = today + '_' + now
os.mkdir(output_folder)
os.chdir(output_folder)
print('\nImages will be saved in folder ' + output_folder)

#%%
"""
Initialize variables not to be changed in the rest of the code

"""
# Initial scale thickness, the length scale for dimensionless equations, should be 1.0.
# so not necessarily referenced in code.
length_0 = 1.0
# Oxygen diffusion coeffient $D_{\text{O}'$ on the grain boundary as defined by MB
# Should always be 1.0 for dimensionless equations, so not necessarily referenced in code.
mobility_O_GB = 1.0
# Aluminium diffusion coeffient $D_{\text{Al}'$ on the grain boundary as defined by MB
# represented here as its ratio to $D_{\text{O}'$.
mobility_Al_GB = 1.0
# Initialize the number of nodes of each type on the grid along the y-axis
# Number of nodes for growth into gas phase
ngridy_gas_0 = 20
# Number of nodes covering initial scale
ngridy_scale_0 = 51
# Number of nodes for growth into metal phase
ngridy_metal_0 = 20
# Fixed total provision of nodes for  gas-scale-metal. This will not be changed in the code.
ngridy_max = ngridy_scale_0 + ngridy_gas_0 + ngridy_metal_0
# Ratio nu of thickness of grain boundary delta to (delta + d)
nu = 0.01
# f' as defined by MB, or the oxygen vacancy concentration at the oxide-metal interface
f_Ox_GB = 0.0001
# We define here the notionally similar quantity for Al, which would be the Al vacancy
# concentration at the gas-scale interface. This is only used to scale the velocity of the
# gas-scale interface
f_Al_GB = 0.0001
# Initial velocity of scale-metal interface (> 0)
# Note - velocities are measured with respect to the fixed lattice of the scale.
# So we imagine that as the scale grows it pushes the interfaces back in opposite directions.
v_2_0 = nu*f_Ox_GB
# Initial velocity of gas-scale interface (< 0)
v_1_0 = - nu*f_Al_GB*mobility_Al_GB
# Distances along the y-axis for the fixed array of nodes to cover the three phases
# The origin (y = 0) is chosen to be the first node of the initial scale, 

#%%

# and y=1 is the initial position of the last node of the scale.
# These distances are set up here and should not be changed by the subsequent code. 
gridy_max = (np.linspace(0,ngridy_max-1,ngridy_max)-ngridy_gas_0)/(ngridy_scale_0-1)
# Increment of distances, a constant
dy = gridy_max[1]-gridy_max[0]
# Initial excess oxygen isotope.
# The gas is treated as having an isotope fraction of 1.
c_GB_0 = np.zeros_like(gridy_max)
c_GB_0[0:ngridy_gas_0+1] = 1.0
np.set_printoptions(precision=6)
# print('gridy_max =\n',np.array2string(gridy_max))
print('gridy_max = \n',gridy_max)
print('\nc_GB_0 = \n',c_GB_0)
# Initial oxygen vacancy concentration
vac_Ox_0 = np.zeros_like(gridy_max)
for ny in range(ngridy_gas_0+1,ngridy_gas_0+ngridy_scale_0):
    vac_Ox_0[ny] =  f_Ox_GB*(gridy_max[ny]-gridy_max[ngridy_gas_0])
print('\nOxygen vacancy concentration (vacancies per O-lattice site):\n',vac_Ox_0)

#%%

fig = plt.figure()
"""
The initial fraction of isotope

"""
ax = fig.add_subplot()
line_c_GB_0 = ax.plot(gridy_max, c_GB_0)
plt.show()

#%%

    
""" 
DEFINE TRACER  DIFFUSION COEFFICIENT AND ITS DERIVATIVE

"""
# Oxygen tracer diffusion coefficient d_GB and its linear y-derivative dd_GB
d_GB_0 = np.zeros_like(gridy_max)
dd_GB_0 = np.zeros_like(gridy_max)
for ny in range(ngridy_gas_0+1,ngridy_gas_0+ngridy_scale_0):   
    d_GB_0[ny] = vac_Ox_0[ny]
    dd_GB_0[ny] = (vac_Ox_0[ny+1]-vac_Ox_0[ny-1])/(2.0*dy)
# Correct the last value of this gradient, else it would be negative and nonsense.     
dd_GB_0[ngridy_gas_0+ngridy_scale_0-1] = dd_GB_0[ngridy_gas_0+ngridy_scale_0-2]
print('\nTracer diffusion coefficient:\n',d_GB_0)
print('\ny-derivative of tracer diffusion coefficient:\n',dd_GB_0)

#%%

def plot_concentrations(new_page, plot_count,fig_number, first_or_second, \
                        graph_layout,number_plots,\
                        title = 'Scale growth',c_1=c_GB_0, c_2=c_GB_0, \
                        label_1 = 'c_1', label_2='c_2'):
    """
    new_page : if zero, open a new page to plot figures on, if   
    graph-layout : tuple to specify the layout of figures on a page, e.g. (2,3).
    first_or_second : if set to 1, plot just one figure, c_1
                      if set to 2 plot only second figure,  c_2
                      if set to anything else, plot both graphs
    """
    print('\n*** ENTERING plot_concentrations ***\n')
    nfigsx, nfigsy = graph_layout
    max_figs_per_page = nfigsx*nfigsy
    # Count the number of figures so far on the current page
    figs_on_page = plot_count%max_figs_per_page
    print('plot_count = ',plot_count, " figs_on_page = ", figs_on_page," max_figs_per_page = ", max_figs_per_page)
    if figs_on_page == 0:
        new_page = plt.figure(figsize =(6,8))
    ax = new_page.add_subplot(nfigsx,nfigsy,figs_on_page+1)
    if first_or_second == 1:
        line_c_1 = ax.plot(gridy_max, c_1,label=label_1)
    elif first_or_second == 2:
        line_c_2 = ax.plot(gridy_max, c_2,label=label_2)
    else:
        line_c_1 = ax.plot(gridy_max, c_1,label=label_1)
        line_c_2 = ax.plot(gridy_max, c_2,label=label_2)
    plt.legend(loc="lower left")
    plt.xlabel('Depth into scale', fontdict=None, labelpad=None        )
    plt.ylabel('Fraction of isotope', fontdict=None, labelpad=None    )
    plt.title(title)
    plot_count += 1
    figs_on_page +=1
    if number_plots == plot_count:
        new_page.tight_layout()
        plt.savefig('Figure_{}.pdf'.format(plot_count))
        plt.show()
        plt.close()
    elif figs_on_page == max_figs_per_page:
        new_page.tight_layout()
        plt.savefig('Figure_{}.pdf'.format(plot_count))
        plt.show()
        plt.close()
    print('\nplot_count =',plot_count, " figs_on_page =", figs_on_page," max_figs_per_page =", max_figs_per_page)
    print('\n*** LEAVING plot_concentrations ***\n')
    return new_page, plot_count

#%%

def update_ranges(direction, ngridy_gas, ngridy_scale, ngridy_metal,c_GB,vac_Ox,d_GB,dd_GB):
    """
    Increment the numbers of nodes in each region, ngridy_gas, ngridy_scale, ngridy_metal,
    according to the direction of scale growth. 
    The overall number of nodes remains fixed. Possible parameter values:
    direction = +1 : The scale grows by 1 node while the metal retreats by 1 node           
    direction = -1: The scale grows by 1 node while the gas retreats by 1 node 
    Then update the other parameters, c_GB, vac_Ox, d_GB, dd_GB, because they depend on the new length 
    of the scale.
    
    """
    print('\n*** ENTERING update_ranges ***\n')
    print('Previous values of ngridy_gas, ngridy_scale, ngridy_metal',ngridy_gas, ngridy_scale, ngridy_metal)
    if direction != 1 and direction != -1:
        print('Argument of update_ranges = ', direction)
        sys.exit('Argument of update_ranges must be -1 or +1')
    if direction == 1:
        c_GB[ngridy_gas+ngridy_scale] = c_GB[ngridy_gas+ngridy_scale-1]
        ngridy_metal = ngridy_metal-1
        if ngridy_metal == 0:
            raise ValueError('RUN OUT OF METAL NODES! Increase ngridy_metal_0 and try again')
        print('One node to be added to scale grid and taken from metal grid')
    if direction == -1:
        ngridy_gas = ngridy_gas-1
        if ngridy_gas == 0:
            raise ValueError('RUN OUT OF GAS NODES! Increase ngridy_gas_0 and try again')
        print('One node to be added to scale grid and taken from gas grid')  
    ngridy_scale = ngridy_scale+1
    for ny in range(ngridy_gas + 1,ngridy_gas + ngridy_scale):
        vac_Ox[ny] = f_Ox_GB*(gridy_max[ny]-gridy_max[ngridy_gas])*((ngridy_scale_0-1)/(ngridy_scale-1))
    for ny in range(ngridy_gas + 1,ngridy_gas + ngridy_scale):   
        d_GB[ny] = vac_Ox[ny]
        dd_GB[ny] = (vac_Ox[ny+1]-vac_Ox[ny-1])/(2.0*dy)
# Correct the last value of this gradient, else it would be negative and nonsense.     
    dd_GB[ngridy_gas+ngridy_scale-1] = dd_GB[ngridy_gas+ngridy_scale-2]
#   print('vac_Ox :\n',vac_Ox,'\n')
#   print('c_GB :\n',c_GB,'\n')
    print('New values of ngridy_gas, ngridy_scale, ngridy_metal = ',ngridy_gas, ngridy_scale, ngridy_metal)
#   print('\nTracer diffusion coefficient:\n',d_GB)
#   print('\ny-derivative of tracer diffusion coefficient:\n',dd_GB,'\n')
    if direction==1:
        print('\n*** LEAVING update_ranges, scale node added in metal direction *** \n')
    if direction==-1:
        print('\n*** LEAVING update_ranges, scale node added in gas direction *** \n')
    return ngridy_gas, ngridy_scale, ngridy_metal, c_GB, vac_Ox, d_GB, dd_GB   

#%%

def tracer_diffusion(timesteps,dt, ngridy_gas_0, ngridy_scale_0, ngridy_metal_0, vac_Ox_0, \
                     c_GB_0, d_GB_0, dd_GB_0,length_0, graph_layout, bleed=0,number_plots=9):
    """  
    Test simple tracer diffusion, ideal solution, vacancy mechanism. Dimensionless units.
    All locally varying quantities are defined on a fixed array, gridy_max[], along the y-axis. 
    timesteps    : maximum number of timesteps to iterate
    dt           : the length of a timestep
    ngridy_gas   : number of nodes in the gas (on the y axis)
    ngridy_scale : number of nodes in the scale    
    ngridy_metal : number of nodes in the metal
    vac_Ox[]     : concentration per site of oxygen vacancies
    c_GB[]       : local grain-boundary ratio of isotope to total oxygen
    d_GB[]       : local oxygen tracer diffusion coefficient
    dd_GB[]      : y-derivative of d_GB[] 
    bleed        : parameter to describe rate of bleeding of c_GB into bulk according to -f(c_GB - c_bulk)
    graph_layout : tuple showing layout of figures on a page, e.g. for 6 figures (2,3) or (3,2).]
    number_plots : the total number of snapshots to plot, including the final iteration. They are equally
                   spaced in time by the value associated with time_stride. 
    """
    print('\n*** STARTING tracer_diffusion ***\n')
# Dimensionless oxygen current J_\text{O}_GB, depends only on time-dependant length
# Note that j_O_GB = f_Ox_GB/length 
# length should start with value 1, but will be incremented in the time loop.
    length = length_0
    ngridy_gas = ngridy_gas_0
    ngridy_scale = ngridy_scale_0
    ngridy_metal = ngridy_metal_0
    c_GB = np.copy(c_GB_0)
    c_bulk = np.copy(c_GB_0)
    vac_Ox = np.copy(vac_Ox_0)
    d_GB = np.copy(d_GB_0)
    dd_GB = np.copy(dd_GB_0) 
    print('timesteps = ',timesteps,'  dt = ',dt,'\ndy =',dy,'\n')
    print('ngridy_gas =',ngridy_gas)
    print('ngridy_scale =',ngridy_scale)
    print('ngridy_gas =',ngridy_metal)
# The inflation factor must be <1 for convergence of the Euler method. 
    print('Inflation factor = ', 2*dt*f_Ox_GB/(dy*dy))
    print('A near-optimal timestep would be ', dy*dy/(2.01*f_Ox_GB))
# The distance moved by the scale-metal interface since the last update of the scale grid.
# When delta_length reaches dl it's time to redefine the next gridpoint as scale. 
    delta_length_1 = 0.
    delta_length_2 = 0.
# Counter for iterations of the time loop between the single node transfers from gas to scale.
    t_count_1 = 0
# Counter for iterations of the time loop between the single node transfers from metal to scale.
    t_count_2 = 0
# Initial velocity of gas-scale interface ( < 0 )
    v_1 = v_1_0
# Initial velocity of scale-metal interface ( > 0 )
    v_2 = v_2_0
    if v_2*dt > dy:
         print('Error, length_increment will exceed grid spacing! v_1 =',v_1)
    length_increment_1 = v_1*dt
    length_increment_2 = v_2*dt
# t_count_1 counts the number of timesteps executed before extension of scale grid into gas
    t_count_1 = 0
# t_count_2 counts the number of timesteps executed before extension of scale grid into metal
    t_count_2 = 0
# Counts how many times a node is added to the scale
    n_nodes_added = 0
# Counts how many times the plotting function 'plot_concentrations' has been called
    plot_count = 0
    nfigx,nfigy = graph_layout
    print('Number of figures allowed per page = ',nfigx*nfigy)
# Make a list of timesteps at which a plot will be made
    time_stride = timesteps//(number_plots-2)
    list_plots = list(range(0,timesteps,time_stride))
    list_plots.append(timesteps)
    print('list_plots =',list_plots)
# Create the first figure for plots of the concentration through the scale
    new_page = plt.figure(figsize = (6,8))
    t=0
    new_page,plot_count = plot_concentrations(new_page,plot_count,'1.{}'.format(n_nodes_added), 3, \
                    graph_layout,number_plots, title='Fig.{}, t = {} '.format(plot_count+1,t), c_1=c_GB, c_2=c_bulk, \
                    label_1='GB', \
                    label_2='bulk')
    """
    The main time loop starts here
    
    """
    for t in range(1,timesteps+1):
        if delta_length_2 < dy-length_increment_2:
            for ny in range(ngridy_gas_0+1,ngridy_gas+ngridy_scale):
                c_GB[ny] = ( c_GB[ny] + (dd_GB[ny]- f_Ox_GB)*dt*((c_GB[ny+1]-c_GB[ny-1])/(2*dy))  + 
                             d_GB[ny]*dt*(c_GB[ny-1]+c_GB[ny+1]-2*c_GB[ny])/(dy*dy))
                temp = bleed*(c_GB[ny]-c_bulk[ny])
                c_GB[ny] = c_GB[ny]-temp*dt
                c_bulk[ny] = c_bulk[ny] + temp*nu*dt
# Apply zero gradient boundary condition at the scale metal interface
            c_GB[ngridy_gas+ngridy_scale-1] = c_GB[ngridy_gas+ngridy_scale-2]
# Increment the notional thickness of the scale
# This is increment of scale thickness due to the velocity of the scale-metal interface
            length_increment_2 = v_2*dt
            length = length+length_increment_2
            delta_length_2 = delta_length_2+length_increment_2
            v_2 = v_2_0/length  
            t_count_2 = t_count_2 + 1
# loop back unless it's time to grab another node of scale from the metal  
        else:
            print('\n* Seems like time to add a scale node from metal *\n')
            print('\nt_count_2 =',t_count_2 )
#           print('\nc_GB :\n',c_GB,'\n')
            print('delta_length_2 =', delta_length_2,  't =', t , 't_count_2 =', t_count_2,'\n')
            delta_length_2 = 0.0
            ngridy_gas, ngridy_scale, ngridy_metal, c_GB, vac_Ox, d_GB, dd_GB =  \
                update_ranges(1, ngridy_gas, ngridy_scale, ngridy_metal, c_GB, vac_Ox, d_GB, dd_GB)
            n_nodes_added = n_nodes_added + 1
            print('\nExpected number of timesteps before incrementing scale grid')
            print(' based on initial increment only = ',int(dy/(v_2*dt)), ' Actual number = ',t_count_2,'\n')
            print(' Number of nodes added to scale :',n_nodes_added)
            t_count_2=0
        if delta_length_1 < dy-length_increment_1:
            length_increment_1 = - v_1*dt
            length = length+length_increment_1
            delta_length_1 = delta_length_1+length_increment_1
            v_1 = v_1_0/length  
            t_count_1 = t_count_1 + 1
        else:
            print('\n* Seems like time to add a scale node from gas *\n')
            print('\nt_count_1 =',t_count_1 )
#           print('\nc_GB :\n',c_GB,'\n')
            print('delta_length_1 =', delta_length_1,  't =', t , 't_count_1 =', t_count_1,'\n')
            delta_length_1 = 0.0
            ngridy_gas, ngridy_scale, ngridy_metal, c_GB, vac_Ox, d_GB, dd_GB =  \
                update_ranges(-1, ngridy_gas, ngridy_scale, ngridy_metal, c_GB, vac_Ox, d_GB, dd_GB)
            n_nodes_added = n_nodes_added + 1
#           print('c_GB[ngridy_gas-1]',c_GB[ngridy_gas-1])
#           print('c_GB[ngridy_gas]',c_GB[ngridy_gas]) 
#           print('c_GB[ngridy_gas+ngridy_scale]',c_GB[ngridy_gas+ngridy_scale])
#           print('c_GB[ngridy_gas+ngridy_scale-1]',c_GB[ngridy_gas+ngridy_scale-1])
#           print('c_GB[ngridy_gas+ngridy_scale-2]',c_GB[ngridy_gas+ngridy_scale-2])
#           print('\nlength = ',length)
            print('\nExpected number of timesteps before adding node to scale')
            print(' based on initial increment only = ',-int(dy/(v_1*dt)), ' Actual number = ',t_count_1,'\n')
            print(' Number of nodes added to scale :',n_nodes_added)
            t_count_1=0
        if t == list_plots[plot_count]:
            new_page,plot_count = plot_concentrations(new_page,plot_count,'1.{}'.format(n_nodes_added), 3, \
                    graph_layout, number_plots,title='Fig.{}, t = {} '.format(plot_count+1,t), c_1=c_GB, c_2=c_bulk, \
                    label_1='GB', \
                    label_2='bulk')
    print('\n*** LEAVING tracer_diffusion ***\n')
    return plot_count, t, t_count_1, t_count_2, ngridy_gas, ngridy_scale, ngridy_metal, c_GB, c_bulk

#%%

"""
MAIN CODE FOR EXECUTION

"""
new_page = plt.figure()
plot_count, t,  t_count_1, t_count_2, ngridy_gas, ngridy_scale, ngridy_metal, c_GB, c_bulk = \
tracer_diffusion(50000, 1.9, ngridy_gas_0, ngridy_scale_0, ngridy_metal_0, vac_Ox_0, c_GB_0, d_GB_0, dd_GB_0,\
                 length_0, (3,2), bleed=0.0,number_plots=20)
os.chdir(working_directory)

print('c_GB =\n',c_GB)

