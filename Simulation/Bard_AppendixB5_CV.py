# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 11:05:10 2024

This code is the beginning step for numerical simulation, this is following 
steps laid out by Allen J. Bard and Larry R. Faulkner, accompanied work 
produced by Lisa I. Stephens and Janine Mauzeroll. 

The first step here is to produce code which solves a Cottrell experiment 
laid out in appendix B.2 within Electrochemical Methods Fundamentals and 
Applications

@author: Elliot
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import math

###############################################################################
# Constants for simulation 
## Not required here
F = 96485.332 # Faradays constant (units = C/mol)
R = 8.314 # Gas Constant (units = J/mol*K)
T = 298.15 # Temperature of the system (units = K)

###############################################################################
# Input parameters of the system (initial params)

modelSize = 50 #Size/number of boxes of the model
timeSteps = 50 # AKA l # Number of time steps

# Model diffusion coefficient for A and B (can also define parameters of this)
#D_MA = Diff*tInterval/(xInterval**0.5)
D_MA = 0.40 # where D_M = D*dt/dx^2
D_MB = 0.40

# Initial concentrations for species A and species B (units of molar) and
# other inherent parameters of redox system
# Typical start will just be A = 1M and B = 0
C_Ai = 1.0e-6
if C_Ai == 0.0: #setting fractional concentrations of A and in the event A == 0
    f_Ai = 0.0  #then just set A to 0
else: 
    f_Ai = C_Ai/C_Ai # setting fractional concentrations of A

C_Bi = 0.0
if C_Bi == 0.0: #setting fractional concentrations of B and in the event B == 0
    f_Bi = 0.0  #then just set f_Bi to 0
else:
    f_Bi = C_Bi/C_Bi
    
# These are not needed for this assignment but will be useful later
Diff = 1e-6 # diffusion coefficient for species A/B, units = cm^2/s
n = 1 # number of electrons transferred in redox
alpha = 0.5 # transfer coeff

# Electrode parameters
# Again not needed but saving for future iterations
Area = 1 # electrode area, units = cm^2 (default = 0.1 cm^2)
k_0 = 0.1 # rate constant, units = m/s (default = 0.1 m/s)

# One dimensional array for time
# Not completely necessary but useful for visualizing time span
tSize = 10 # AKA t_k # Length of time, units = seconds
tInterval = tSize / timeSteps # delta t # time interval
sampleTime = np.arange(start = 0.0, stop = tSize, step = tInterval)+tInterval/2

# One dimensional array for the physical x-axis
# yet again, not necessary but helps with visualization
JMax = 4.2*(timeSteps**0.5) # useful for optimizing processing power
xInterval = ((Diff * tSize) / (D_MA * timeSteps))**0.5 # delta x 
xSize =  modelSize * xInterval # units = cm
x = np.arange(start = 0.0, stop = xSize, step = xInterval) 

# Setting arrays for A and B with initial values
# using 2D array (time = columns, rows = x space)
f_A_array = np.full(shape = (modelSize, timeSteps), fill_value = f_Ai)
f_B_array = np.full(shape = (modelSize, timeSteps), fill_value = f_Bi)
f_A_array_holder = np.full(shape = (modelSize, timeSteps), fill_value = f_Ai)
f_B_array_holder = np.full(shape = (modelSize, timeSteps), fill_value = f_Bi)

# Arrays which are looking at time/x-space from electrode to end of meshing
f_A_erf = np.full(shape = (modelSize, timeSteps), fill_value = 1.0)#analyticalA
f_B_erf = np.full(shape = (modelSize, timeSteps), fill_value = 0.0)#analyticalB
Chi_array = np.full(shape = modelSize, fill_value = 0.0) # array for chi

# Creating arrays for data storage of current, potential, concentrations
# these should all correspond to specific time step
# for this assignment only Z, Z_Cott, t/tk are really used
measured_I = np.full(shape = timeSteps, fill_value = 0.0) #current (Amps)
measured_Z = np.full(shape = timeSteps, fill_value = 0.0)#dimensionless current
re_Cott = np.full(shape = timeSteps, fill_value = 0.0) #dimensionless current
Z_Cott = np.full(shape = timeSteps, fill_value = 0.0) #dimensionless current
t_over_tk = np.arange(start = 0.0, stop = 1, step = 1/timeSteps)
iteration = np.arange(start = 1, stop = timeSteps, step = 1)

###############################################################################
# Assignment Requirements
"""

C_A(j,0) = C_Ai -> k = 0 iteration all boxes should contain a 1
C_A(1,1) = 0 -> When k > 0 box 1 should contain 0
C_A(j_max,0) = C_Ai -> Final box in model needs to be fixed at 1

OK Actually reading through problem!
    - Work through first 10 iterations (time steps?)
    - l = 50 (timeSteps)
    - D_M = 0.4
    - Calc Z(k) for each iteration
    - Calc Z_Cott(k) for each iteration
    - Compare Z(k) and Z_Cott(k)
    - Calculate X (chi) values corresponding to first 12 boxes (model size)
    - Plot conc profiles f_A and f_B vs X (chi) for t/t_k = 0.2
        - THIS IS JUST ASKING THAT you plot up to t = 0.2*t_k
    - Derive functions describing f_A and f_B vs X (chi) and t/t_k 
      from (5.2.13)
        - In this case t/t_k = k/l
    - Draw analytical curves on graphs of concentration profiles
    - Comment on agreement between model and known solution

"""

# Starting the simulation

for t in range(0, timeSteps):
    
    # for each iteration just go ahead and calculate Z_cot
    Z_Cott[t] = 1/((np.pi*((t+0.5)/timeSteps))**0.5)
    
    # instantaneous start only happens for t = 0
    if t == 0:
        # Instant conversion of bulk at electrode surface from A -> B
        f_B_array[0, t] = f_A_array[0, t]
        f_B_array_holder[0, t] = f_A_array_holder[0, t]
        f_A_array[0,t] = 0
        f_A_array_holder[0,t] = 0

        # not req. but sanity check with provided resources
        measured_I[t] = n*F*Area*(Diff**0.5)*C_Ai*(timeSteps**0.5)/(
                        (tSize**0.5)*(D_MA**0.5))
        
        # calculate dimensionless current initial
        measured_Z[t] = ((timeSteps/D_MA)**0.5)

        # not req. but sanity check with provided resources
        re_Cott[t] = n*F*Area*C_Ai*(Diff**0.5)/(
                     (np.pi**0.5)*(sampleTime[0]**0.5))
        
        #analytical results calculaton, can't divide by 0 so manual punch (t=0)
        f_A_erf[0,t] = 0 # f_Ai * math.erf((0)/(2*((D_MA*t)**0.5)))
        f_B_erf[0,t] = 1 #- (f_Ai * math.erf((0)/(2*((D_MA*t)**0.5))))
        
    else:
        
        # making sure array lines up with previous iteration to perform calcs
        f_A_array[:,t] = f_A_array[:,t-1]
        f_B_array[:,t] = f_B_array[:,t-1]
        
        # start cycling through boxes for each speciic time interval
        for k in range(modelSize):
            
            # calculate analytical results
            f_A_erf[k,t] = f_Ai * math.erf((k)/(2*((D_MA*t)**0.5)))
            f_B_erf[k,t] = 1 - (f_Ai * math.erf((k)/(2*((D_MA*t)**0.5))))
            
            # Calculation of Chi only needs 1 iteration as doing it more would
            # just take up computational power
            if t == 1:
                Chi_array[k] = (k)/((D_MA*timeSteps)**0.5)
            
            # Diffusion into the first box
            if k == 0:
                f_A_array_holder[k,t] = f_A_array[k,t] + D_MA*(
                    f_A_array[k+1,t] - f_A_array[k,t])
                f_B_array_holder[k,t] = f_B_array[k,t] + D_MB*(
                    f_B_array[k+1,t] - f_B_array[k,t])
            
            # End boundary condition of C_A(j_max, 0) = C_Ai
            elif k == modelSize-1:
                f_A_array_holder[k] = f_Ai
            
            # Diffusion beyond the first box
            else:
                f_A_array_holder[k,t] = f_A_array[k,t] + D_MA*((
                    f_A_array[k+1,t] - 2*f_A_array[k,t]) + f_A_array[k-1,t])
                f_B_array_holder[k,t] = f_B_array[k,t] + D_MB*((
                    f_B_array[k+1,t] - 2*f_B_array[k,t]) + f_B_array[k-1,t])
        
        # Replace values with newly calculated values of 
        # fractional concentration
        f_A_array[:,t] = f_A_array_holder[:,t]
        f_B_array[:,t] = f_B_array_holder[:,t]
        
        # calculate dimensionless current due to changes
        measured_Z[t] = ((timeSteps/D_MA)**0.5) * f_A_array_holder[0,t]
        
        # Again, immediate electrolysis of species A to B at surface
        f_A_array[0,t] = 0
        f_B_array[0,t] = f_Ai
        
        # Calculate the measured current due to change
        measured_I[t] = ((n*F*Area*(Diff**0.5)*C_Ai)*(
            (f_A_array[1,t-1])*((D_MA*timeSteps)**0.5) / (tSize**0.5)))
        
        re_Cott[t] = n*F*Area*C_Ai*(Diff**0.5)/(
            (np.pi**0.5)*(sampleTime[t]**0.5))

# Setting up array for Z/Z_Cott once simulation has finished
R_array = measured_Z/Z_Cott

###############################################################################
# Plotting
# Plot time in terms of dt/2 so x axis will be time array but points will be
# time point - tInterval/2

# Setting warning in the event inputted values don't allow for the t/tk <=0.2
if timeSteps > modelSize:
    plotIndex = modelSize
    print("Warning: the amount of time steps exceeds the size of the model ",
          "and therefore will be the max plot amount.")
else : 
    plotIndex = np.where(t_over_tk <= 0.2)[-1][-1]
print(plotIndex)

## Plot of Z/Z_cot versus t/t_k
plt.figure()
plt.title("Plot of Z/$Z_{Cot}$ against normalized time (t/$t_{k}$) \n" +
          " for the first ten iterations")
plt.plot(t_over_tk[:10], R_array[:10], linewidth = 2, markersize = 5, 
          marker = 'o')
plt.xlabel('t/$t_{k}$')
plt.ylabel('R (Z/$Z_{Cot}$)')
plt.show()

## Plot of percent error between Z and Z_cot versus t/t_k
plt.figure()
plt.title("Plot of percent error between Z and $Z_{Cot}$ against " +
          "normalized time (t/$t_{k}$) \n" +
          " for the first ten iterations")
plt.plot(t_over_tk[:10], 100*(measured_Z[:10]-Z_Cott[:10])/Z_Cott[:10], 
         linewidth = 2, markersize = 5, 
          marker = 'o')
plt.xlabel('t/$t_{k}$')
plt.ylabel('Percent Difference ((Z-$Z_{Cot}$)/$Z_{Cot}$)')
plt.show()

## Plot of normalized concentrations (sim/analytical) vs Chi for t/t_k=0.2
plt.figure()
plt.title("Plot of $f_{A}$ and $f_{B}$ versus Chi ($\\chi$) for" + 
          " t/$t_{k}$ = 0.2")
plt.plot(Chi_array[:12] , 
          f_A_array[:12, plotIndex], 
          linewidth = 2, label = 'Simulated $f_{A}$')
plt.scatter(Chi_array[:12] , 
          f_A_erf[:12, plotIndex], label = 'Analytical $f_{A}$')
plt.plot(Chi_array[:12] , 
          f_B_array[:12, plotIndex] , 
          linewidth = 2, label = 'Simulated $f_{B}$')
plt.scatter(Chi_array[:12] , 
          f_B_erf[:12, plotIndex], label = 'Analytical $f_{B}$')
plt.xlabel('Chi ($\\chi$)')
plt.ylabel('Normalized concentration')
plt.legend()
plt.show()

## Plot of percent errors for f_A and f_B (sim/analytical) vs Chi for t/t_k=0.2
plt.figure()
plt.title("Plot of $f_{A}$ and $f_{B}$ versus Chi ($\\chi$) for" + 
          " t/$t_{k}$ = 0.2")
plt.plot(Chi_array[:12] , 
          100*(f_A_array[:12, plotIndex]-
               f_A_erf[:12, plotIndex])/f_A_erf[:12, plotIndex], 
          linewidth = 2, label = '($f_{A sim.}$-$f_{A anal.}$)/$f_{A anal.}$')
plt.plot(Chi_array[:12] , 
          100*(f_B_array[:12, plotIndex]-
               f_B_erf[:12, plotIndex])/f_B_erf[:12, plotIndex], 
          linewidth = 2, label = '($f_{B sim.}$-$f_{B anal.}$)/$f_{B anal.}$')
plt.xlabel('Chi ($\\chi$)')
plt.ylabel('Percent Difference (($f_{sim.}$-$f_{anal.}$)/$f_{anal.}$)')
plt.legend()
plt.show()

###############################################################################