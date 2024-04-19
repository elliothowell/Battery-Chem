# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 11:05:10 2024

This code is the beginning step for numerical simulation, this is following 
steps laid out by Allen J. Bard and Larry R. Faulkner, accompanied work 
produced by Lisa I. Stephens and Janine Mauzeroll. 

The first step here is to produce code which solves a simulation of a cyclic 
voltammorgram laid out in appendix B.5 within Electrochemical Methods 
Fundamentals and Applications

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

# Potential Parameters!!!!
Einorm = 1 #dimensionless starting potential
Efnorm = -1 #dimensionless ending potential

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
Enorm = np.full(shape = timeSteps, fill_value = 0.0) #checking what potentials we go thru
measured_Z = np.full(shape = timeSteps, fill_value = 0.0)#dimensionless current
re_Cott = np.full(shape = timeSteps, fill_value = 0.0) #dimensionless current
Z_Cott = np.full(shape = timeSteps, fill_value = 0.0) #dimensionless current
t_over_tk = np.arange(start = 0.0, stop = 1, step = 1/timeSteps)
iteration = np.arange(start = 1, stop = timeSteps, step = 1)


testZ = np.full(shape = timeSteps, fill_value = 0.0)

###############################################################################
# Assignment Requirements
"""
Reqs:
    - l = 50 (timeSteps)
    - D_M = 0.45
    - a = 0.5
    - Quasireversible!
    - diff coeff of oxidized and reduced forms equal
    - cast dimensionless intrinsic rate param in terms of function psi 
    which is defined in (6.5.5)
    - carry out calculations for psi = 20, 1, and 0.1
    - compare peak splitting in simulated voltammograms with values in table
    6.5.2:
        psi         dE_p (E_pa - E_pc) 
                    (mV)
        ------------------------------
        20          61
        7           63
        6           64
        5           65
        4           66
        3           68
        2           72
        1           84
        0.75        92
        0.50        105
        0.35        121
        0.25        141
        0.10        212
        ------------------------------
    
    upperLambda = k^o / (D_O^(1-a) * D_R^(a) * F * v)^1/2
    D_O=D_R=D   = k^o / (D * F * v)^1/2
    
    (6.5.5) psi = upperLambda * pi^-1/2
                = ((D_O/D_R)^(a/2) * k^o) / (pi * D_O * f * v)^1/2

"""
'''
Some discussion on potential sweeps in modelling (Bard B.4.3):
    - since this is a CV we want to have a cyclic process, so going from Einorm 
    to Efnorm and then back to Einorm?
    1V, to -1V, back to 1V
        - first cycle is sweep from 1V to -1V, then using all of the data sweep
        back from -1V to 1V
    - using this we shouldn't have to change t/tk correct??
    
    

'''


# Starting the simulation

### THIS SHOULD BE MOVED/ALTERED ###
psi = 20
psiConst = (psi * math.sqrt(math.pi * (Einorm - Efnorm)))

for t in range(0, timeSteps):
    
    # this is basically the reversal step for potential change
    if t == timeSteps/2:
        Ehold = Einorm
        Einorm = Efnorm
        Efnorm = Ehold

    # calculating what potential we at, going from Ei->Ef->ei
    if t < timeSteps/2:
        Enorm[t] = Einorm + (Efnorm - Einorm) * (t/((timeSteps-1)/2))
        
    else:
        Enorm[t] = Einorm + (Efnorm - Einorm) * ((t-(timeSteps-1)/2)/((timeSteps-1)/2))
        #print((t-(timeSteps-1)/2)/((timeSteps-1)/2))
        
    # current meas.
    testZ[t] = psiConst * (math.exp(-alpha * Enorm[t]) - math.exp((1-alpha) * Enorm[t]))
        
    print(math.exp((1-alpha) * Enorm[t]))
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
        

# Setting up array for Z/Z_Cott once simulation has finished
#R_array = measured_Z/Z_Cott

###############################################################################
# Plotting
# Plot time in terms of dt/2 so x axis will be time array but points will be
# time point - tInterval/2

# Setting warning in the event inputted values don't allow for the t/tk <=0.2
# if timeSteps > modelSize:
#     plotIndex = modelSize
#     print("Warning: the amount of time steps exceeds the size of the model ",
#           "and therefore will be the max plot amount.")
# else : 
#     plotIndex = np.where(t_over_tk <= 0.2)[-1][-1]
# print(plotIndex)

## Plot of Z/Z_cot versus t/t_k
plt.figure()
plt.title("Plot of Z against normalized time (t/$t_{k}$)")
plt.plot(t_over_tk, testZ, linewidth = 2, markersize = 5, 
          marker = 'o')
plt.xlabel('t/$t_{k}$')
plt.ylabel('Z')
plt.show()

# ## Plot of percent error between Z and Z_cot versus t/t_k
# plt.figure()
# plt.title("Plot of percent error between Z and $Z_{Cot}$ against " +
#           "normalized time (t/$t_{k}$) \n" +
#           " for the first ten iterations")
# plt.plot(t_over_tk[:10], 100*(measured_Z[:10]-Z_Cott[:10])/Z_Cott[:10], 
#          linewidth = 2, markersize = 5, 
#           marker = 'o')
# plt.xlabel('t/$t_{k}$')
# plt.ylabel('Percent Difference ((Z-$Z_{Cot}$)/$Z_{Cot}$)')
# plt.show()

# ## Plot of normalized concentrations (sim/analytical) vs Chi for t/t_k=0.2
# plt.figure()
# plt.title("Plot of $f_{A}$ and $f_{B}$ versus Chi ($\\chi$) for" + 
#           " t/$t_{k}$ = 0.2")
# plt.plot(Chi_array[:12] , 
#           f_A_array[:12, plotIndex], 
#           linewidth = 2, label = 'Simulated $f_{A}$')
# plt.scatter(Chi_array[:12] , 
#           f_A_erf[:12, plotIndex], label = 'Analytical $f_{A}$')
# plt.plot(Chi_array[:12] , 
#           f_B_array[:12, plotIndex] , 
#           linewidth = 2, label = 'Simulated $f_{B}$')
# plt.scatter(Chi_array[:12] , 
#           f_B_erf[:12, plotIndex], label = 'Analytical $f_{B}$')
# plt.xlabel('Chi ($\\chi$)')
# plt.ylabel('Normalized concentration')
# plt.legend()
# plt.show()

# ## Plot of percent errors for f_A and f_B (sim/analytical) vs Chi for t/t_k=0.2
# plt.figure()
# plt.title("Plot of $f_{A}$ and $f_{B}$ versus Chi ($\\chi$) for" + 
#           " t/$t_{k}$ = 0.2")
# plt.plot(Chi_array[:12] , 
#           100*(f_A_array[:12, plotIndex]-
#                f_A_erf[:12, plotIndex])/f_A_erf[:12, plotIndex], 
#           linewidth = 2, label = '($f_{A sim.}$-$f_{A anal.}$)/$f_{A anal.}$')
# plt.plot(Chi_array[:12] , 
#           100*(f_B_array[:12, plotIndex]-
#                f_B_erf[:12, plotIndex])/f_B_erf[:12, plotIndex], 
#           linewidth = 2, label = '($f_{B sim.}$-$f_{B anal.}$)/$f_{B anal.}$')
# plt.xlabel('Chi ($\\chi$)')
# plt.ylabel('Percent Difference (($f_{sim.}$-$f_{anal.}$)/$f_{anal.}$)')
# plt.legend()
# plt.show()

###############################################################################