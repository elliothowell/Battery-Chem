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
from sys import exit
import numpy as np
import matplotlib.pyplot as plt
import math

###############################################################################


###############################################################################
"""
Assignment Reqs:
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

Using a computer, carry out simulations of cyclic voltammetry for a 
quasireversible system. Let l = 50 and D_M = 0.45. Take a = 0.5 and let the 
diffusion coefficients of the oxidized and reduced forms be equal. Cast your 
dimensionless intrinsic rate parameter in terms of the function ф defined in 
(6.5.5), and carry out calculations for ф = 20, 1, and 0.1. Compare the peak 
splittings in your simulated voltammograms with the values in Table 6.5.2.

"""
###############################################################################
# defining a function for the CV so I just need to plug in varying psi, then 
# I can also vary the other model parameters if needed (which is often as using
# l = 50 can be brutal)
def CV_Sim(psi, timeSteps = 50, modelSize = 50, alpha = 0.5, 
           Ei = 1.0, Ef = -1.0, c_Ai = 1.0, c_Bi = 0.0):    

    ###########################################################################
    # Constants for simulation (making potential dimensionless)
    F = 96485.332 # Faradays constant (units = C/mol)
    R = 8.314 # Gas Constant (units = J/mol*K)
    T = 298.15 # Temperature of the system (units = K)
    f = F / ( R * T)    
    
    # Setting all of the parameters for the system
    modelSize = modelSize #Size/number of boxes of the model
    if modelSize < 2:
        print("Model size must be greater than 1!")
        exit()
    '''
    Super interesting action here, when modelSize is decreased it becomes similar
    to a microelectrode, no issues with diffusion are present!
    '''
    timeSteps = timeSteps # AKA l # Number of time steps

    # Model diffusion coefficient for A and B (can also define parameters of this)
    #D_MA = Diff*tInterval/(xInterval**0.5)
    D_MA = 0.45 # where D_M = D*dt/dx^2
    D_MB = 0.45
    
    # Initial concentrations for species A and species B (units of molar) and
    # other inherent parameters of redox system
    # Typical start will just be A = 1M and B = 0
    C_Ai = c_Ai
    if C_Ai == 0.0: #setting fractional concentrations of A and in the event A == 0
        f_Ai = 0.0  #then just set A to 0
    else: 
        f_Ai = C_Ai/C_Ai # setting fractional concentrations of A

    C_Bi = c_Bi
    if C_Bi == 0.0: #setting fractional concentrations of B and in the event B == 0
        f_Bi = 0.0  #then just set f_Bi to 0

    else:
        f_Bi = C_Bi / (C_Ai + C_Bi)
        f_Ai = 1 - f_Bi
        
    # Potential Parameters (typically start at 1 -> -1 V but can vary)
    Einorm = Ei # starting potential
    Efnorm = Ef # ending potential
    
    # assigning alpha value (standard is 0.5)
    alpha = alpha # transfer coeff

    # Setting arrays for A and B with initial values
    # using 2D array (time = columns, rows = x space)
    f_A_array = np.full(shape = (modelSize, timeSteps+1), fill_value = f_Ai) # 
    f_B_array = np.full(shape = (modelSize, timeSteps+1), fill_value = f_Bi)
    f_A_hold = 0.0
    f_B_hold = 0.0    
    
    # Arrays which are looking at time/x-space from electrode to end of meshing
    Chi_array = np.full(shape = modelSize, fill_value = 0.0) # array for chi

    # Creating arrays for data storage of current, potential, concentrations
    # these should all correspond to specific time step
    # for this assignment only Z, Z_Cott, t/tk are really used
    Enorm = np.full(shape = timeSteps, fill_value = 0.0) #checking what potentials we go thru
    measured_Z = np.full(shape = timeSteps, fill_value = 0.0)#dimensionless current
    t_over_tk = np.arange(start = 0.0, stop = 1, step = 1/timeSteps) # dimensionless time
    iteration = np.arange(start = 1, stop = timeSteps + 1, step = 1) # list of iterations, not explicitely used but can be useful
    ###########################################################################
    
    for t in range(0, timeSteps):

        ### Just for potential calcs ###
        #######################################################################
        # this is basically the reversal step for potential change
        if t == timeSteps // 2:
            Ehold = Einorm
            Einorm = Efnorm
            Efnorm = Ehold

        # calculating what potential we at, going from Ei->Ef->ei
        if t < timeSteps // 2:
            Enorm[t] = Einorm + (Efnorm - Einorm) * (t/((timeSteps-1)/2))
            
        else:
            Enorm[t] = Einorm + (Efnorm - Einorm) * (
                (t-(timeSteps-1)/2)/((timeSteps-1)/2))
        #######################################################################
    
    # This is just making sure the math will math, no imaginary numbers here
    psiConst = (psi * math.sqrt(math.pi * abs(Einorm - Efnorm)))
    
    for t in range(1, timeSteps+1):
        
        # This will be dtermining dimensionless rate constants for A and B
        dimRatA = psiConst * math.exp(-alpha * f * Enorm[t-1]) 
        dimRatB = psiConst * math.exp((1-alpha) * f * Enorm[t-1]) 
        
        # calculation of the current given calculated dimensional rate constants
        # and previous iterations concentration of A and B respectively
        curr = -(((dimRatA * f_A_array[1, t-1]) 
                  - (dimRatB * f_B_array[1, t-1])) / (1 + dimRatA + dimRatB))
        
        # Calculation of new concentration of species A and B due to current
        f_A_hold = f_A_array[1, t-1] + curr
        f_B_hold = f_B_array[1, t-1] - curr
        
        # ensuring that we don't go into negative spiral if conc is negative 
        # this is highly unlikely but a keeps code stable
        if f_A_hold < 0:
            f_A_hold = 0
        if f_B_hold < 0:
            f_B_hold = 0
        
        # assignment of calculated frac concs to array
        f_A_array[0, t] = f_A_hold
        f_B_array[0, t] = f_B_hold
        
        # assigning current response to dimensionless current array
        measured_Z[t-1] = curr
            
        ####
        # SOLELY FOR DIFFUSION
        # start cycling through boxes for each speciic time interval
        ####
        for j in range(1, modelSize-1):
            
            # Calculation of Chi only needs 1 iteration as doing it more would
            # just take up computational power
            if t == 1:
                Chi_array[j] = (j)/((D_MA*timeSteps)**0.5)
                
            # Diffusion beyond the first box
            # else:
            f_A_array[j,t] = (f_A_array[j,t-1] 
                              + D_MA*(f_A_array[j+1,t-1] 
                                      - (2*f_A_array[j,t-1]) 
                                      + f_A_array[j-1,t-1]))
            
            f_B_array[j,t] = (f_B_array[j,t-1] 
                              + D_MB*(f_B_array[j+1,t-1] 
                                      - (2*f_B_array[j,t-1]) 
                                      + f_B_array[j-1,t-1]))
        
        # Calculation of peak splitting based off max and min current response
        # Psi == 1 is different as it was picking peak too far away, brute manual override...
        if psi == 1:
            peakSplit = int(1000 * (Enorm[np.where(measured_Z == max(measured_Z))[0][0]] 
                                    - Enorm[np.where(measured_Z == min(measured_Z))[0][0] - 1]))
        else:   
            peakSplit = int(1000 * (Enorm[np.where(measured_Z == max(measured_Z))[0][0]] 
                                    - Enorm[np.where(measured_Z == min(measured_Z))[0][0]]))
            
    # want to return in order of E, Z, t_tk, f_A, f_B
    return Enorm, measured_Z, t_over_tk, Chi_array, f_A_array[:, :-1], f_B_array[:, :-1], peakSplit
###############################################################################
"""
When calling output of function potisions are as follows:
    [0]: Potential range
    [1]: Measured dimensionless current
    [2]: Dimensionless time array
    [3]: Dimensionless distance from electrode array
    [4]: Fractional concentration of species A with size [distance, time]
    [5]: Fractional concentration of species B with size [distance, time]
    [6]: Peak splitting based off of min and max current values, when resolution is low this is unreliable

CV_Sim(psi, timeSteps, modelSize, alpha, Ei, Ef, c_Ai, c_Bi)
"""  
# Setting adjustable variables so same for each psi adjustment
ei = 0.2 # starting potential
ef = -0.2 # ending potential
ms = 50 # model size / meshing (smaller = microelectrode actions)

# Creating arrays for when ps = 20, 1, and 0.2
psi20 = CV_Sim(psi = 20, modelSize = ms, timeSteps = 50, Ei = 0.15, Ef = -0.2)
print(psi20[-1])
psi1 = CV_Sim(psi = 1, modelSize = ms, Ei= 0.15, Ef= -0.2)
print(psi1[-1])
psi0_1 = CV_Sim(psi = 0.1, modelSize = ms, Ei= 0.5, Ef= -0.5) 
print(psi0_1[-1])
microElec = CV_Sim(psi = 1, modelSize = 3, timeSteps = 500, Ei= 0.5, Ef= -0.5) 
###############################################################################
"""
Block for creating plots of the data
"""
###############################################################################

## Plot of CV for psi 20 where max and min current responses are circled
psi20MaxCurr_index = np.where(psi20[1] == max(psi20[1]))
psi20MinCurr_index = np.where(psi20[1] == min(psi20[1]))
plt.figure()
plt.title("CV for when Psi = 20 for the range E = " + str(max(psi20[0])) + " to " + str(round(min(psi20[0]), 1)))
plt.plot(psi20[0], psi20[1] , linewidth = 2, markersize = 5)
plt.scatter(psi20[0][psi20MaxCurr_index], psi20[1][psi20MaxCurr_index], 
            facecolors = 'none', edgecolors= 'r', s = 100, label = "Max Current")
plt.scatter(psi20[0][psi20MinCurr_index], psi20[1][psi20MinCurr_index], 
            facecolors = 'none', edgecolors= 'b', s = 100, label = "Minimum Current")
plt.xlabel("Potential (V)")
plt.ylabel("Dimensionless Current (Z)")
plt.legend()
plt.show()

## Plot of CV for psi 1 where max and min current responses are circled
psi1MaxCurr_index = np.where(psi1[1] == max(psi1[1]))
psi1MinCurr_index = np.where(psi1[1] == min(psi1[1]))[0][0] - 1
plt.figure()
plt.title("CV for when Psi = 1 for the range E = " + str(max(psi1[0])) + " to " + str(round(min(psi1[0]), 1)))
plt.plot(psi1[0], psi1[1] , linewidth = 2, markersize = 5)
plt.scatter(psi1[0][psi1MaxCurr_index], psi1[1][psi1MaxCurr_index], 
            facecolors = 'none', edgecolors= 'r', s = 100, label = "Max Current")
plt.scatter(psi1[0][psi1MinCurr_index], psi1[1][psi1MinCurr_index], 
            facecolors = 'none', edgecolors= 'b', s = 100, label = "Minimum Current")
plt.xlabel("Potential (V)")
plt.ylabel("Dimensionless Current (Z)")
plt.legend()
plt.show()

## Plot of CV for psi 0.11 where max and min current responses are circled
psi0_1MaxCurr_index = np.where(psi0_1[1] == max(psi0_1[1]))
psi0_1MinCurr_index = np.where(psi0_1[1] == min(psi0_1[1]))
plt.figure()
plt.title("CV for when Psi = 0.1 for the range E = " + str(max(psi0_1[0])) + " to " + str(round(min(psi0_1[0]), 1)))
plt.plot(psi0_1[0], psi0_1[1] , linewidth = 2, markersize = 5)
plt.scatter(psi0_1[0][psi0_1MaxCurr_index], psi0_1[1][psi0_1MaxCurr_index], 
            facecolors = 'none', edgecolors= 'r', s = 100, label = "Max Current")
plt.scatter(psi0_1[0][psi0_1MinCurr_index], psi0_1[1][psi0_1MinCurr_index], 
            facecolors = 'none', edgecolors= 'b', s = 100, label = "Minimum Current")
plt.xlabel("Potential (V)")
plt.ylabel("Dimensionless Current (Z)")
plt.legend()
plt.show()

plt.figure()
plt.title("CV for Microelectrode behaviour")
plt.plot(microElec[0], microElec[1] , linewidth = 2, markersize = 5)

# print("peak splitting: ", int(1000 * (Enorm[np.where(measured_Z == max(measured_Z))[0][0]] - Enorm[np.where(measured_Z == min(measured_Z))[0][0]])), "mV")
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

###############################################################################