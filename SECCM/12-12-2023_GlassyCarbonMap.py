# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 16:31:53 2023

OK SO, asc file is very difficult to deal with, just going to move along with .mat file

Core points: 
    Trace_x_x_x_1 == Voltage/Potential
    Trace_x_x_x_2 == Current 
    Trace_x_x_1,2,3_x = Rep of CV 

@author: Elliot
"""

import os
import scipy.io as sp
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io as sp
import scipy.signal as ss
import numpy as np

path = 'C:\\Users\\Elliot\\SynologyDrive\\Research - Elliot Howell\\McGill - AC Corrosion\\Scanny\\12-12-2023_GlassyCarbon'

file = '12-12-2023_GlassySample_Map.mat'

os.chdir(path)

data = sp.loadmat(file)

dataKeys = list(data.keys())

i = 20
j = 1
n = 1


while i < (len(dataKeys)-1):
    # print("i: ", i)
    # print("i data key: ", dataKeys[i])
    # print("i+1 data key: ", dataKeys[i+1])
    # print("i-1 data key: ", dataKeys[i-1])
    if (int(dataKeys[i-1][-1]) - int(dataKeys[i-2][-1]) == 0) and (int(dataKeys[i][-1]) - int(dataKeys[i-1][-1]) == -1):
        
        graphName = "Cyclic Voltammograms of Landing " + str(n)
        
        current1 = data[dataKeys[i]][:,1]
        voltage1 = data[dataKeys[i+1]][:,1]
        smoothCurrent1 = ss.savgol_filter(current1,20,3)
        
        current2 = data[dataKeys[i+2]][:,1]
        voltage2 = data[dataKeys[i+3]][:,1]
        smoothCurrent2 = ss.savgol_filter(current2,20,3)
        
        current3 = data[dataKeys[i+4]][:,1]
        voltage3 = data[dataKeys[i+5]][:,1]
        smoothCurrent3 = ss.savgol_filter(current3,20,3)
        
        # smoothCurrent = ss.savgol_filter(current3,20,3)
        
        plt.figure()
        
        plt.ylim([-0.6e-9, 0.6e-9])
        plt.title(graphName)
        plt.xlabel("Voltage (V vs. Ag/AgCl)")
        plt.ylabel("Current (A)")
        # fig, ax = plt.subplots()
        plt.plot(voltage1, smoothCurrent1, label = "Replicate 1", linewidth = 1, markersize = 0.5)
        plt.plot(voltage2, smoothCurrent2, label = "Replicate 2", linewidth = 1, markersize = 0.5)
        plt.plot(voltage3, smoothCurrent3, label = "Replicate 3", linewidth = 1, markersize = 0.5)
        plt.legend()
        #plt.savefig(graphName+'.tif', dpi = 300)
        plt.show()
        
        # print("n: ", n)
        n +=1
        # print("Trace: ", dataKeys[i])
        

    i+=1

# dataKeys
# Out[78]: 
# ['__header__',
#  '__version__',
#  '__globals__',
#  'Trace_1_1_1_2',
#  'Trace_1_2_1_2',
#  'Trace_1_3_1_1',
#  'Trace_1_3_1_2',
#  'Trace_1_4_1_1',
#  'Trace_1_4_1_2',
#  'Trace_1_4_2_1',
#  'Trace_1_4_2_2',
#  'Trace_1_4_3_1',
#  'Trace_1_4_3_2',
#  'Trace_1_5_1_1',
#  'Trace_1_5_1_2',
#  'Trace_1_5_2_1',
#  'Trace_1_5_2_2',
#  'Trace_1_5_3_1',
#  'Trace_1_5_3_2',
#  'Trace_1_6_1_2',
#  'Trace_1_7_1_1', --> ITEM 20!
#  'Trace_1_7_1_2',
#  'Trace_1_8_1_1',
#  'Trace_1_8_1_2',
#  'Trace_1_9_1_1',
#  'Trace_1_9_1_2',
#  'Trace_1_10_1_2',
#  'Trace_1_11_1_1',
#  'Trace_1_11_1_2',
#  'Trace_1_12_1_1',
#  'Trace_1_12_1_2',
#  'Trace_1_13_1_1',
#  'Trace_1_13_1_2',
#  'Trace_1_14_1_2',
#  'Trace_1_15_1_1',
#  'Trace_1_15_1_2',
#  'Trace_1_16_1_1',
#  'Trace_1_16_1_2',
#  'Trace_1_17_1_1',
#  'Trace_1_17_1_2',
#  'Trace_1_18_1_2',
#  'Trace_1_19_1_1',
#  'Trace_1_19_1_2',
#  'Trace_1_20_1_1',
#  'Trace_1_20_1_2',
#  'Trace_1_21_1_1',
#  'Trace_1_21_1_2',
#  'Trace_1_22_1_2',
#  'Trace_1_23_1_1',
#  'Trace_1_23_1_2',
#  'Trace_1_24_1_1',
#  'Trace_1_24_1_2',
#  'Trace_1_25_1_1',
#  'Trace_1_25_1_2',
#  'Trace_1_26_1_2',
#  'Trace_1_27_1_1',
#  'Trace_1_27_1_2',
#  'Trace_1_28_1_1',
#  'Trace_1_28_1_2',
#  'Trace_1_29_1_1',
#  'Trace_1_29_1_2',
#  'Trace_1_30_1_2',
#  'Trace_1_31_1_1',
#  'Trace_1_31_1_2',
#  'Trace_1_32_1_1',
#  'Trace_1_32_1_2',
#  'Trace_1_33_1_1',
#  'Trace_1_33_1_2',
#  'Trace_1_34_1_2',
#  'Trace_1_35_1_1',
#  'Trace_1_35_1_2',
#  'Trace_1_36_1_1',
#  'Trace_1_36_1_2',
#  'Trace_1_37_1_1',
#  'Trace_1_37_1_2',
#  'Trace_1_38_1_2',
#  'Trace_1_39_1_1',
#  'Trace_1_39_1_2',
#  'Trace_1_40_1_1',
#  'Trace_1_40_1_2',
#  'Trace_1_41_1_1',
#  'Trace_1_41_1_2',
#  'Trace_1_42_1_2',
#  'Trace_1_43_1_1',
#  'Trace_1_43_1_2',
#  'Trace_1_44_1_1',
#  'Trace_1_44_1_2',
#  'Trace_1_45_1_1',
#  'Trace_1_45_1_2',
#  'Trace_1_46_1_2',
#  'Trace_1_47_1_1',
#  'Trace_1_47_1_2',
#  'Trace_1_48_1_1',
#  'Trace_1_48_1_2',
#  'Trace_1_49_1_1',
#  'Trace_1_49_1_2',
#  'Trace_1_50_1_2',
#  'Trace_1_51_1_1',
#  'Trace_1_51_1_2',
#  'Trace_1_52_1_1',
#  'Trace_1_52_1_2',
#  'Trace_1_53_1_1',
#  'Trace_1_53_1_2']