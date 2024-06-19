# -*- coding: utf-8 -*-
"""
Created on Fri May 24 16:07:34 2024

@author: Elliot
"""

import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import datetime as dt
import numpy as np
import scipy.signal as ss
from General_txt_file_reader import read_file

os.chdir("/Users/elliothowell/SynologyDrive/Research - Elliot Howell/UQAM - Battery/Half Cells/2024-05-06_LiSbF6_LFP_HalfCells/353EH1/")

Cov10_353EH1 = read_file([r"353EH1_FullCycling_01_GCPL_C02.txt"])

Cov5_353EH1 = read_file([r"353EH1_FullCycling_02_GCPL_C02.txt"])

Cov2_353EH1 = read_file([r"353EH1_FullCycling_03_GCPL_C02.txt"])

C_353EH1 = read_file([r"353EH1_FullCycling_04_GCPL_C02.txt"])

C2_353EH1 = read_file([r"353EH1_FullCycling_05_GCPL_C02.txt"])

FCov10_353EH1_1of2 = read_file([r"353EH1_FullCycling_06_GCPL_C02.txt"])

FCov10_353EH1_2of2 = read_file([r"353EH1_FullCycling_restart_06_GCPL_C02.txt"])

'''


Here is the list of header names, for ease of calling specific names we can just have:
    headers[7] = time (s)
    headers[11] = potential (V)
    headers[19] = current (mA)
    headers[21] = Q discharge
    headers[22] = Q charge
    
    mode == 1 and ox/red == 1 -> charge
    mode == 1 and ox/red == 0 -> discharge
    
'''

# C/10 first cycles analysis

plt.figure('chgCurve')
plt.title('Charge Curves for 353EH1')
plt.xlabel('Capacity (mAh)')
plt.ylabel('Voltage (V)')

xChg = Cov10_353EH1[0][[Cov10_353EH1[0][:,0] == Cov10_353EH1[0][:,1]][0],22]
yChg = Cov10_353EH1[0][[Cov10_353EH1[0][:,0] == Cov10_353EH1[0][:,1]][0],11]

plt.scatter(xChg, yChg, label = 'testing')
    
plt.show()


'''
Going to try to automate going through cycles, etc.

so we know that:
    mode == 1 and ox/red == 1 -> charge
    mode == 1 and ox/red == 0 -> discharge
    
Since we have that we can cycle through Ns (index 6) which 
'''
toAnalyze = Cov10_353EH1[0]

chgCyc = 1
dchgCyc = 1

for i in np.unique(toAnalyze[:,6]):
    print(i)
    if np.any(toAnalyze[toAnalyze[:,6] == i][:,0] == 1) & np.any(toAnalyze[toAnalyze[:,6] == i][:,1] == 1):
        # should be on charge cycle
        print("charge")    
    elif np.any(toAnalyze[toAnalyze[:,6] == i][:,0] == 1) & np.any(toAnalyze[toAnalyze[:,6] == i][:,1] == 0):
        # should be on discharge cycle
        print("discharge")


