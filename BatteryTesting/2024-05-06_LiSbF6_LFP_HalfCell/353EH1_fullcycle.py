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

os.chdir("C:/Users/Elliot/SynologyDrive/Research - Elliot Howell/UQAM - Battery/Half Cells/2024-05-06_LiSbF6_LFP_HalfCells/353EH1/")

Cov10_353EH1 = pd.read_csv("353EH1_FullCycling_01_GCPL_C02.txt", 
                           sep = '\t', skiprows = 78, encoding_errors = 'ignore')
Cov5_353EH1 = pd.read_csv("353EH1_FullCycling_02_GCPL_C02.txt", 
                           sep = '\t', skiprows = 78, encoding_errors = 'ignore')
Cov2_353EH1 = pd.read_csv("353EH1_FullCycling_03_GCPL_C02.txt", 
                           sep = '\t', skiprows = 78, encoding_errors = 'ignore')
C_353EH1 = pd.read_csv("353EH1_FullCycling_04_GCPL_C02.txt", 
                           sep = '\t', skiprows = 78, encoding_errors = 'ignore')
C2_353EH1 = pd.read_csv("353EH1_FullCycling_05_GCPL_C02.txt", 
                           sep = '\t', skiprows = 78, encoding_errors = 'ignore')
FCov10_353EH1_1of2 = pd.read_csv("353EH1_FullCycling_06_GCPL_C02.txt", 
                           sep = '\t', skiprows = 78, encoding_errors = 'ignore')
FCov10_353EH1_2of2 = pd.read_csv("353EH1_FullCycling_restart_06_GCPL_C02.txt", 
                           sep = '\t', skiprows = 78, encoding_errors = 'ignore')

headers = list(Cov10_353EH1)
'''
Here is the list of header names, for ease of calling specific names we can just have:
    headers[7] = time (s)
    headers[11] = potential (V)
    headers[19] = current (mA)
    
'''

plt.figure()

