# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:10:40 2024

Ok finally making a concise python file to go through both biologic files and metrohm files

Options I want to include:
    - Plot over blank
    - If I want to just plot last cycles or all
    - Conversion of potential scale
    - 

e.g. plot(figure = file, blank = y/n, allCycles = y/n, potConv = float)

@author: Elliot
"""

import os
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.colors as colors

Cmap = plt.get_cmap('tab20')
#locator = MultipleLocator(base = 0.5, offset= 0.0)


def read_file(file):
    with open(file, 'rb') as check:
        for line in check:
            if line.split()[0] == b'Index,Time':
                fileType = 'metrohm'
            else:
                fileType = 'biologic'
    
    if fileType == 'biologic':
        with open(file, 'rb') as file:
            for line in file:
                check = line.split()
                if (len(check) > 0) and (check[0].decode('utf-8') == 'Nb'):
                    rowSkip = int(check[-1]) 
                    break
        
    
    elif fileType == 'metrohm':
        defg
    
    else: 
        print("This file doesn't seem to be metrohm or biologic...")
        return
    
    return

def metrohmRead(met_file):
    read_file = pd.read_csv(met_file) #reading in csv file
    read_file["Cycle Number"] = np.nan
    
    potential = read_file['Potential applied (V)']
    direction = potential.diff()
    cycle_num = 1

    for i in range(min(read_file["Index"])-1, max(read_file["Index"])):
        if i == min(read_file["Index"]) or max(read_file["Index"]):
            read_file.loc[i, "Cycle Number"] = cycle_num
            continue
        if potential[i] == 0.0 and direction[i-1] < 0 and direction[i+1] > 0:
            cycle_num += 1
        if potential[i] == 0.0 and direction[i-1] > 0 and direction[i+1] > 0:
            cycle_num += 1
        read_file.loc[i, "Cycle Number"] = cycle_num
            
    return read_file

def biologicRead(bio_file):
    
    
    return

def folder_process(folder_path):
    files = [f for f in os.listdir(folder_path) if f.endswith('.csv' or '.txt' or '.mpt')]
    
    return