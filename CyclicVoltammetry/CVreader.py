# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:10:40 2024

Ok finally making a concise python file to go through both biologic files and metrohm files

Options I want to include:
    - Plot over blank (maybe future iterations...)
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

Cmap = plt.get_cmap('tab10')
#locator = MultipleLocator(base = 0.5, offset= 0.0)

columnNames = ['time(s)', 'potential(V)', 'current(A)', 'cycle number'] #column names I want

def metrohmRead(met_file):
    read_file = pd.read_csv(met_file) #reading in csv file
    read_file["Cycle Number"] = np.nan
    
    potential = read_file['Potential applied (V)']
    direction = potential.diff()
    cycle_num = 1
    
    # this is just adding a column for cycle number as metrohm doesn't do it (:
    for i in range(len(read_file)):
        if i == 0 or i == len(read_file)-1:
            read_file.loc[i, "Cycle Number"] = cycle_num
            continue
        elif potential[i] == 0.0 and (direction[i-1] < 0 and direction[i+1] > 0):
            cycle_num += 1
            read_file.loc[i, "Cycle Number"] = cycle_num
        elif potential[i] == 0.0 and (direction[i-1] > 0 and direction[i+1] > 0):
            cycle_num += 1
            read_file.loc[i, "Cycle Number"] = cycle_num
        else:
            read_file.loc[i, "Cycle Number"] = cycle_num
    
    #extract only columns I want
    read_file = read_file[['Time (s)', 'Potential applied (V)', 'WE.Current (A)', 'Cycle Number']]
    
    #renaming columns to general form
    for i in range(len(read_file.columns)):
        read_file = read_file.rename(columns={read_file.columns[i]:columnNames[i]})
    
    read_file[columnNames[-1]] = read_file[columnNames[-1]].astype(int)
    
    return read_file

def biologicRead(bio_file):
    
    #finding number for row skip
    with open(bio_file, 'rb') as file:
        for line in file:
            check = line.split()
            if (len(check) > 0) and (check[0].decode('utf-8') == 'Nb'):
                rowSkip = int(check[-1]) 
                break
    
    #reading in file based on rowSkip
    read_file = pd.read_csv(bio_file, encoding = 'unicode_escape', encoding_errors= 'ignore', 
                          skiprows = rowSkip-1, sep = '\t')
    
    # extracting only columns I want
    read_file = read_file[['time/s', 'Ewe/V', '<I>/mA', 'cycle number']]
    
    # this is a rough hard code but want to translate mA to A
    read_file['<I>/mA'] = read_file['<I>/mA']/1000
    
    # changing column names to general form
    for i in range(len(read_file.columns)):
        read_file = read_file.rename(columns={read_file.columns[i]:columnNames[i]})
    
    read_file[columnNames[-1]] = read_file[columnNames[-1]].astype(int)
    
    return read_file

# this just takes in the file you want to read and decides if its biologic or metrohm
# which then calls which specific function you need
def read_file(file):
    with open(file, 'rb') as check:
        for line in check:
            if len(line.split()) > 0:
                if line.split()[0] == b'Index,Time':
                    fileType = 'metrohm'
                    break
                elif line.split()[0] == b'EC-Lab':
                    fileType = 'biologic'
                    break
                else:
                    continue
            else:
                continue
    
    if fileType == 'biologic':
        cvFile = biologicRead(file)
    
    elif fileType == 'metrohm':
        cvFile = metrohmRead(file)
    
    else: 
        print("This file doesn't seem to be metrohm or biologic...")
        return
    
    return cvFile

# kind of want this to return both stacked reps and 
def cv_plotAndSave(folderPath, fcPot, saveFigs):
    
    cvFiles = [f for f in os.listdir(folderPath) if f.endswith('.csv') or f.endswith('.txt') or f.endswith('.mpt')]
    figFolder = os.path.join(folderPath, 'newestFigures')
    if not os.path.exists(figFolder):
        os.makedirs(figFolder)
    
    # in here want to make plots for last cycle as well as stacked cycles
    for file in cvFiles:
        
        file_path = os.path.join(folderPath, file)
        
        #make a dataframe of the file
        fileDF = read_file(file_path)
        
        #first making the stacked cycles
        figStack, axStack = plt.subplots(figsize = (10, 6))
        
        # fig.figure(figsize=(10, 6))
        
        for i in range(1, max(fileDF[columnNames[-1]].unique()) + 1):
            
            if fcPot is not None:
                potential = fileDF.loc[fileDF[columnNames[-1]] == i, columnNames[1]] - fcPot
            else: 
                potential = fileDF.loc[fileDF[columnNames[-1]] == i, columnNames[1]]
            
            # multiplying by 1*10^6 to convert to uA
            current = fileDF.loc[fileDF[columnNames[-1]] == i, columnNames[2]] * (1*10**6)

            axStack.plot(potential, current, label = 'Cycle ' + str(i), c = Cmap(i-1))
            
        axStack.set_xlabel(f'Potential (vs {"Fc/$Fc^+$" if fcPot is not None else "Ag/AgCl"})') 
        axStack.set_ylabel('Current ($\mu$A)')
        axStack.legend()
        axStack.set_title(f'All cycles for {file}')
        
        if saveFigs:    
            figStack.savefig(os.path.join(figFolder, f'All_Cycles_{file}.tif'), 
                        dpi = 300, bbox_inches = 'tight')
        
        figLast, axLast = plt.subplots(figsize = (10,6))
        
        axLast.plot(potential, current, c = Cmap(0))
        axLast.set_xlabel(f'Potential (vs {"Fc/$Fc^+$" if fcPot is not None else "Ag/AgCl"})') 
        axLast.set_ylabel('Current ($\mu$A)')
        axLast.set_title(f'Last cycle for {file}')
        
        if saveFigs:    
            figLast.savefig(os.path.join(figFolder, f'Last_Cycle_{file}.tif'), 
                        dpi = 300, bbox_inches = 'tight')
        
# def plotAllCycles(dataFrame):
    

# def plotLastCycle(dataFrame):


