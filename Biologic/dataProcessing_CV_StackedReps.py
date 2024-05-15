# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 12:57:35 2023

@author: Elliot
"""

import os
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# # Folder path containing the .mpt files ** windows **
folder_path = r'C:\Users\Elliot\SynologyDrive\Research - Elliot Howell\Durbis CV Measurements\TPB-OMe\DCM\To Use\vsFc'
# # Folder path containing the .mpt files ** mac **
#folder_path = r'/Users/elliothowell/SynologyDrive/Research - Elliot Howell/UQAM - Battery/LiSbF6 Salt/15-02-2024 Blank ECDMC/Bottle A'

os.chdir(folder_path)

# Get a list of files ending with '.mpt' in the folder ** WINDOWS **
txt_files = [f for f in os.listdir(folder_path) if (f.endswith('.mpt') or f.endswith('.txt'))]

# User inputs to determine Reference electrode, analyte, etc.
refElec = input("What is the reference electrode used?\n")
analyte = input("And what is the analyte of interest?\n")

fcCond = input("Do you want to compare potential to Fc? (Y = Yes / N = No)\n")

if fcCond == 'Y':
    refElec = "Fc/$Fc^+$"
    fcHalf = input("Ferrocene half potential?\n")
    fcHalf = float(fcHalf)

# Process each .mpt file
for file_name in txt_files:
    file_path = os.path.join(folder_path, file_name)
    
    #holder for scan rate number
    scanRate = 0 
    scanUnit = ''
    n = 0
    with open(file_path, 'rb') as file:
        for line in file:
            check = line.split()
            n +=1
            if scanRate == 0:
                if (len(check) > 0) and (check[0].decode('utf-8') == 'dE/dt'):
                    #setting scan rate
                    scanRate = check[1]
            if scanUnit == '':    
                if (len(check) > 0) and (check[1].decode('utf-8') == 'unit'):
                    scanUnit = check[2]
            if (len(check) > 0) and (check[0].decode('utf-8') == 'Nb'):
                rowSkip = int(check[-1]) 
            if (len(check) > 0) and (check[0].decode('utf-8') == 'Acquisition'):
                CVtime = re.findall(r'[^.]+', check[-1].decode('utf-8'))[0]
                CVtimeDot = re.sub(':', '.', CVtime)
   
    # Read the EC-Lab CV data
    cv_data = pd.read_csv(file_path, encoding = 'unicode_escape', encoding_errors= 'ignore', 
                          skiprows = rowSkip-1, sep = '\t')
    # print(file_name)
    # print(rowSkip)
    # print(cv_data.head())
    
    # Breaking the data into replicates (i.e. list of replicate dataframes)
    repList = []
    for i in range(len(np.unique(cv_data['cycle number']))):
        repList.append(cv_data.loc[cv_data['cycle number'] == i+1])
    
    # Checking if anything will actually be plotted, if not just continue
    if len(repList) == 1:
        continue
    
    # now creating plot of stacked reps with first cycle included
    title = analyte + ' versus ' + refElec + ' at a scan rate of ' + scanRate.decode('utf-8') + " " + scanUnit.decode('utf-8')
    plt.figure('CV for ' + file_name + 'with first cycle')
    plt.title('CV of ' + title + '\n' + '\nfile: ' + file_name + '\nAt ' + CVtime)
    plt.xlabel('Potential (V vs. ' + refElec + ')')
    plt.ylabel('Current (\u03BCA)')
    
    # Getting file name information
    file_name_list = file_name.split(sep = '_')
    
    #plot of each rep
    for i in range(len(repList)):
        #if i == 0:
        #    continue
        if fcCond =='Y':
            potential = repList[i]['Ewe/V'] - fcHalf
        else: 
            potential = repList[i]['Ewe/V']
        
        current = repList[i]['<I>/mA'] * 1000
        
        if i == 0:
            plt.plot(potential, current, label = "Initial Scan")
        else:
            plt.plot(potential, current, label = 'Scan ' + str(i))
    

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig('firstCycle_fig_' + CVtimeDot + '_' + file_name[:-4] + '.tif', dpi = 300, bbox_inches = 'tight')
    plt.show()
    
    # now creating plot of stacked reps without first cycle
    title = analyte + ' versus ' + refElec + ' at a scan rate of ' + scanRate.decode('utf-8') + " " + scanUnit.decode('utf-8')
    plt.figure('CV for ' + file_name)
    plt.title('CV of ' + title + '\n' + '\nfile: ' + file_name + '\nAt ' + CVtime)
    plt.xlabel('Potential (V vs. ' + refElec + ')')
    plt.ylabel('Current (\u03BCA)')
    
    #plot of each rep
    for i in range(len(repList)):
        if i == 0:
            continue
        if fcCond == 'Y':
            potential = repList[i]['Ewe/V'] - fcHalf
        else: 
            potential = repList[i]['Ewe/V']
        current = repList[i]['<I>/mA'] * 1000
        
        plt.plot(potential, current, label = 'Scan ' + str(i))
    
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig('fig_' + CVtimeDot + '_' + file_name[:-4] + '.tif', dpi = 300, bbox_inches = 'tight')
    plt.show()
    
