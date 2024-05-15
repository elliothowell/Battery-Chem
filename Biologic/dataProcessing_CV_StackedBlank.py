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
# folder_path = r'C:\Users\Elliot\SynologyDrive\Research - Elliot Howell\Durbis CV Measurements\TPB-OMe\DCM\To Use\vsFc'
# # Folder path containing the .mpt files ** mac **
#folder_path = r'/Users/elliothowell/SynologyDrive/Research - Elliot Howell/UQAM - Battery/LiSbF6 Salt/15-02-2024 Blank ECDMC/Bottle A'

figure = "C:/Users/Elliot/SynologyDrive/Research - Elliot Howell/Durbis CV Measurements/TPB-OMe/DCM/To Use/21-03-2024_100mvs_-500-1000mV_TPBOMe_3rep_C02.txt"
blank = "C:/Users/Elliot/SynologyDrive/Research - Elliot Howell/Durbis CV Measurements/TPB-OMe/DCM/To Use/21-03-2024_100mvs_-500-1000mV_Blk_3rep_postTPBOMe_C02.txt"
saveLoc = r'C:\Users\Elliot\SynologyDrive\Research - Elliot Howell\Durbis CV Measurements\TPB-OMe\DCM\To Use'

n = 0
figScanRate = 0
figScanUnit = ''
with open(figure, 'rb') as file:
    for line in file:
        check = line.split()
        n +=1
        if figScanRate == 0:
            if (len(check) > 0) and (check[0].decode('utf-8') == 'dE/dt'):
                #setting scan rate
                figScanRate = check[1]
        if figScanUnit == '':    
            if (len(check) > 0) and (check[1].decode('utf-8') == 'unit'):
                figScanUnit = check[2]
        if (len(check) > 0) and (check[0].decode('utf-8') == 'Nb'):
            figRowSkip = int(check[-1]) 

n = 0
blkScanRate = 0
blkScanUnit = ''
with open(figure, 'rb') as file:
    for line in file:
        check = line.split()
        n +=1
        if blkScanRate == 0:
            if (len(check) > 0) and (check[0].decode('utf-8') == 'dE/dt'):
                #setting scan rate
                blkScanRate = check[1]
        if blkScanUnit == '':    
            if (len(check) > 0) and (check[1].decode('utf-8') == 'unit'):
                blkScanUnit = check[2]
        if (len(check) > 0) and (check[0].decode('utf-8') == 'Nb'):
            blkRowSkip = int(check[-1]) 

figure_data = pd.read_csv(figure, encoding = 'unicode_escape', encoding_errors= 'ignore', 
                      skiprows = figRowSkip-1, sep = '\t')

blank_data = pd.read_csv(blank, encoding = 'unicode_escape', encoding_errors= 'ignore', 
                      skiprows = blkRowSkip-1, sep = '\t')

figure_data = figure_data.loc[figure_data['cycle number'] == len(np.unique(figure_data['cycle number']))]

blank_data = blank_data.loc[blank_data['cycle number'] == len(np.unique(blank_data['cycle number']))]

# User inputs to determine Reference electrode, analyte, etc.
refElec = input("What is the reference electrode used?\n")
analyte = input("And what is the analyte of interest?\n")

fcCond = input("Do you want to compare potential to Fc? (Y = Yes / N = No)\n")

if fcCond == 'Y':
    refElec = "Fc/$Fc^+$"
    fcHalf = input("Ferrocene half potential?\n")
    fcHalf = float(fcHalf)
    blankPot = blank_data['Ewe/V'] - fcHalf
    figPot = figure_data['Ewe/V'] - fcHalf
else:
    blankPot = blank_data['Ewe/V']
    figPot = figure_data['Ewe/V']

    # now creating plot of stacked reps without first cycle
title = analyte + ' versus ' + refElec + ' at a scan rate of ' + figScanRate.decode('utf-8') + " " + figScanUnit.decode('utf-8')
plt.figure()
plt.title('CV stacked over blank')
plt.xlabel('Potential (V vs. ' + refElec + ')')
plt.ylabel('Current (\u03BCA)')

blankCurr = blank_data['<I>/mA'] * 1000

figCurr = figure_data['<I>/mA'] * 1000

plt.plot(blankPot, blankCurr, label = "Blank")
plt.plot(figPot, figCurr, label = "Compound")

plt.savefig(saveLoc + '/' + analyte + 'overBlk' + '.tif', dpi = 300, bbox_inches = 'tight')
# # Process each .mpt file
# for file_name in txt_files:
#     file_path = os.path.join(folder_path, file_name)
    
#     #holder for scan rate number
#     scanRate = 0 
#     scanUnit = ''
#     n = 0
#     with open(file_path, 'rb') as file:
#         for line in file:
#             check = line.split()
#             n +=1
#             if scanRate == 0:
#                 if (len(check) > 0) and (check[0].decode('utf-8') == 'dE/dt'):
#                     #setting scan rate
#                     scanRate = check[1]
#             if scanUnit == '':    
#                 if (len(check) > 0) and (check[1].decode('utf-8') == 'unit'):
#                     scanUnit = check[2]
#             if (len(check) > 0) and (check[0].decode('utf-8') == 'Nb'):
#                 rowSkip = int(check[-1]) 
#             if (len(check) > 0) and (check[0].decode('utf-8') == 'Acquisition'):
#                 CVtime = re.findall(r'[^.]+', check[-1].decode('utf-8'))[0]
#                 CVtimeDot = re.sub(':', '.', CVtime)
   
#     # Read the EC-Lab CV data
#     cv_data = pd.read_csv(file_path, encoding = 'unicode_escape', encoding_errors= 'ignore', 
#                           skiprows = rowSkip-1, sep = '\t')
#     # print(file_name)
#     # print(rowSkip)
#     # print(cv_data.head())
    
#     # Breaking the data into replicates (i.e. list of replicate dataframes)
#     repList = []
#     for i in range(len(np.unique(cv_data['cycle number']))):
#         repList.append(cv_data.loc[cv_data['cycle number'] == i+1])
    
#     # Checking if anything will actually be plotted, if not just continue
#     if len(repList) == 1:
#         continue
    
#     # now creating plot of stacked reps with first cycle included
#     title = analyte + ' versus ' + refElec + ' at a scan rate of ' + scanRate.decode('utf-8') + " " + scanUnit.decode('utf-8')
#     plt.figure('CV for ' + file_name + 'with first cycle')
#     plt.title('CV of ' + title + '\n' + '\nfile: ' + file_name + '\nAt ' + CVtime)
#     plt.xlabel('Potential (V vs. ' + refElec + ')')
#     plt.ylabel('Current (\u03BCA)')
    
#     # Getting file name information
#     file_name_list = file_name.split(sep = '_')
    
#     #plot of each rep
#     for i in range(len(repList)):
#         #if i == 0:
#         #    continue
#         if fcCond =='Y':
#             potential = repList[i]['Ewe/V'] - fcHalf
#         else: 
#             potential = repList[i]['Ewe/V']
        
#         current = repList[i]['<I>/mA'] * 1000
        
#         if i == 0:
#             plt.plot(potential, current, label = "Initial Scan")
#         else:
#             plt.plot(potential, current, label = 'Scan ' + str(i))
    

#     plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#     plt.savefig('firstCycle_fig_' + CVtimeDot + '_' + file_name[:-4] + '.tif', dpi = 300, bbox_inches = 'tight')
#     plt.show()
    
#     # now creating plot of stacked reps without first cycle
#     title = analyte + ' versus ' + refElec + ' at a scan rate of ' + scanRate.decode('utf-8') + " " + scanUnit.decode('utf-8')
#     plt.figure('CV for ' + file_name)
#     plt.title('CV of ' + title + '\n' + '\nfile: ' + file_name + '\nAt ' + CVtime)
#     plt.xlabel('Potential (V vs. ' + refElec + ')')
#     plt.ylabel('Current (\u03BCA)')
    
#     #plot of each rep
#     for i in range(len(repList)):
#         if i == 0:
#             continue
#         if fcCond == 'Y':
#             potential = repList[i]['Ewe/V'] - fcHalf
#         else: 
#             potential = repList[i]['Ewe/V']
#         current = repList[i]['<I>/mA'] * 1000
        
#         plt.plot(potential, current, label = 'Scan ' + str(i))
    
#     plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#     plt.savefig('fig_' + CVtimeDot + '_' + file_name[:-4] + '.tif', dpi = 300, bbox_inches = 'tight')
#     plt.show()
    
