# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 13:13:47 2024

@author: Elliot
"""

import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import datetime as dt
import numpy as np
import scipy.signal as ss

# # Folder path containing the .mpt files ** windows **
folderPath = 'C:\\Users\\Elliot\\SynologyDrive\\Research - Elliot Howell\\UQAM - Battery\\23-10-2023_HalfCells'

## Folder path containing the .mpt files ** mac **
#folderPath = '/Users/elliothowell/SynologyDrive/Research - Elliot Howell/UQAM - Battery/23-10-2023_HalfCells'
os.chdir(folderPath)

# Get a list of files ending with '.mpt' in the folder ** WINDOWS **
csvFiles = [f for f in os.listdir(folderPath) if f.endswith('.csv')]

# Creating a empty list to store dataFrames and names of 
dataList = []
fileNames = []

# List of active material weight
actMatWeight = {'LFP_3-1.csv': 4.00816, 'LFP_3-2.csv': 4.3376, 'LFP_3-3.csv': 5.1376}

# Creating 3 counting variables that start at zero
i = 0
j = 0
k = 0
secInHour = 3600

# Adding each data file to a list/creating dataFrames
for fileName in csvFiles:
    #creating file path and names for each data set
    filePath = os.path.join(folderPath, fileName)
    fileNames.append(fileName[:-4])
    cvData = pd.read_csv(filePath)
    dataList.append(cvData)
    
    #convert total time to seconds and will also convert cycle time to seconds
    dataList[-1]['Total Seconds'] = (pd.to_timedelta(dataList[-1]['Total Time']).dt.total_seconds())
    dataList[-1]['Cycle Seconds'] = (pd.to_timedelta(dataList[-1]['Time']).dt.total_seconds())

#going through each data frame (each battery) and creating plots
for i in range(len(dataList)):
    
    # creating empty list for cycles and x/y variables
    cycleList = []
    ccChgList = []
    ccDchgList = []
    xChg = []
    yChg = []
    xDchg = []
    yDchg = []
    
    #making a plot for ceach plot I want to make (these should go before each for loop)
    # plt.figure('totCap'+str(i))
    # plt.figure('dQdE'+str(i))
    
    ###########################################################################
    #This should result in plots for voltage against total capacity 
    plt.figure('capCurve'+str(i))
    plt.title('Total Capacity Curves for sample '+fileNames[i])
    plt.xlabel('Total Capacity (mAh/g)') 
    plt.ylabel('Voltage (V)')
    
    
    #variables for totCap
    integration = []
    totCap = 0
    for j in range(len(dataList[i]['DataPoint']) - 1):
        
        avgCurr = 0.5*((dataList[i].loc[j+1, 'Current(A)']) + (dataList[i].loc[j, 'Current(A)']))*1000
        
        timeDiff = ((dataList[i].loc[j+1, 'Total Seconds']) - (dataList[i].loc[j, 'Total Seconds']))/secInHour
    
        totCap = (avgCurr * timeDiff) / (actMatWeight[csvFiles[i]]/1000)
        
        integration.append(totCap)
        
        dataList[i].loc[j, 'totCap'] = np.cumsum(integration)[-1]
        
    totCap = dataList[i]['totCap']
    voltage = dataList[i]['Voltage(V)']
    
    plt.plot(totCap, voltage)
    plt.savefig('totCapCurve_'+fileNames[i]+'.tif', dpi = 300, bbox_inches = 'tight')
    plt.show()
    
    ###########################################################################
    #Creating cycle loop to make specific lists
    for j in range(max(dataList[i]['Cycle Index'])):  
        # this only needs to be performed this one time as this is performed first for each dataframe
        cycleList.append(dataList[i][dataList[i]['Cycle Index'] == j+1])
        
        ccChgList.append(cycleList[j][cycleList[j]['Step Type'] == 'CC Chg'])
        
        ccDchgList.append(cycleList[j][cycleList[j]['Step Type'] == 'CC DChg'])
        
    ###########################################################################
    #creating figure for charge curve, overlapping cycles
    plt.figure('chgCurve'+str(i))
    plt.title('Charge Curves for sample '+fileNames[i])
    plt.xlabel('Capacity (mAh/g)') #mAh as multiplied cap by 1000 (1000 mA/A)
    plt.ylabel('Voltage (V)')
    #This should result in stacked cycles for charge/discharge
    for j in range(max(dataList[i]['Cycle Index'])):  
        # # this only needs to be performed this one time as this is performed first for each dataframe
        # cycleList.append(dataList[i][dataList[i]['Cycle Index'] == j+1])
        
        # ccChgList.append(cycleList[j][cycleList[j]['Step Type'] == 'CC Chg'])
        
        # Multiplying cap by 1000 to turn to mAh and dividing by actMatWeight index /1000 to find mAh/g
        xChg = (ccChgList[j]['Capacity(Ah)'] * 1000)/(actMatWeight[csvFiles[i]]/1000)
        yChg = ccChgList[j]['Voltage(V)']
        
        plt.plot(xChg, yChg, label = 'Charge cycle '+str(j+1))
        
    plt.legend(loc = 'lower right')
    plt.savefig('chargeCurve_'+fileNames[i]+'.tif', dpi = 300, bbox_inches = 'tight')
    plt.show()
    
    ###########################################################################
    #creating figure for charge curve for HOURS, overlapping cycles
    plt.figure('chgCurve'+str(i))
    plt.title('Charge Curves for sample '+fileNames[i])
    plt.xlabel('Hours') #mAh as multiplied cap by 1000 (1000 mA/A)
    plt.ylabel('Voltage (V)')
    #This should result in stacked cycles for charge/discharge
    for j in range(max(dataList[i]['Cycle Index'])):  
        # # this only needs to be performed this one time as this is performed first for each dataframe
        # cycleList.append(dataList[i][dataList[i]['Cycle Index'] == j+1])
        
        # ccChgList.append(cycleList[j][cycleList[j]['Step Type'] == 'CC Chg'])
        
        # Multiplying cap by 1000 to turn to mAh and dividing by actMatWeight index /1000 to find mAh/g
        xChg = (ccChgList[j]['Cycle Seconds']) / secInHour
        yChg = ccChgList[j]['Voltage(V)']
        
        plt.plot(xChg, yChg, label = 'Charge cycle '+str(j+1))
        
    plt.legend(loc = 'lower right')
    plt.savefig('chargeCurveHours_'+fileNames[i]+'.tif', dpi = 300, bbox_inches = 'tight')
    plt.show()
    
    ###########################################################################
    #creating figure for discharge curve
    plt.figure('dchgCurve'+str(i))
    plt.title('Discharge Curves for sample '+fileNames[i])
    plt.xlabel('Capacity (mAh/g)') #mAh as multiplied cap by 1000 (1000 mA/A)
    plt.ylabel('Voltage (V)')
    #This should result in stacked cycles for charge/discharge
    for j in range(max(dataList[i]['Cycle Index'])):  
        
        # ccDchgList.append(cycleList[j][cycleList[j]['Step Type'] == 'CC DChg'])
        
        # Multiplying cap by 1000 to turn to mAh
        xDchg = (ccDchgList[j]['Capacity(Ah)'] * 1000)/(actMatWeight[csvFiles[i]]/1000)
        yDchg = ccDchgList[j]['Voltage(V)']

        plt.plot(xDchg, yDchg, label = 'Discharge cycle '+str(j+1))
        
    plt.legend()
    plt.savefig('dischargeCurve_'+fileNames[i]+'.tif', dpi = 300, bbox_inches = 'tight')
    plt.show()
    
    ###########################################################################
    #creating figure for discharge curve against HOURS
    plt.figure('dchgCurve'+str(i))
    plt.title('Discharge Curves for sample '+fileNames[i])
    plt.xlabel('Hours') #mAh as multiplied cap by 1000 (1000 mA/A)
    plt.ylabel('Voltage (V)')
    #This should result in stacked cycles for charge/discharge
    for j in range(max(dataList[i]['Cycle Index'])):  
        
        # ccDchgList.append(cycleList[j][cycleList[j]['Step Type'] == 'CC DChg'])
        
        # Multiplying cap by 1000 to turn to mAh
        xDchg = ccDchgList[j]['Cycle Seconds'] / secInHour
        yDchg = ccDchgList[j]['Voltage(V)']

        plt.plot(xDchg, yDchg, label = 'Discharge cycle '+str(j+1))
        
    plt.legend()
    plt.savefig('dischargeCurveHours_'+fileNames[i]+'.tif', dpi = 300, bbox_inches = 'tight')
    plt.show()
    
    ###########################################################################
    # #creating figure for diff cap curve
    # plt.figure('dQdE curve'+str(i))
    # plt.title('dQdE Curves for sample '+fileNames[i])
    # plt.xlabel('Voltage (V)') 
    # plt.ylabel('Differential Capacity (mAh/V)')
    # #This should result in stacked cycles for differential capacity
    # for j in range(max(dataList[i]['Cycle Index'])):  
        
    #     voltage =  cycleList[j]['Voltage(V)']
    #     dQdE = cycleList[j]['dQ/dV(mAh/V)']
    #     smoothdQdE = ss.savgol_filter(dQdE,12,2)

    #     plt.plot(voltage, smoothdQdE, label = 'Cycle '+str(j+1))
        
    # plt.legend(loc = 'lower left')
    # plt.savefig('diffCapCurve_'+fileNames[i]+'.tif', dpi = 300, bbox_inches = 'tight')
    # plt.show()
    
    ###########################################################################
    #This should result in plots for manually calculated differential cap
    plt.figure('capCurve'+str(i))
    plt.title('Differential Capacitance Curves for sample '+fileNames[i])
    plt.xlabel('Voltage (V)') 
    plt.ylabel('Differential Capacitance (dQ/dV)')
    for j in range(max(dataList[i]['Cycle Index'])):    
        cycleData = cycleList[j].reset_index()
        
        slopeList = []
        # Integrating current/time to determine dQ/dE
        for k in range(len(cycleData)-1): 
            
            #rise is tot cap difference and run is voltage difference
            rise = (cycleData['totCap'][k+1] - cycleData['totCap'][k])
            run = abs((cycleData['Voltage(V)'][k+1] - cycleData['Voltage(V)'][k]))
            
            # aFinding slope
            slope = rise / run
            slopeList.append(slope)
            
            
        #window = int(len(cycleData)/5)
        smoothSlope = ss.savgol_filter(slopeList, 50, 3)
        dQdE = smoothSlope
        voltage = cycleData['Voltage(V)'][:-1]
        
        plt.plot(voltage, dQdE, label = 'Cycle '+str(j+1))
        
    plt.legend(loc = 'lower left')
    plt.savefig('difCapCurve_'+fileNames[i]+'.tif', dpi = 300, bbox_inches = 'tight')
    plt.show()
    
    ###########################################################################
    #This should result in plots for current/voltage against time 
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    plt.title('Applied Current and Voltage response in respect to time for sample '+fileNames[i])
    
    ax1.set_xlabel('Hours')
    ax1.set_ylabel('Voltage (V)')
    ax1.tick_params(axis = 'y')
    
    ax2.set_ylabel('Current (mA)')
    ax2.tick_params(axis = 'y', labelcolor = 'tab:red')
    current = dataList[i]['Current(A)'] * 1000
    time = dataList[i]['Total Seconds'] / secInHour
    ax2.plot(time, current, color = 'tab:red')
    
    for j in range(max(dataList[i]['Cycle Index'])):    
        
        # current in mA
        # current = cycleList[j]['Current(A)'] * 1000
        voltage = cycleList[j]['Voltage(V)']
        time = cycleList[j]['Total Seconds'] / secInHour
        
        ax1.plot(time, voltage, label = 'Cycle '+str(j+1))
        # ax2.plot(time, current, label = 'Cycle '+str(j+1))
    
    ax1.legend(loc = 'upper center', bbox_to_anchor=(0.5, -0.15), ncol = j+1)
    plt.savefig('currVolTimCurve_'+fileNames[i]+'.tif', dpi = 300, bbox_inches = 'tight')
    plt.show()
    
    ###########################################################################
    #This should result in plots for voltage against total capacity 
    # plt.figure('capCurve'+str(i))
    # plt.title('Total Capacity Curves for sample '+fileNames[i])
    # plt.xlabel('Total Capacity (mAh/g)') 
    # plt.ylabel('Voltage (V)')
    # for j in range(max(dataList[i]['Cycle Index'])):    
    #     cycleData = cycleList[j].reset_index()
    #     # cycleData = cycleList[j].reset_index()
        
    #     integration = []
    #     totCap = 0
    #     # Integrating current/time to determine dQ/dE
    #     for k in range(len(cycleData)-1): 
            
    #         # dividing current values by 1000 to get mA from A
    #         avgCurr = 0.5*((cycleData.loc[k+1, 'Current(A)'])/1000 + (cycleData.loc[k, 'Current(A)'])/1000) 
    #         # and dividing specific times by 3600 to get hour diff
    #         timeDiff = ((cycleData.loc[k+1, 'Total Seconds']/secInHour) - (cycleData.loc[k, 'Total Seconds']/secInHour))
    #         totCap = (avgCurr * timeDiff) / (actMatWeight[i]/1000)
    #         integration.append(totCap)
            
    #         cycleData.loc[k, 'totCap'] = np.cumsum(integration)[-1]
        
    #     totCap = cycleData['totCap'] / (actMatWeight[i]/1000)
    #     voltage = cycleData['Voltage(V)']
        
    #     plt.plot(totCap, voltage, label = 'Cycle '+str(j+1))
        
    # plt.legend()
    # plt.savefig('totCapCurve_'+fileNames[i]+'.tif', dpi = 300, bbox_inches = 'tight')
    # plt.show()
    
    