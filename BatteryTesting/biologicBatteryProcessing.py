# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 16:13:42 2024

Going to try to make a script to just read in whatever file I have, C/20 etc. 
and export those 5 cycles

Each file will be its own C-rate, this will be correlated through my lab note-
book, all names should be self-explanatory

@author: Elliot
"""

import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import datetime as dt
import numpy as np
from matplotlib import cm
 
blues = cm.get_cmap("Blues", 5)
reds = cm.get_cmap("Reds", 5)

def bio_file_reader(bio_file):
    
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
    
    #adding column which indicates which cycle we are in
    read_file = read_file.assign(cycle = lambda x: x.Ns.apply(lambda y: y // 2 if y % 2 == 0 else (y + 1) // 2))
    
    #adding column which says if we are in "charge" or "discharge"
    conditions = [
    (read_file['mode'] == 1) & (read_file['ox/red'] == 1),
    (read_file['mode'] == 1) & (read_file['ox/red'] == 0),
    (read_file['mode'] == 3)]

    # Define choices
    choices = ['charge', 'discharge', 'rest']
    
    # Create new column using np.select
    read_file['state'] = np.select(conditions, choices, default='unknown')
    
    return read_file

def plot_battery(battery_file):
    
    #file_path = os.path.join(folderPath, file)
    
    #make a dataframe of the file
    fileDF = bio_file_reader(battery_file)
    
    #list for column names
    colNames = fileDF.columns
    print(colNames[25])
    ### [11] == potential, [19] == current, [21] == discharge cap, 
    ### [22] == charge cap, [24] == cycle, [25] == state
    
    #making charge/discharge DF
    chargeDF = fileDF[fileDF[colNames[25]] == 'charge']
    dischargeDF = fileDF[fileDF[colNames[25]] == 'discharge']
    
    
    #first making the figure
    figCap, axCap = plt.subplots(figsize = (10, 6), layout = 'constrained')
    
    for i in range(1, max(fileDF['cycle'])+1):
        charge = chargeDF.loc[chargeDF[colNames[24]] == i, colNames[22]]
        potCharge = chargeDF.loc[chargeDF[colNames[24]] == i, colNames[11]]
        discharge = dischargeDF.loc[dischargeDF[colNames[24]] == i, colNames[21]]
        potDischarge = dischargeDF.loc[dischargeDF[colNames[24]] == i, colNames[11]]
    
        axCap.plot(charge, potCharge, label = 'Charge Cycle ' + str(i), c = blues(i))
        axCap.plot(discharge, potDischarge, label = 'Discharge Cycle ' + str(i), c = reds(i))
        
    axCap.set_xlabel('Capacity (mAh)') 
    axCap.set_ylabel('Potential (V)')
    figCap.legend(loc = 'outside right upper')
    axCap.set_title('Charge/Discharge Capacities')
    
    
file = '/Users/elliothowell/SynologyDrive/Research - Elliot Howell/UQAM - Battery/Coin Cells/2024-05-06_LiSbF6_LFP_HalfCells/354EH2/354EH2_FullCycling_06_GCPL_C03.txt'

