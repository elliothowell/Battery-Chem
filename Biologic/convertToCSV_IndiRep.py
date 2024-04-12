# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 12:50:13 2024

@author: Elliot
"""

"""
This file is just for converting CV files to txt/csv files white have just two columns
First column == Voltage
Second column == current
"""
import os
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# # Folder path containing the .mpt files ** windows **
folder_path = r'C:\Users\Elliot\SynologyDrive\Research - Elliot Howell\Durbis CV Measurements\TMPA\27-02-2024'
# # Folder path containing the .mpt files ** mac **
#folder_path = r'/Users/elliothowell/SynologyDrive/Research - Elliot Howell/UQAM - Battery/LiSbF6 Salt/15-02-2024 Blank ECDMC/Bottle A'

os.chdir(folder_path)

#Creating new directory for CSV files, if theres already a CSV_files folder we don't make one
new_path = folder_path +  r'\CSV_Files'
if not os.path.exists(new_path):
    os.makedirs(new_path)

# Get a list of files ending with '.mpt' in the folder ** WINDOWS **
txt_files = [f for f in os.listdir(folder_path) if (f.endswith('.mpt') or f.endswith('.txt'))]

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
        csv_name = file_name.split('.txt')[0] + '_CSV_Rep' + str(i+1)
        repList[i][['Ewe/V', '<I>/mA']].to_csv(new_path + '/' + csv_name + '.csv', header = False, index = False)
        

    
    
    
    