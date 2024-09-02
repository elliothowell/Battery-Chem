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

def plot_battery(batteryDF):
    
