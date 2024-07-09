# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 16:44:15 2024

@author: Elliot
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.signal import find_peaks

def extract_last_cycle(data):
    potential = data['Potential applied (V)']
    start_potential = potential.iloc[0]
    threshold = 0.5 * (potential.max() - potential.min())  # 50% of the range as threshold
    peaks, _ = find_peaks(potential)
    troughs, _ = find_peaks(-potential)
    turning_points = np.sort(np.concatenate((peaks, troughs)))
    cycle_starts = [i for i in turning_points if abs(potential.iloc[i] - start_potential) < threshold]
    if len(cycle_starts) > 1:
        start = cycle_starts[-2]
        end = cycle_starts[-1]
        return data.iloc[start:end+1]
    else:
        return None

def plot_and_save_last_cycle(folder_path):
    last_cycle_folder = os.path.join(folder_path, 'lastCycles')
    if not os.path.exists(last_cycle_folder):
        os.makedirs(last_cycle_folder)
    
    csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]
    
    for file in csv_files:
        file_path = os.path.join(folder_path, file)
        data = pd.read_csv(file_path)
        last_cycle_data = extract_last_cycle(data)
        
        if last_cycle_data is not None:
            last_cycle_filename = f'lastCycle_{file}'
            last_cycle_data.to_csv(os.path.join(last_cycle_folder, last_cycle_filename), index=False)
            
            plt.figure()
            plt.plot(last_cycle_data['Potential applied (V)'], last_cycle_data['WE.Current (A)'], label='Compound')
            
            parts = file.split('_')
            parts[1] = 'post' + parts[1]
            blank_file_name_pattern = '_'.join(parts[:-1]) + '_'
            blank_files = [f for f in csv_files if f.startswith(blank_file_name_pattern)]
            
            if blank_files:
                blank_data = pd.read_csv(os.path.join(folder_path, blank_files[0]))
                last_blank_cycle = extract_last_cycle(blank_data)
                if last_blank_cycle is not None:
                    plt.plot(last_blank_cycle['Potential applied (V)'], last_blank_cycle['WE.Current (A)'], label='Blank', linestyle='--')
            
            plt.xlabel('Potential applied (V)')
            plt.ylabel('WE.Current (A)')
            plt.title(f'Last Cycle Plot for {file}')
            plt.legend()
            plt.grid(True)
            #plt.savefig(os.path.join(folder_path, f'Last_Cycle_{file}.png'))
            plt.show()
        else:
            print(f"Not enough cycle starts found in {file}, skipping.")

# Usage
folder_path = r'C:/Users/Elliot/SynologyDrive/Research - Elliot Howell/Durbis CV Measurements/Metrohm/2024-07-04'
plot_and_save_last_cycle(folder_path)