# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 10:52:16 2024

@author: Elliot
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.signal import savgol_filter

def apply_savgol_filter(data, window_length=2, polyorder=1):
    """ Apply Savitzky-Golay filter to smooth the data. """
    if window_length > len(data):
        window_length = len(data) // 2 * 2 + 1  # Make it the nearest odd number
    return savgol_filter(data, window_length, polyorder)

def extract_cycles(data):
    potential = data['Potential applied (V)']
    
    direction = potential.diff()
    cyclicStarts = []
    cyclicEnds = []
    cycleNum = 1
    # print(range(1, len(potential)))
    for i in range(1, len(potential)-1):
        # print(i)
        if cycleNum == 1:
            cyclicStarts.append(0)
            cycleNum += 1
        if potential[i] == 0.0 and direction[i-1] < 0 and direction[i+1] > 0:
            cyclicStarts.append(i)
            cyclicEnds.append(i)
            cycleNum += 1
        elif potential[i] == 0.0 and direction[i-1] > 0 and direction[i+1] > 0:
            cyclicStarts.append(i)
            cyclicEnds.append(i)
            cycleNum += 1
        elif i == len(potential)-2:
            cyclicEnds.append(i+1)
        else:
            continue

    if len(cyclicStarts) < 2:
        cyclicStarts = []
        cyclicEnds = []
        cycleNum = 1
        for i in range(1, len(potential) - 1):
            if cycleNum == 1:
                cyclicStarts.append(0)
                cycleNum += 1
            if potential[i] == np.min(potential):
                cyclicStarts.append(i)
                cyclicEnds.append(i)
            elif i == len(potential)-2:
                cyclicEnds.append(i+1)
            else:
                continue
            
    cycles = []
    for i in range(0, len(cyclicStarts)):
        cycles.append(data.iloc[cyclicStarts[i]:cyclicEnds[i]])
        #print(cyclicStarts[i])
        
    return cycles

def extract_last_cycle(data):
    cycles = extract_cycles(data)
    if cycles:
        return cycles[-1]
    return None

def plot_and_save_cycles(folder_path, plot_all_cycles, fc_potential, save_figures):
    # figures_folder = os.path.join(folder_path, 'figures')
    # if save_figures and not os.path.exists(figures_folder):
    #     os.makedirs(figures_folder)
    
    last_cycle_folder = os.path.join(folder_path, 'lastCycles')
    if not os.path.exists(last_cycle_folder):
        os.makedirs(last_cycle_folder)
    
    fc_folder = os.path.join(folder_path, 'fc_figures')
    if not os.path.exists(fc_folder):
        os.makedirs(fc_folder)
    
    csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]
    
    for file in csv_files:
        parts = file.split('_')
        parts[1] = 'post' + parts[1]
        blank_file_name_pattern = '_'.join(parts[:4]) + '_'
        blank_files = [f for f in csv_files if f.startswith(blank_file_name_pattern)]
        
        file_path = os.path.join(folder_path, file)
        data = pd.read_csv(file_path)
        last_cycle = extract_last_cycle(data)
        
        if last_cycle is not None:
            # Save the last cycle to a CSV file
            last_cycle_file_path = os.path.join(last_cycle_folder, f'last_cycle_{file}')
            last_cycle.to_csv(last_cycle_file_path, index=False)
        
        if plot_all_cycles:
            cycles = extract_cycles(data)
        else:
            cycles = [last_cycle] if last_cycle is not None else []
        
        plt.figure(figsize=(10, 6))
        num = 1
        for cycle in cycles:
            smoothed_current = (1*10**6) * apply_savgol_filter(cycle['WE.Current (A)'])
            potential_vsFc = cycle['Potential applied (V)'] - fc_potential
            if plot_all_cycles:    
                plt.plot(potential_vsFc, smoothed_current, label='Cycle ' + str(num), alpha=0.5)
            else:
                plt.plot(potential_vsFc, smoothed_current, alpha=0.5)
            num += 1
        
        # if blank_files:
        #     blank_data = pd.read_csv(os.path.join(folder_path, blank_files[0]))
        #     last_blank_cycle = extract_last_cycle(blank_data)
        #     if last_blank_cycle is not None:
        #         smoothed_current_blank = (1*10**6) * apply_savgol_filter(last_blank_cycle['WE.Current (A)'])
        #         plt.plot(last_blank_cycle['Potential applied (V)'], smoothed_current_blank, label='Blank Last Cycle', linestyle='--', color='red')
        
        # plt.xticks(np.arange(min(potential_vsFc), max(potential_vsFc) , step = 0.25))
        plt.xlabel('Potential applied (V vs. Fc/$Fc^+$)')
        plt.ylabel('Current ($\mu$A)')
        plt.title(f'{"All Cycles" if plot_all_cycles else "Last Cycle"} versus Ferrocene for {file}')
        if plot_all_cycles:
            plt.legend()
        if save_figures:
            plt.savefig(os.path.join(fc_folder, f'{"All_Cycles" if plot_all_cycles else "Last_Cycle"}_versus_Fc_{file}.png'))
        plt.show()

# Usage
folder_path = r'C:\Users\Elliot\SynologyDrive\Research - Elliot Howell\Durbis CV Measurements\0 - To use\DCM Solvent\14-PTH'
user_input_cycles = input("Plot all cycles? (y/n): ").strip().lower()
plot_all_cycles = True if user_input_cycles == 'y' else False
# user_input_blanks = input("Plot only if blank is associated? (y/n): ").strip().lower()
# plot_only_with_blanks = True if user_input_blanks == 'y' else False
user_input_save = input("Save figures? (y/n): ").strip().lower()
save_figures = True if user_input_save == 'y' else False
user_input_FcPot = float(input("Ferrocene Potential: "))
fc_potential = user_input_FcPot if user_input_FcPot != 0 else 0.00

plot_and_save_cycles(folder_path, plot_all_cycles, fc_potential, save_figures)