# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 15:29:06 2024

@author: Elliot
"""

import pandas as pd
import matplotlib.pyplot as plt

def plot_battery_data(file_path):
    # Load the data, skipping initial non-data lines
    data = pd.read_csv(file_path, delimiter = '\t', skiprows=77, encoding_errors='ignore')
    
    # Calculate cycle numbers
    # Assuming a new cycle starts with a new charge (mode 1 and ox/red 1) after a discharge
    conditions = (data['mode'] == 1) & (data['ox/red'].shift(1) != data['ox/red'])
    data['cycle number'] = conditions.cumsum()

    # Filter for charge and discharge cycles
    charging_data = data[(data['mode'] == 1) & (data['ox/red'] == 1)]
    discharging_data = data[(data['mode'] == 1) & (data['ox/red'] == 0)]
    
    # Plot 1: Time vs Ewe/V and <I>/mA
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(data['time/s'], data['Ewe/V'], 'g-')
    ax2.plot(data['time/s'], data['<I>/mA'], 'b-')
    
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Ewe (V)', color='g')
    ax2.set_ylabel('Current (mA)', color='b')
    plt.title('Time vs Ewe and Current')
    plt.show()

# Plot 2: Ewe/V vs Q charge/discharge without connecting cycles
    plt.figure()
    for cycle, group in charging_data.groupby('cycle number'):
        plt.plot(group['Q charge/mA.h'], group['Ewe/V'], label=f'Cycle {cycle}')
    for cycle, group in discharging_data.groupby('cycle number'):
        plt.plot(group['Q discharge/mA.h'], group['Ewe/V'], label=f'Cycle {cycle}')
    
    plt.xlabel('Capacity (mA.h)')
    plt.ylabel('Ewe (V)')
    plt.title('Ewe vs Capacity')
    plt.legend()
    plt.show()

    # Plot 3: Last Q charge/discharge value for each cycle
    last_charges = charging_data.groupby('cycle number')['Q charge/mA.h'].last()
    last_discharges = discharging_data, data.groupby('cycle number')['Q discharge/mA.h'].last()
    
    plt.figure()
    plt.plot(last_charges.index, last_charges, 'ro-', label='Last Charge Capacity')
    plt.plot(last_discharges.index, last_discharges, 'bo-', label='Last Discharge Capacity')
    plt.xlabel('Cycle Number')
    plt.ylabel('Capacity (mA.h)')
    plt.title('Last Capacity per Cycle')
    plt.legend()
    plt.show()

# Example usage
file_path = 'C:/Users/Elliot/SynologyDrive/Research - Elliot Howell/UQAM - Battery/Coin Cells/2024-05-06_LiSbF6_LFP_HalfCells/353EH1/353EH1_FullCycling_06_GCPL_C02.txt'
plot_battery_data(file_path)