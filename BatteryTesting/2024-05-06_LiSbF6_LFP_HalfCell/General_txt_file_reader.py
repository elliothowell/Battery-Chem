# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 11:43:20 2024

@author: Danny
"""

import os
import re

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

"""
Function to read general txt files it will automatically skip header and read only the actual data
"""

def read_file(files): # NEed to define a file to read first
    
    dataRaw=[]
    dataTemp=[]
    tempHeader=[]
    temp=[]
    
    for i,j in enumerate(files):      
        with open(j, 'r', encoding='utf-8', errors='replace') as f:
            for j in f.readlines():
                dataRaw.append(j)
  
        for m,n in enumerate(dataRaw):
            if re.fullmatch("^[0-9eEX\-\t\.\+\s\,]+\n", n):
                entry_i= m
                header_i= m-1
                break 
        data_trace= np.genfromtxt(dataRaw, skip_header= entry_i)
        dataTemp.append(data_trace)
   
    #Optional code to store the header name to identify which column of data you want to read
        if header_i>-1:
           temp=dataRaw[header_i].split("\t")
           if "\n" in temp:
               temp.remove("\n")
               
           tempHeader.append(temp)
              
        elif header_i==-1:
            len_header= len(dataRaw[entry_i].split("\t"))
            for n in range(len_header):
                temp.append(f"col {n+1}")
            tempHeader.append(temp)
            temp=[]
    
    dataTemp = dataTemp[0]
    tempHeader = tempHeader[0]
    
    return dataTemp, tempHeader, dataRaw

#files=[r"/Users/elliothowell/SynologyDrive/Research - Elliot Howell/UQAM - Battery/Half Cells/2024-05-06_LiSbF6_LFP_HalfCells/353EH1/353EH1_FullCycling_01_GCPL_C02.txt"]
#data, header, raw=read_file(files)

