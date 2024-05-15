# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 11:05:17 2024

Reader Functions

@author: Elliot
"""
import os 
import h5py
import random
import time
import re
import numpy as np

#If an HDF5 file hasn't been created, then call this function!

def HDF5_creator(dir,fn):
    f = h5py.File(dir+'\\'+fn, "w")
    print('Your HDF5 was successfully created in '+dir+'\\'+fn+'.')
    
###############################################################################
'''
HEKA ASC notebook file import.

This takes everything you just created in the cell above and adds it to the 
master HDF5 file you specify.

Make sure the file is opened in 'a' mode and not 'w' mode. 'w' mode deletes all
other HDF5 entries.
'''

def Add2HDF5(dir, fileName, data, notebook, notebook_header, metaD):

    #Putting something here that checks if the file exists would be nice.
    
    #This will get us the column we want from the numpy array called 'data' 
    #created earlier.
    def getColumnData(seriesNumber,sweepNumber,columnNumber):
        newList=[]
        for iiii in range(len(data[seriesNumber][sweepNumber])):
            newList.append(data[seriesNumber][sweepNumber][iiii][columnNumber])
        return newList
    
    with h5py.File(dir,'a') as f: #Make sure this is an 'a' and not 'w'. Only use 'w' once!!!
        
        #Build a framework to eventually put variables in.
        #If the filename doesn't exist, create a new entry
        #If it already exists, 
        if f.__contains__(fileName): #make sure the filename hasn't already been created
            metaData = f[fileName]
            dataDirectory= f[fileName]['Data']
            notebookDirectory = f[fileName]['Notebook']
        else:
            uniqueID=random.randint(0,1000000) #This is currently not used but could be used for your project!
            metaData = f.create_group(fileName)
            dataDirectory= metaData.create_group('Data')
            notebookDirectory = metaData.create_group('Notebook')

        
        #Add the data into the HDF5!
        for i in range(len(data)): #Series
            #Check to see if the data has alreadyb been added. We don't need to add it again.
            if f.__contains__(fileName+'/Data/Serie '+str(i)):
                print('This has already been added and was not overwritten.')
            else:
                series=dataDirectory.create_group('Serie '+str(i))
                
                for j in range(len(data[i])): #Sweep
                    # print(j) #good for troubleshooting
                    for k in range(len(data[i][j][0])): #Columns
                        series.create_dataset('Sweep'+str(j)+'Column'+str(k),
                                              data=getColumnData(i,j,k))

        #Add the notebook into the HDF5
        for i in range(len(notebook)):
            if f.__contains__(fileName+'/Notebook/Serie '+str(i)):
                myFavoriteColor='The most intense color for humans.'
            else:
                notebookSeries=notebookDirectory.create_group('Serie '+str(i)) 

                for j in range(len(notebook[i])): #Sweep
                    # print(j) #good for troubleshooting
                    if notebook:
                        for l in range(len(notebook[i][j])): #Notebook placements, I suppose just like columns.
                            notebookSeries.attrs[notebook_header[0][l].strip()+str(j)] = notebook[i][j][l] 

        #Add the metadata into the HDF5!
        for i in range(len(metaD)):
            metaData.attrs[metaD[i][0]] = metaD[i][1]

        #Check to see if the file was successfully added.
        flag=0
        for i in f.keys():
            print(i)
            if i == fileName:
                print(fileName+' was successfully added to the HDF5 file specified.')
                flag=1
        if flag == 0:
            print(fileName+' was not added to the HDF5 file specified.')
    
###############################################################################
'''
This creates 'data', and 'notebook' which contains the data from your file and 
the metadata about individual scans. 'data' and 'notebook' have corresponding 
positions, so scan 5's data can be found in data[4] and scan 5's metadata in 
notebook[4] 'data' is currently a list of lists of lists. 

The top level is each series (64 long for test.asc), the next level is each 
sweep (4 long for series 64 in test.asc, usually only 1 long), and finally the 
data is stored in relatively the same format you would see if you opened 
test.asc in Excel to look at a sweep.

If 'data' is consistently stored in this format, the function used after this 
one will work properly no matter what is stored in it.
'''

#As for specifically whats going in for HEKA file reading, please refer to Danny!!!
def Reader_HEKA_JonFormat(dir,fn, directory):
    os.chdir(dir)
    st = time.time()
    #Separate the entire file into an array that has comma separated values.
    
    #This is Danny's code for reading the ASC files coming out of the HEKA.
    with open(directory+fn, 'r') as f:
        data_raw=[]
        for i in f.readlines():
            data_raw.append(i.split(','))
     
    """
    The data is comprised of each line of the HEKA.asc files as a list of items 
    separated by a comma.
    
    However, the data is not organized in a way that is convenient for us
    
    Also the values of the measurement are in a weird string format that is not 
    appropriate for data treatment.
    
    The data in the HEKA.asc is organized this way each 
        - Series_X_X: Each series correspond to a bunch of measurement and for 
        mapping it is one coordinate of the map
            - Sweeps from measurement: Each sweeps contains the data from the 
            measurement and is usually comprised of 5 columns:
                (index, time, I, time, E)
            you can have multiple sweeps per series for example OCP and the 
            then LSV
            - Notebook information: Additionnal data from HEKA that give us 
                                    information from the experiment
            There is only one entry of the notebook per serie and it is always 
            the last element of the serie
            
    To organize the data we will read trhough the data variable line per line
    and based on the values read of the first item of the line we can determine 
    what kind of information is stored in the line
    
    We are reading the first item of the data line so it will be data[0] or
    n[0] in this case
   
    We identify using the command re.search which looks for a cluster of 
    character in a list of string
    
        - If it's a serie -> you can identify using "Serie"
        - If it's a sweep -> you can identify using "Sweep_"
        - If it's a notebook entry -> you can identify using "sweep #"
        - A data entry for a measurement -> you can identify using " "
    
    However both the sweep and notebook data entry have " " as an identifier so 
    we need an additional identifier
    
    To do we add a tracker 'sweep_tracker' and 'meta_tracker'
    
    They start at 0 to indicate that they are inactive
    
    If we identify a sweep it will also trigger the change to make 
    sweep_tracker active and notebook_tracker inactive
    
    Conversely, if we identify a notebook it will also trigger the change to 
    make sweep_tracker inactive and notebook_tracker active
    
    Based on these identifier and tracker we can make an algorithm that will 
    place the data in a list of list that have the data conveniently placed in 
    the following order
    
    Level 1 -> Master list: each element is a serie
     level 2 -> Serie list: each element is a sweep
         level 3 -> Sweep list: each element the data form a measurement
             level 4 -> Column: each element consist of the type of the 
             attribute of the measurement (potential, time, current, ...) 
    """
    data=[] #Permanent storage for organized data from the HEKA files to be eventually exported
    data_header=[]
    data_serie=[] #Temporary container for storing "series_data" that will be appended to the "master_data" and reset after a serie is completed
    data_sweep=[]  #Temporary container for storing "sweep_data" that will be appended to the "seies_data" and reset after a serie is completed
    data_line=[]
    
    serie=[]
    test=[]
    notebook=[]
    notebook_line=[]
    notebook_header=[] #Permanent storage of notebook header to identify notebook data
    notebook_sweep=[]
    
    sweep_tracker=0 #Tracker for storing a data entry in the "sweep_data"
    notebook_tracker=0 #Tracker for storing a data entry in the "notebook_data"
    
    for n in data_raw: #Going through the HEKA file stored in data line per line
    
    #Identifying a sweep header means two thing we just started a new serie or we ended a sweep and start a new sweep         
        if re.search('Sweep_', n[0]): 
            sweep_tracker=1 #If we identify a "Sweep" then we activate the "sweep_tracker"
            notebook_tracker=0
            if len(data_sweep)>0: #If we just started the serie then sweep_data is empty then len(sweep_data=0)
                data_serie.append(data_sweep) #If len(sweep_data)>0 it means we ended a sweep and starting a new one so we need to append the sweep data to the serie
                data_sweep=[] #Additionnaly we need to clean up the sweep_data for the next sweep 
        if re.search('Sweep #', n[0]):
            if len(notebook_header)<1:
                notebook_header.append(n)
     
    #Identifying a " " means a data entry from a sweep or the notebook 
        if re.search(' ', n[0][0]):
            if sweep_tracker==1: # if its a sweep then we append to 'sweep_data'
                for element in n:
                    data_line.append(float(element))
                
                data_sweep.append(data_line)
                data_line=[]
                                   
            if notebook_tracker==1: # if its a notebook then we append to 'notebook_data'
                for element in n[:-1]:
                    
                    if element[-1]=="m":
                        notebook_line.append(float(element[:-1])*1e-3)
                    if element[-1]=="µ":
                        notebook_line.append(float(element[:-2])*1e-6)
                    if element[-1]=="n":
                        notebook_line.append(float(element[:-1])*1e-9)
                    if re.search("[0-9]", element[-1]):
                        notebook_line.append(float(element))
                    if re.search("NAN", element):
                        notebook_line.append("14114")           
                
                if n[-1][-2]=="m":
                    notebook_line.append(float(n[-1][:-2])*1e-3)
                if n[-1][-2]=="µ":
                    notebook_line.append(float(n[-1][:-3])*1e-6)
                if n[-1][-2]=="n":
                    notebook_line.append(float(n[-1][:-2])*1e-9)  
                if re.search("[0-9]", n[-1][-2]):
                    notebook_line.append(float(n[-1][:-1]))
                if re.search("NAN", n[-1]):
                    notebook_line.append("14114")
                    
                notebook_sweep.append(notebook_line)
                notebook_line=[]
                
        if re.search('Serie', n[0]): 
            serie.append(n[0])
        
    # Identifying a line skip (\n) as the first element always mean end of a sweep and the start of a new sweep           
        if re.search('\n', n[0]):
            
            if len(data_sweep)>0:
                data_serie.append(data_sweep) #Since we are starting the notebook section that means we reached the end of a sweep then we need to append the 'sweep_data' to the 'serie_data'
                data_sweep=[] # Clear out for new sweep   
                sweep_tracker=0 #If we identify a "Sweep" then we activate the "sweep_tracker"
                notebook_tracker=1
                       
            if len(notebook_sweep)>0: # If it's the end of a sweep and the
                data.append(data_serie)
                data_serie=[] # Clear out for new serie
                notebook.append(notebook_sweep)
                notebook_sweep=[]              
    
    X_pos=[]
    Y_pos=[]
    header_i= {'x':'', 'y':'', 'z':'', 'scanrate':''}
    technique= []
    
    for i,j in enumerate(notebook_header[0]):
        if re.search("X-pos", j):
            header_i['x']= i
            
        if re.search("Y-pos", j):
            header_i['y']= i
            
        if re.search("Z-pos", j):
            header_i['z']= i
            
        if re.search("Rate", j):
            header_i['scanrate']= i
        
    for m,n in enumerate(notebook):
        
        X_pos.append(np.abs(np.round((n[0][header_i['x']]-notebook[0][0][header_i['x']])*1e6,0)))
        Y_pos.append(np.abs(np.round((n[0][header_i['y']]-notebook[0][0][header_i['y']])*1e6,0)))   
    
        if n[0][header_i['scanrate']]!= "14114":
            technique.append('Voltametry')
            
        if n[0][header_i['scanrate']]== "14114":
            if len(data[m][0][0])==5:
                technique.append('OCP')
                
            if len(data[m][0][0])==3:
                technique.append('tip move')
            
    OCP= []
    
    Tn= len(dict((i, Y_pos.count(i)) for i in technique))      
    
    for i,j in enumerate(data):
            
        if technique[i]== 'OCP':
            OCP.append(j[0][-1][-1])
            
    # plt.scatter(X_pos[0:-1:Tn], Y_pos[0:-1:Tn], s=250, c=OCP, cmap='viridis')
    # plt.colorbar()
    
    et = time.time()
    print('HEKA_JonFormat_DBR() took '+str(round(et-st,2))+' s to run.')

    return data, notebook, technique, notebook_header
            

        