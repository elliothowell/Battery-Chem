# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 11:35:35 2024

HDF5 Files/creation and such

Opening/creating file modes:
r:          Readonly, file must exist (default)
r+:         Read/write, file must exist
w:          Create file, truncate if exists
w- or x:    Create file, fail if exists
a:          Read/write if exists, create otherwise

@author: Elliot
"""

import h5py
import numpy as np
import os


# Line by line:
    # this process keeps hdf5 file open though so to delete will need to 
    # manually close...
os.chdir(r'C:\Users\Elliot\Documents\Python Scripts\myDataBase')
f = h5py.File('mydataset.hdf5', 'w')
dset = f.create_dataset('mydataset', (100,), dtype = 'i')
# f.close()
# to close file simply use: f.close()

# displaying how to append group:
f = h5py.File('mydataset.hdf5', 'a')
grp = f.create_group('subgroup')

dset2 = grp.create_dataset('anotha_one', (50,), dtype = 'f')
dset2.name

dset3 = f.create_dataset('subgroup2/dataset_three', (10,), dtype = 'i')
dset3.name

dataset_three = f['subgroup2/dataset_three']

# Or using with statement (this way the file is closed after this process):
# with h5py.File('testFile.hdf5', 'w') as f:
#     dset = f.create_dataset('mydataset', (100,), dtype = 'i')
    
    