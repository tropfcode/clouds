"""
PROGRAM DESCRIPTION:
   The purpose of this program is to filter out all PC-TAU HGG data that contains clouds.
Saved data is of shape (samples,42) where samples is the number of gridboxes with 
clouds and 42 is number of PC-TAU histogram bins. The saved data is in chronological order.

PROGRAM IMPLEMENTATION:
   Provide as command line input the start and end years of HGG data of interest.
Using MPI, the program runs processes in parallel for each year.
   Open a new netcdf file createing a single variable of shape (samples,42).
Next, for every 3-hour period in a year, obtain the PC-TAU and ntotal data.
Divide PC-TAU data by ntotal to normalize, then filter out the PC-TAU
histograms with cloud data and save to the netcdf file.
   Each process saves it's own file.
"""

import os
import gc
import sys
import time
import numpy as np
from mpi4py import MPI
from netCDF4 import Dataset

# MPI variables
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Command line variables
year_start = int(sys.argv[1])
year_end = int(sys.argv[2])

# Obtain year for each process
years = [str(i) for i in range(year_start, year_end+1)]
year = years[rank]

# Create netcdf file to save flattned PC-TAU histograms with clouds
basepath = '/discover/nobackup/dtropf/hgg/'
dataset = Dataset(year+'.nc','w')
n_pctau = dataset.createDimension('samples', None)
pctau_flat = dataset.createDimension('features', 42)
var = dataset.createVariable('data', np.float32,('samples','features'))

# Go through every 3-hour HGG file for a year
# Obtain PC-TAU data, normalize, and filter out data containing clouds
t1 = time.time()
cycle = 0
count = 0
for path in sorted(os.listdir('{}{}/'.format(basepath,year))):
   for i, fpath in enumerate(sorted(os.listdir('{}{}/{}'.format(basepath,year,path)))):
      path2 = '{}{}/{}/{}'.format(basepath,year,path,fpath)
      if '.nc' not in path2:
         continue
      with Dataset(path2) as data:
         pctaudist = data['n_pctaudist'][:,:,:]
         ntotal = data['n_total'][:]
         normed = pctaudist/(ntotal*1.0)
         # Replace masked values with zeros
         rmv = np.where(pctaudist.mask,0,normed)
         # Filter out all histograms containing at least one non-zero value (has clouds)
         rmv2 = rmv[:,:,np.unique(np.nonzero(rmv)[2])]
         # Properly reshape, flatten, and save to file
         var[count:count+rmv2.shape[2],:] = np.reshape(np.transpose(rmv2,(2,0,1)),(rmv2.shape[2],42))
         #var[count:count+rmv2.shape[2],:] = rmv2.flatten().reshape((rmv2.shape[2],42))[:,:]
         count += rmv2.shape[2]
         del pctaudist, ntotal, normed, rmv, rmv2
         gc.collect()
      cycle +=1
dataset.close()
t2 = time.time()
print('DONE Total Time FOR RANK', rank, t2-t1)
