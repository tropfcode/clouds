"""
PROGRAM DESCRIPTION:
   Count the number of gridboxes that have clouds, no clouds, or
no data from every 3-hour HGG datafile for a single year. 
Output the results in time resolution of one month.

PROGRAM IMPLEMENTATION:
   Have a single CPU cycle through every 3-hour file of a single month
from a single year.
   For each 3-hour file, obtain the 'n_pctaudist' variable which is a 
3D numpy array of shape (6,7,42252). Cycle through the 3rd dimension
of size 42252 to determine if a given gridbox (2D array of shape (6,7))
is one of three things:
   1. No Data: All values in the 2D array are masked
   2. Cloud Data: No masked values and at least one non-zero value
   3. No Cloud Data: No masked values and all zero

   Utilize MPI to do 12 months in parallel. Save results as 2D numpy array
of shape (n_files, 3).
"""

import os
import gc
import sys
import time
import numpy as np
from mpi4py import MPI
from netCDF4 import Dataset

# MPI VARIABLES
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Command line arguments
year = sys.argv[1]

# Get month worth of paths for each processor for given year
base_path = '/discover/nobackup/dtropf/{}/'.format(year)
data_path = base_path + os.listdir(base_path)[rank] + '/'
data_files = [data_path+f for f in os.listdir(data_path)]

# Initalize array to hold counts of bad gridboxes, no clouds, and data
counts_array = np.zeros((len(data_files), 3), dtype=np.int32)

t1 = time.time()
# Compute counts for every 3-hour file for given month
for i,dataf in enumerate(data_files):
   nocloud = 0
   with Dataset(dataf) as data:
      pctaudist = data['n_pctaudist'][:,:,:]
      pctaudist1d = np.sum(np.sum(pctaudist, axis=0), axis=0)
      cloud_count = np.nonzero(pctaudist1d)[0].shape[0]
      bad_count = np.where(pctaudist1d.mask==True)[0].shape[0]
      for j in range(pctaudist.shape[2]):
         if np.all(pctaudist[:,:,j]==0):
            nocloud += 1
      #nocloud = np.where(pctaudist1d==0)[0].shape[0]
      #cloud_count = np.unique(np.nonzero(pctaudist[:,:,:])[2]).shape[0]
      #bad_count = np.unique(np.where(pctaudist.mask==True)[2]).shape[0]
      counts_array[i,0] = int(cloud_count)
      counts_array[i,1] = int(bad_count)
      counts_array[i,2] = int(nocloud)
   del pctaudist,pctaudist1d,cloud_count,bad_count,nocloud
   gc.collect()
t2 = time.time()
print('TIME AND RANK: {} {}'.format(t2-t1,rank))
np.save('counts_{}_{}'.format(year, rank), counts_array)
