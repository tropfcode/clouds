"""
PROGRAM DESCRIPTION:
   For a description for this please see the description of
hproduct.py for the creation of the H1 Equal Angle dataset.
The only difference here is that the final output data is
in Equal Area format. This is seen in the final_array variable.
"""

import gc
import os
import sys
import numpy as np
from mpi4py import MPI
import analysispy as ap
from netCDF4 import Dataset

# MPI Variables
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Input arguments
start_year = int(sys.argv[1])
wspath = sys.argv[2]
output = sys.argv[3]

# Files must have these as hourly time stamps each file
hourlist = ['0000', '0300', '0600', '0900', '1200', '1500', '1800', '2100']
timelist = ['{}.{}.GPC'.format(i,h) for i in range(1,32) for h in hourlist]

# Assign year based on MPI rank
year = start_year + rank
print(year)

# Obtain file paths
base_path = '/discover/nobackup/dtropf/hgg/{}/'.format(year)
month_paths = [base_path+path+'/' for path in sorted(os.listdir(base_path))]

# Testing
count_list = []
if rank == 0:
   print(len(month_paths))
   print(month_paths)
   print('\n')


# Load weather state data
with Dataset(wspath) as data:
   centroids = data['centroids'][:,:]
weather_states = np.zeros((centroids.shape[0],6,7))
for i in range(centroids.shape[0]):
   weather_states[i,:,:] = centroids[i,:].reshape(6,7)

# Create empty array to hold results
final_array = np.zeros((12,248,41252), dtype=np.int)

# Create h1 dataset
for i, mpath in enumerate(month_paths):
   print('MONTH: {}'.format(mpath))
   count = 0 # Used to keep track of files in month directory
   for j, dfile in enumerate(sorted(os.listdir(mpath))):
      dpath = mpath+dfile

      # Check if files are being read in correct order
      # and that files are proper type
      if timelist[j] not in dpath[-31:-20] or '.nc' not in dpath[-3:]:
         if rank == 0:
            print('BAD FILE: {} {} {}'.format(i,j,dpath))
            print('\n')
         break

      # Print path of each file for root process
      if rank == 0:
         print(dpath)
         print('\n')

      # Load data
      with Dataset(dpath) as data:
         pctaudist = data['n_pctaudist'][:,:,:]
         ntotal = data['n_total'][:]
         sqlonbeg = data['sqlon_beg'][:]
         sqlonend = data['sqlon_end'][:]
         eqlat = data['eqlat_index'][:]

      # Normalize data, assign each gridbox to a WS, and transform to Equal-Angle
      pctaudist_norm = ap.fnorm(pctaudist, ntotal,pctaudist.mask)
      wsnums = ap.findws(pctaudist_norm, weather_states[:,:,:], pctaudist.mask,
                         sqlonbeg, sqlonend, eqlat, centroids.shape[0])
      final_array[i,j,:] = wsnums[:]
      count += 1
   count_list.append(count)
np.save('{}_{}'.format(year,output), final_array)
print('rank countlist', rank, count_list)
