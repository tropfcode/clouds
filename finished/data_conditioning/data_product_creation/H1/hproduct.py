"""
PROGRAM DESCRIPTION:
   Assign every gridbox of HGG data to a WS in Equal-Angle format.

PROGRAM IMPLEMENTATION:
   Supply as command line arguments the start year for gridboxes
desired to be assigned, the path to the weather states (WS), and
name for outputfile of gridbox assignments.
   Each process of the program will work on a single year of data
HGG data at a time, with rank 0 working on the command line
argument start_year. 
   Every process will load the provided WS and then load a single
3-hour HGG file at a time in temporal order -- starting with Hour
0000 Jan 01 and ending on 2100 Dec 31 of the same year. From the
HGG data, relevant variables are loaded and then given as arguments
to Fortran subroutines to assign every gridbox to a WS and also
convert the data to Equal-Angle (see hproduct_analysis.f90 for details).
   The resulting assignments are saved as a 4D numpy array of shape
(12,248,360,180). This can be thought of as a 3D array of shape
(248,360,180) for each month with 248 maps of shape (360,180),
one for each of the 3-hour periods in the month.
   Results from program are N 4D numpy arrays where N are the number
of years determined by user.

NOTE: Every month is represented with the same shape array (248,360,180).
The 248 comes from 8*31=248, which is 8 files for each day of the month.
Not every month has 31 days, so this means all gridboxes beyond the
max number of files for that month will be populated with 0.
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
final_array = np.zeros((12,248,360,180), dtype=np.int)

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
      final_array[i,j,:,:] = wsnums[:,:]
      count += 1
   count_list.append(count)
np.save('{}_{}'.format(year,output), final_array)
print('rank:{} countlist:{}'.format(rank, count_list))
