"""
PROGRAM DESCRIPTION:
  Create Weather State Maps (wsmaps) from centroids derived from k-means clustering.
Weather state maps are the geographic distribution of a weather state (centroid).

PROGRAM IMPLEMENTATION:
   Using MPI, y processes are spawned where y is the number of years of HGG data chosen
for the creation of wsmaps. Choosing k centroids derived from k-means clustering,
create an empty 3D matrix of shape (t,k,360*180) where t is the number of months
chosen to analyze in a year. The third dimension is of static size 
due to the conversion of HGG data from Equal-Area to Equal-Angle.
   Each 2D slice of the 3D matrix of shape (k,360*180) can be thought of as k maps 
flattened to a vector of length 360*180. Each cell of the vector corresponds directly 
to a gridbox in the raw HGG data. In the same vein, each vector in the
2D array is dedicated to one of the k centroids. The cells in each vector counts how many
times the gridbox in that geographic position from raw HGG data resulted in a minimal
euclidean distance amongst all the k centroids.
   For example, if the euclidean distance between the first gridbox of the first three hour file 
for year 1, month 1, of HGG data and the kth centroid was minimal, then cell (0,k-1,0) in python or
(1,k,1) in Fortran of process 1 would be incremented by 1. No other cells would be augmented from the analysis
of that specific gridbox of the raw HGG data.
   This process happens for all gridboxes across all the selected years by the user.
All years are processed in parallel.
   For each process, the number of times a centroid resulted in a minimal euclidean distance, counts,
and the sum of squared errors, sse, are stored. These are 1D vectors of length k.
   The calculations of euclidean distance, sse, and counts are done using fortran subroutines
found in analysis.f90. The subroutines can be called from python code by first compiling the fortran
code using the f2py program.
   The output of this program are y netCDF files with each file containing the following:
      1. A 3D matrix of shape (t,k,360*180) which are the wsmaps for each month in a specific year
      2. A vector of length k called counts which are the total gridbox membership for each centroid
across the total year (this can also be calculated by summing across the first and third dimension of the
3D wsmap matrix)
      3. A vector of length k called sse which is the summed squared error associated with each centroid.

"""

import os
import gc
import sys
import time
import analysispy
import numpy as np
from mpi4py import MPI
from netCDF4 import Dataset

# Start clock to time program
t1 = time.time()

# Initialize MPI variables
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# Initalize other variables
year_start = int(sys.argv[1])
year_end = int(sys.argv[2])
base_path = '/discover/nobackup/dtropf/'
data_path = sys.argv[3]
save_id = str(sys.argv[4])
save_path = os.getcwd()+'/'
years = [str(i) for i in range(year_start, year_end+1)]
year_path = base_path+years[rank]

# Load centroids and reshape for analysis
data = Dataset(data_path)
pctau_histograms = data['centroids'][:,:]
data.close()
centroids = np.zeros((pctau_histograms.shape[0],6,7))
for i in range(pctau_histograms.shape[0]):
    centroids[i,:,:] = pctau_histograms[i,:].reshape(6,7)

# Initialize arrays to hold counts of each weather state and weather state map data
k = int(centroids.shape[0])
print(centroids.shape, np.sum(centroids))
all_data = np.zeros((12,k+1,360*180))
totalcounts = np.zeros(k+1)
sse = np.zeros((k))

# Loop over all data in the year according to rank of process
# First normalize the gridbox data
# If Gridbox has no data, skip
# If Gridbox has no clouds, do not calculate euclidean distance
# Else calculate euclidean distance against every centroid and find smallest distance
for i, f1 in enumerate(sorted(os.listdir(year_path))):
    month_path = year_path+'/'+f1
    for f2 in sorted(os.listdir(month_path)):
        data_path = month_path+'/'+f2
        with Dataset(data_path) as data:
            pctaudist = data['n_pctaudist'][:,:,:]
            sqlonbeg = data['sqlon_beg'][:]
            sqlonend = data['sqlon_end'][:]
            eqlat = data['eqlat_index'][:]
            ntotal = data['n_total'][:]
            pctaudist_norm = analysispy.fnorm(pctaudist, ntotal, pctaudist.mask)
            ws_num , counts, tmp_sse = analysispy.findws(pctaudist_norm, centroids, pctaudist.mask, sqlonbeg, sqlonend, eqlat, k)
            for j in range(k+1):
                all_data[i,j,:] = all_data[i,j,:] + ws_num[:,:,j].reshape(360*180)
            totalcounts[:] = totalcounts[:] + counts[:]
            print(type(tmp_sse), tmp_sse)
            sse += tmp_sse
            del pctaudist, sqlonbeg, eqlat, ntotal
            gc.collect()

t2 = time.time()
print('TOTAL TIME FOR '+str(rank))

# Create netCDF4 file
dataset = Dataset(save_path+'/{}_{}.nc'.format(save_id, str(years[rank])), 'w', format='NETCDF4_CLASSIC')
dataset.description = 'year start: {} \nyear end: {}'.format(str(year_start), str(year_end))

# Add dimensions for variables
dim = dataset.createDimension('dim', pctau_histograms.shape[1])
kdim = dataset.createDimension('kdim', pctau_histograms.shape[0])
timedim = dataset.createDimension('timedim', all_data.shape[0])
wsdim = dataset.createDimension('wsdim', all_data.shape[1])
mapdim = dataset.createDimension('mapdim', all_data.shape[2])
ssedim = dataset.createDimension('ssedim', k)

# Save centroids, wsmaps, and counts
var_initcentroids = dataset.createVariable('init_centroids', np.float32, ('kdim', 'dim'))
var_wsmaps = dataset.createVariable('wsmaps', np.float32, ('timedim', 'wsdim', 'mapdim'))
var_counts = dataset.createVariable('counts', np.int32, ('wsdim'))
var_sse = dataset.createVariable('sse', np.float32, ('ssedim'))
var_initcentroids[:,:] = pctau_histograms[:,:]
var_wsmaps[:,:,:] = all_data[:,:,:]
var_counts[:] = totalcounts[:]
var_sse[:] = np.asarray(sse)[:]
dataset.close()
for i in range(12):
    print(i, totalcounts[i]/np.sum(totalcounts))
