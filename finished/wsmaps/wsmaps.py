"""
PROGRAM DESCRIPTION:
   Create Equal Angle Weather State (WS) distribution maps for each
WS over some time period.

PROGRAM IMPLEMENTATION:
   This program is meant to run on the NCCS Discover computing system using
SLURM. To run this program, use it's associated submit.sh file.
   For each year that maps are to be made, make a separate parallel process
using MPI. Each process then loads the same centroids used for WS determination.
   Given k centroids, an empty matrix of shape (12,k+1,360*180) is made for each
process. The first dimension is number of months in year, second is number of
WS types plus no clouds (hence k+1), and finally the number of gridboxes in each
map. So, this matrix can be thought of as a stack where each layer is represents a
3-hour period for a given month hold k+1 flattened maps. Matricies for
Sum Squared Error (SSE) and total counts are intialized as well.
   For each 3-hour period availble to the respective processes, the n_pctaudist,
sqlon_beg, sqlon_end, eqlat_index, and n_total varibles are loaded into memory.
Using methods found in the analysis.f90 file (converted to python using f2py),
the PC-TAU data is normalized. Following normalization, another analysis.f90
method is used to determine the WS designation of each PC-TAU gridbox
by calculating the minimal euclidean distance between every centroid.
During this process, the maps are also converted from Equal Area to Equal Angle.
   After all gridboxes have been designated a WS, each process saves its own
netcdf file with the WS maps, total counts, SSE, and the centroids used for
designating each gridbox.
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
base_path = '/discover/nobackup/dtropf/hgg/'
data_path = sys.argv[3]
save_id = str(sys.argv[4])
save_path = os.getcwd()+'/'
years = [str(i) for i in range(year_start, year_end+1)]
year_path = base_path+years[rank]
print(data_path)

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
# If Gridbox has no data, skip
# If Gridbox has no clouds, do not calculate euclidean distance
# Else calculate euclidean distance against every centroid and find smallest distance
for i, f1 in enumerate(sorted(os.listdir(year_path))):
    month_path = year_path+'/'+f1
    for f2 in sorted(os.listdir(month_path)):
        data_path = month_path+'/'+f2
        if '.nc' not in data_path:
            print('SOMETHING WRONG WITH RANK {} FILE {}'.format(rank,data_path))
            continue
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
