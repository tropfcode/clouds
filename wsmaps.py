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
centroid_path = '/discover/nobackup/dtropf/mpi4py/array_batch/centroids/'
save_path = '/discover/nobackup/dtropf/mpi4py/array_batch/wsmaps/'+sys.argv[3]
years = [str(i) for i in range(year_start, year_end)]
year_path = base_path+years[rank]
all_data = np.zeros((12,12,360*180))

# Load centroids and reshape for analysis
centroids = np.zeros((12,6,7))
pctau_histograms = np.load(centroid_path+sys.argv[4])
for i in range(pctau_histograms.shape[0]):
    centroids[i,:,:] = pctau_histograms[i,:].reshape(6,7)

# Loop over all data in the year according to rank of process
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
            ws_num = analysispy.findws(pctaudist_norm, centroids, pctaudist.mask, sqlonbeg, sqlonend, eqlat)
            for j in range(12):
                all_data[i,j,:] = all_data[i,j,:] + ws_num[:,:,j].reshape(360*180)

            del data, pctaudist, sqlonbeg, eqlat, ntotal
            gc.collect()
t2 = time.time()
print('TOTAL TIME FOR '+str(rank))
np.save(save_path+str(years[rank]), all_data)
