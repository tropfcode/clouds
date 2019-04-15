"""
Old method for obtaining initialization centroids that was completely random
temporally and geographically.

This is replaced by init_centroids.py 
"""

import os
import sys
import time
import random
import numpy as np
from netCDF4 import Dataset

def init_centroids(k, base_path, year_start, year_end=None):
    """        
    """
    oldyears = []
    centroids = np.zeros((k, 42))
    if year_end is None:
        year_end = year_start + 1
    years = [i for i in range(year_start, year_end+1) if i != 2005]
    for i in range(k):
        t1 = time.time()
        print('OBTAINING CENTROID', i)
        if year_end-year_start <= 10:
            year = years[i]
        else:
            year = random.randint(year_start, year_end) 
            if len(oldyears) < 14:
                while year in oldyears:
                    year = random.randint(year_start, year_end)
        oldyears.append(year)
        print('YEAR', year)
        with Dataset(base_path+str(year)+'.nc') as data:
            nrows = data['data'].shape[0]
            chunks = int(np.floor(nrows/k))
            row = i*chunks
            centroids[i,:] = data['data'][row,:]
        print('OBTAINED', year, row)
        t2 = time.time()
        print('TIME FOR CENTROID', i, t2-t1,'\n')
    print('DONE')
    return centroids
k=int(sys.argv[1])
init_path = '/discover/nobackup/dtropf/data/ncdata/' 
centroids = init_centroids(k, init_path, 1984, 2012)
for i in range(k):
    print(i, np.sum(centroids[i,:]), np.amax(centroids[i,:]), np.amin(centroids[i,:]))
print('CENTROID SHAPE', centroids.shape)
np.save(os.getcwd()+'/k{}'.format(k), centroids)

# Create netcdf version
dataset = Dataset('{}/k{}.nc'.format(os.getcwd(), k), 'w', format='NETCDF4_CLASSIC')

# Create dimensions
kdim = dataset.createDimension('kdim', k)
pctaudim = dataset.createDimension('pctaudim', 42)

# Create variables
var_initcentroids = dataset.createVariable('init_centroids', np.float32, ('kdim', 'pctaudim'))

# Save data and close file
var_initcentroids[:,:] = centroids[:,:]
dataset.close()
