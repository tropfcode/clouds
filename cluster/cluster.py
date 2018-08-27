#!/usr/bin/env python
import gc
import sys
import time
import random
import logging
import numpy as np
import cluster_label
from mpi4py import MPI
from netCDF4 import Dataset

def init_centroids(base_path, year_start, year_end=None):
    """
    
    """
    k = 11
    centroids = np.zeros((11, 42))
    if year_end is None:
        year_end = year_start + 1
    for i in range(k):
        t1 = time.time()
        print('OBTAINING CENTROID', i)
        year = random.randint(year_start, year_end)
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

# MPI variables
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Setup logging
logging.basicConfig(filename='{}.log'.format(string(rank)), level=logging.DEBUG)
logging.getLogger().addHandler(logging.StreamHandler()) # Have logging also print to stdout

# Paths to data
init_path = '/discover/nobackup/dtropf/cluster/' # .nc files
data_path = '/discover/nobackup/dtropf/cluster/transformed_data/transformed_data/data/' #.npy files

# Initialize variables
ndim = 42
ncent = 11

# Command line arguments
niter = int(sys.argv[1]) 
year_start = int(sys.argv[2])
year_end = int(sys.argv[3])

centroids = np.zeros((ncent, ndim))
years = [i for i in range(year_start, year_end)]
vectors = np.load(data_path+str(years[rank])+'.npy')
nvec = vectors.shape[0]

# Start program
t1 = time.time()

# Initialize centroids from even distribution
if rank == 0:
    centroids = init_centroids(init_path, year_start, year_end)
    logging.info('SUM FOR INITALIZED CENTROIDS '+str(np.sum(centroids)))
else:
    centroids = np.empty((ncent, ndim))

# Broadcast initalize centroids to all nodes
comm.Bcast(centroids, root=0)

# Check Broadcast done properly
if rank != 0:
    logging.info('CENTROID SUM FOR RANK '+str(rank)+''+str(np.sum(centroids)))


for i in range(niter):
    t3 = time.time()
    
    # K-means clustering
    clusters, labels = cluster_label.cluster_label(vectors, centroids, nvec, ncent, ndim)
    
    # Gather clusters and labels from each task to task rank=0
    recvclusters = comm.gather(clusters, root=0)
    recvlabels = comm.gather(labels, root=0)

    # Calculate new centroids
    if rank ==0:
        cluster_sum = np.sum(recvclusters, axis=0)
        label_sum = np.sum(recvlabels, axis=0)
        
        for j in range(ncent):
            centroids[j,:] = cluster_sum[j,:]/label_sum[j]

    # Synchronize and Broadcast new centroids
    comm.Barrier()
    comm.Bcast(centroids, root=0)
    t4 = time.time()
    logging.info('TOTAL TIME FOR ITERATION '+str(i)+' OF RANK '+str(rank), str(t4-t3))
    comm.Barrier()

# Save resulting centroids and labels
for i in range(size):
    if rank == 0:
        logging.info('DONE WITH RANK 0 SAVING CENTROIDS')
        np.save('centroids_{}_{}_{}iter'.format(str(year_start), str(year_end), str(niter)), centroids)
        np.save('labels_{}_{}_{}iter'.format(str(year_start), str(year_end), str(niter)), label_sum)
    elif i == rank:
        t2 = time.time()
        logging.info('DONE WITH ENTIRE PROGRAM RANK '+str(rank)+''+str(t2-t1))
