"""
PROGRAM DESCRIPTION:
   Obtain k centroids from raw HGG data and save in numpy and
netCDF file formats.

PROGRAM IMPLEMENTATION:
   The number of centroids desired to be obtained, k, is provided as
command line input. Randomly choose k years from the set of years
available (1984 to 2012) to obtain centroids. Obtain the HGG profile
data for each year (from hgg_profile_data program diagnostics.py) 
and randomly choose which raw HGG data file to obtain the centroid from.
   So far, all steps in the centroid choosing process are temporally
random. Finally, uniformly choose the centroid spatially by dividing
the number of gridboxes from the randomly chosen HGG data file by k.
Step through the file by an amount equal to that calculation to finally
choose the centroid. This helps prevent biasing of geographic regions with
fewer gridboxes (such as poles vs. equator).
"""

import os
import sys
import math
import time
import random
import numpy as np
from netCDF4 import Dataset

def init_centroids(k):
   """
   Obtain k centroids from raw HGG data that act as initialization centroids
   for k-means clustering. The centroids are sampled geographically from pole to
   pole and randomly in time (see PROGRAM IMPLEMENTATION above)

   PARAMETERS
   ----------
   k: int
      Number of centroids to obtain
   """
   # Initalize variables
   centroids = np.zeros((k,42))
   base_path = '/discover/nobackup/dtropf/data/ncdata/'
   data_profile_path = '/discover/nobackup/dtropf/hgg_data_analysis/'

   # Randomly choose which years to obtain centroids from
   year_choices = [i for i in range(1984,2013)]
   years = [year_choices.pop(year_choices.index(random.choice(year_choices)))
            for i in range(k)]
   data_paths = [base_path+'{}.nc'.format(year) for year in years]
   print('Random years that were obtained: ', years)

   # Load profile data for selected years
   # Randomly choose which 3-hour time period within the year to obtain centroid
   profile_data = []
   for year in years:
      with Dataset(data_profile_path+'hgg_profile_{}.nc'.format(year)) as data:
         fnum = random.randint(0,data['data_profile'].shape[0])
         if fnum == 0:
            fstart = 0
         else:
            fstart = np.sum(data['data_profile'][:fnum,0])
         fnext = data['data_profile'][fnum,0]
         profile_data.append((fstart,fnext))

   # Obtain k centroids
   for i, path in enumerate(data_paths):
      t1 = time.time()
      print('OBTAINING CENTROID {}'.format(i))
      with Dataset(path) as data:
         fstart = profile_data[i][0]
         fnext = profile_data[i][1]
         step_size = int(math.floor(fnext/k))
         index = fstart+step_size*i
         centroids[i,:] = data['data'][index,:]
      t2 = time.time()
      print('Time {} for year {} row {} fstart {} fnext {} step_size {} step_size*i {} centroid {}'
            .format(t2-t1, path[-7:-3], index, fstart, fnext, step_size, step_size*i, i),'\n')
      print('CENTROID {}'.format(i), centroids[i,:])
   print('DONE')
   return centroids

def main():
   # Obtain command line argument
   k = int(sys.argv[1])
   centroids = init_centroids(k)
   
   # Print centroid information
   for i in range(k):
      sum = np.sum(centroids[i,:])
      max = np.amax(centroids[i,:])
      min = np.amin(centroids[i,:])
      print('CENTROID {} SUM MAX MIN: {} {} {}'
            .format(i, sum, max, min))

   # Save data in .npy format
   np.save(os.getcwd()+'/k{}'.format(k), centroids)

   # Create netcdf version
   dataset = Dataset(os.getcwd()+'/k{}.nc'.format(k),
                     'w', format='NETCDF4_CLASSIC')

   # Create dimensions
   kdim = dataset.createDimension('kdim', k)
   pctaudim = dataset.createDimension('pctaudim', 42)

   # Create variables
   var_initcentroids = dataset.createVariable('init_centroids', np.float32, ('kdim', 'pctaudim'))

   # Save data and close file
   var_initcentroids[:,:] = centroids[:,:]
   dataset.close()

if __name__ == '__main__':
   print('START PROGRAM')
   main()
   print('END PROGRAM')
