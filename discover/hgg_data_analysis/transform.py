"""
PROGRAM DESCRIPTION:
   Transform results from diagnostic.py into single file.

PROGRAM IMPLEMENTATION:
   Load the 12 files from diagnostic2.py and put them into
a single file in chronological order. Save the resulting file in
netCDF format with two variables: 
   1. month_counts: 1D array with file counts for each month
   2. data_profile: 2D array of shape (n_files, 3). n_files is the total
number of 3-hour HGG data files in a single year. The second dimension covers
the counts of each category of gridbox (clouds, no clouds, no data).
"""

import time
import numpy as np
from netCDF4 import Dataset

years = [i for i in range(1985,2013)]
for year in years:
   t1 = time.time()
   dataset = Dataset('hgg_profile_{}.nc'.format(year), 'w', format='NETCDF4_CLASSIC')
   data_list = [np.load('counts_{}_{}.npy'.format(year,i)) for i in range(12)]
   nfiles = sum([data.shape[0] for data in data_list])
   print('NUMBER OF FILES FOR YEAR {}: '.format(nfiles), nfiles)
   
   # Count number of files from each month
   mcount = [data.shape[0] for data in data_list]

   # Arrange all data into single numpy array
   index = 0
   alldata = np.zeros((nfiles, 3), dtype=np.int32)
   for data in data_list:
      alldata[index:index+data.shape[0],:] = data[:,:]
      index += data.shape[0]

   # Create netcdf dimensions
   filesdim = dataset.createDimension('nfiles', nfiles)
   cloudsdim = dataset.createDimension('ndatatypes', 3)
   monthsdim = dataset.createDimension('nmonths', 12)

   # Create netcdf variables
   var_monthcount = dataset.createVariable('month_counts', np.int32, ('nmonths'))
   var_dataprofile = dataset.createVariable('data_profile', np.int32, ('nfiles', 'ndatatypes'))

   # Place data in variables and close netcdf file
   var_monthcount[:] = np.asarray(mcount, dtype=np.int32)
   var_dataprofile[:,:] = alldata[:,:]
   dataset.close()

   # Print message
   t2 = time.time()
   print('DONE WITH YEAR {} TIME {}'.format(year, t2-t1))
