#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=test.txt
#SBATCH --ntasks=35
#SBATCH --ntasks-per-node=4
#SBATCH --time=2:00:00
#SBATCH --qos=giss

niter=300
ncent=10
startyear=1983
init_path='/discover/nobackup/dtropf/from_laptop/final_k10/k10_v5.nc'
data_path_base='/discover/nobackup/dtropf/data/ncdata/'

cd /discover/nobackup/dtropf/from_laptop/tmp_for_permenant/
module purge
module load comp/intel-18.0.0.128 mpi/impi-18.0.0.128 other/wrf-deps
mpif90 module.f90 netmpi.f90 -L/usr/local/other/SLES11.3/wrf-deps/intel-18.0.0.128/lib -lnetcdff -L/usr/local/other/SLES11.3/wrf-deps/intel-18.0.0.128/lib -lnetcdf -lnetcdf -lhdf5_hl -lhdf5 -lz
mpirun -ppn 4 -np 35 ./a.out $niter $ncent $startyear $init_path $data_path_base
