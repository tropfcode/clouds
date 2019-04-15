#!/bin/bash

#SBATCH --job-name=newv3_9fortran
#SBATCH --output=output_newv3_9%a_.txt
#SBATCH --array=9
#SBATCH --constraint=hasw
#SBATCH --ntasks=29
#SBATCH --ntasks-per-node=4
#SBATCH --time=120:00

niter=300
ncent=$SLURM_ARRAY_TASK_ID

cd /discover/nobackup/dtropf/from_laptop/
module purge
module load comp/intel-18.0.0.128 mpi/impi-18.0.0.128 other/wrf-deps
mpif90 module.f90 cluster.f90 -L/usr/local/other/SLES11.3/wrf-deps/intel-18.0.0.128/lib -lnetcdff -L/usr/local/other/SLES11.3/wrf-deps/intel-18.0.0.128/lib -lnetcdf -lnetcdf -lhdf5_hl -lhdf5 -lz
mpirun -ppn 4 -np 29 ./a.out $niter $ncent
