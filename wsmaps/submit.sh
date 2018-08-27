#!/bin/bash

#SBATCH --job-name=wsmap
#SBATCH --output=arrayjob_%A_%a_iter.txt
#SBATCH --array=50,100,200,300
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=10
#SBATCH --nodes=1
#SBATCH --time=300:00
#SBATCH --qos=giss
#SBATCH --mem=5G

cd /discover/nobackup/dtropf/mpi4py/array_batch/wsmaps/
module purge
module load other/SSSO_Ana-PyD/SApd_2.4.0_py2.7
f2py -c analysis.f90 -m analysispy
id=${SLURM_ARRAY_TASK_ID}iter
year_start=1984
year_end=1994
mpirun -ppn 10 -np 10 python ./wsmaps.py $year_start $year_end $id centroids_${year_start}_${year_end}_${id}.npy
