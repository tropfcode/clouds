#!/bin/bash

#SBATCH --job-name=counts
#SBATCH --output=counts_%a.txt
#SBATCH --array=0,1,2,3,4
#SBATCH --ntasks=12
#SBATCH --ntasks-per-node=12
#SBATCH --time=35:00

cd /discover/nobackup/dtropf/hgg_data_analysis/
module purge
module load other/SSSO_Ana-PyD/SApd_2.4.0_py2.7
mpirun -np 12 python ./diagnostic2.py "$(($SLURM_ARRAY_TASK_ID+2008))"
