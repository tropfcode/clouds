#!/bin/bash

#SBATCH --job-name=arrayjob
#SBATCH --output=arrayjob_%A_%a_iter.txt
#SBATCH --array=50,100,200,300
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=2
#SBATCH --nodes=5
#SBATCH --time=300:00
#SBATCH --qos=giss
#SBATCH --mem=30G

cd /discover/nobackup/dtropf/mpi4py/array_batch/
module purge
module load other/SSSO_Ana-PyD/SApd_2.4.0_py2.7
f2py -c cluster_label.f90 -m cluster_label
mpirun -ppn 2 -np 10 python ./cluster.py $SLURM_ARRAY_TASK_ID 1984 1994
