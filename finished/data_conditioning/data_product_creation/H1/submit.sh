#!/bin/bash

#SBATCH --job-name=h1
#SBATCH --output=out.txt
#SBATCH --ntasks=3
#SBATCH --qos=debug

module purge
module load other/SSSO_Ana-PyD/SApd_2.4.0_py2.7
cd /discover/nobackup/dtropf/dataproduct/h1/
f2py -c hproduct_analysis.f90 -m analysispy
mpirun -np 3 python ./hproduct.py 2000 ./finalcentroids_k10_v5 h1
