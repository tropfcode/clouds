#!/bin/bash

#SBATCH --job-name=h1
#SBATCH --output=out_eqarea.txt
#SBATCH --ntasks=2
#SBATCH --qos=debug


module purge
module load other/SSSO_Ana-PyD/SApd_2.4.0_py2.7
cd /discover/nobackup/dtropf/dataproduct/equal_area/
f2py -c hproduct_analysis.f90 -m analysispy
mpirun -np 2 python ./hproduct.py 2005 ./finalcentroids_k10_v5 equal_area
