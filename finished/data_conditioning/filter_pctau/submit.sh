#!/bin/bash

#SBATCH --job-name=makehdf2015_2017
#SBATCH --output=makehdf2015_2017%a.txt
#SBATCH --constraint=hasw
#SBATCH --ntasks=3
#SBATCH --qos=debug

cd /discover/nobackup/dtropf/hgg_data_code/
module purge
module load other/SSSO_Ana-PyD/SApd_2.4.0_py2.7
mpirun -np 3 python ./filter.py 2015 2017
