#!/bin/bash
#SBATCH --job-name=wsmaps_largev1
#SBATCH --output=wsmaps_fortran_largev1_%a.txt
#SBATCH --array=16,17,18,19,20
#SBATCH --ntasks=29
#SBATCH --ntasks-per-node=15
#SBATCH --time=60:00

cd /discover/nobackup/dtropf/analysis/wsmaps_new/
module purge
module load other/SSSO_Ana-PyD/SApd_2.4.0_py2.7
f2py -c analysis.f90 -m analysispy
v=1
k=$SLURM_ARRAY_TASK_ID
year_start=1984
year_end=2012
saveid=wsmaps_k${k}_v${v}
data_path=/discover/nobackup/dtropf/from_laptop/v${v}_large/out${k}.nc
mpirun -ppn 15 -np 29 python ./wsmaps.py $year_start $year_end $data_path $saveid
