#!/bin/bash
#SBATCH --job-name=wsmaps_k11
#SBATCH --output=wsmaps_k11_%a.txt
#SBATCH --array=6
#SBATCH --ntasks=35
#SBATCH --ntasks-per-node=5
#SBATCH --qos=debug

k=11
v=$SLURM_ARRAY_TASK_ID
year_start=1983
year_end=2017
saveid=wsmaps_k${k}_v${v}_final
data_path=/discover/nobackup/dtropf/from_laptop/k11_set2/out${v}.nc

module purge
module load other/SSSO_Ana-PyD/SApd_2.4.0_py2.7
cd /discover/nobackup/dtropf/current/wsmaps_new/
f2py -c analysis.f90 -m analysispy
mpirun -ppn 5 -np 35 python ./wsmaps.py $year_start $year_end $data_path $saveid
