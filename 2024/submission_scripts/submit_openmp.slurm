#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 16
#SBATCH --time=00:01:00
#SBATCH -o /scratch/l/lcl_uotcscd92/lcl_uotcscd92s2111/openmp.txt

set -x
source teach_setup
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./bin/openmp/main_openmp