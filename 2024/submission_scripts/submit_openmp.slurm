#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 16
#SBATCH --output=2024_openmp.out
#SBATCH --time=00:01:00

source teach_setup
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./bin/openmp/main_openmp