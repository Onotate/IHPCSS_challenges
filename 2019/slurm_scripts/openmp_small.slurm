#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=2019_omp_small
#SBATCH --ntasks-per-node 4
#SBATCH --time=00:01:00

source teachsetup

set -x
./run.sh ${1} openmp small ${2}
