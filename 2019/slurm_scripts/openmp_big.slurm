#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=2019_opm_big
#SBATCH --ntasks-per-node 16
#SBATCH --time=00:10:00

source teachsetup

set -x
./run.sh ${1} openmp big ${2}
