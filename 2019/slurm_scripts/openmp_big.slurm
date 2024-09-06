#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=2019_opm_big
#SBATCH --cpus-per-task=16
#SBATCH --time=00:20:00

source teachsetup

set -x
./run.sh ${1} openmp big ${2}
