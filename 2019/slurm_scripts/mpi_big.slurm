#!/bin/bash

#SBATCH --nodes=4
#SBATCH --job-name=2019_mpi_big
#SBATCH --ntasks-per-node 16
#SBATCH --time=00:06:00

source teachsetup

set -x
./run.sh ${1} mpi big ${2}
