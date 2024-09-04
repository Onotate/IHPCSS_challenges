#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --job-name 2019_hybrid_small

source teachsetup

set -x
./run.sh ${1} hybrid_cpu small ${2}
