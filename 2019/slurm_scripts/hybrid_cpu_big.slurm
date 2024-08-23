#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks-per-node 2
#SBTACH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --job-name 2019_hybrid_big
set -x
./run.sh ${1} hybrid_cpu big ${2}
