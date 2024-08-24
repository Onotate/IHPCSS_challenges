#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --job-name 2019_serial_small

source teachsetup

set -x
./run.sh ${1} serial small ${2}
