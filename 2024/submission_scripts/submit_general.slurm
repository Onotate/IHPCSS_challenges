#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node 1
#SBATCH --time=00:00:30

set -x
./bin/serial/main_serial