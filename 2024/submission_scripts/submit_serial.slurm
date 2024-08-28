#!/bin/bash

# The name of the output file
#SBATCH -o 2024_serial.out

# Jobs are capped at 30 seconds (Your code should run for ~10 seconds anyway)
#SBATCH -t 00:01:00

# The number of nodes (at most 2)
#SBATCH -N 1

# The number of MPI processes per node
#SBATCH --ntasks-per-node=1

source teach_setup

set -x
./bin/serial/main_serial
