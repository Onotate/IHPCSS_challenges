#!/bin/bash

# The name of the output file
#SBATCH -o 2024_mpi.out

# Jobs are capped at 1 minute
#SBATCH -t 00:01:00

# The number of nodes (at most 2)
#SBATCH -N 1

# The number of MPI processes per node
#SBATCH --ntasks-per-node=16

# The number of OpenMP threads per MPI process
#SBATCH --cpus-per-task=1

# Load the modules needed
source teach_setup

# Testing this flag, something to do with Infiniband
# OpenMPI issue: https://github.com/open-mpi/ompi/issues/10693
export OMPI_MCA_btl=^openib

# Execute the program
# Make sure to run the executable in the parent directory and target correct binary
mpirun -n $SLURM_NTASKS --bind-to none --report-bindings ./bin/mpi/main_mpi
