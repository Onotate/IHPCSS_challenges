#!/bin/bash

# The name of the output file
#SBATCH -o 2024_hybrid.out

# Jobs are capped at 1 minute
#SBATCH -t 00:01:00

# The number of nodes (at most 2)
#SBATCH -N 3

# The number of MPI processes per node
#SBATCH --ntasks-per-node=2

# The number of OpenMP threads per MPI process
#SBATCH --cpus-per-task=8

# The number of OpenMP threads. If using MPI, it is the number of OpenMP threads
# per MPI process
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Load the modules needed
source teach_setup

# Testing this flag, something to do with Infiniband
# OpenMPI issue: https://github.com/open-mpi/ompi/issues/10693
export OMPI_MCA_btl=^openib

# Execute the program
# Make sure to run the executable in the parent directory and target correct binary
mpirun -n $SLURM_NTASKS --bind-to none --report-bindings ./bin/hybrid/main_hybrid
