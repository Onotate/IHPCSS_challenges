#!/bin/bash


/*
 * This script generates a set of SLURM job scripts for the hybrid version of the program.
 * The number of nodes, ntasks-per-node, and cpus-per-task are varied.
 * The generated job scripts are saved in the current directory.
 */

num_nodes=(4)
ntasks_per_node=(1 2 4 8 16)
cpus_per_task=(1 2 4 8 16)

for nodes in "${num_nodes[@]}"; do
  for ntasks in "${ntasks_per_node[@]}"; do
    total_tasks=$((nodes * ntasks))  # Calculate total number of tasks

    # Ensure that the condition ntasks * nodes == 8 is satisfied
    if [[ $total_tasks -eq 8 ]]; then
      # Loop through all possible values for cpus-per-task
      for cpus in "${cpus_per_task[@]}"; do
        job_script="slurm_job_${nodes}_${ntasks}_${cpus}.sh"
        echo "Running job script with nodes=$nodes, ntasks=$ntasks, cpus=$cpus"
    sbatch <<EOL >"$job_script"
#!/bin/bash
#SBATCH -o 2024_hybrid_big_${nodes}_${ntasks}_${cpus}.out
#SBATCH -t 01:00:00
#SBATCH -N $nodes
#SBATCH --ntasks-per-node=$ntasks
#SBATCH --cpus-per-task=$cpus

# The number of OpenMP threads. If using MPI, it is the number of OpenMP threads per MPI process
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK

# Load the required modules
source teach_setup

# OpenMPI issue workaround
export OMPI_MCA_btl=^openib

# Execute the program
mpirun -n \$SLURM_NTASKS ./bin/hybrid_cpu_small
EOL

     done
    fi
  done
done
