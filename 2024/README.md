# IHPCSS Programming challenge 2024

## Presentation slides
Tips and guidelines for the programming challenge are given in the presentation slides, available on the Moodle page at https://www.hpc-training.org/moodle/course/view.php?id=93#section-6.


### Submit a job
If you do not use GPUs yet, you can submit your job to CPU nodes using the submission script: `submission_cpu_node.slurm`. However, if you use GPU programming, you must use the `submission_gpu_node.slurm` submission script to submit to a GPU node.

To submit a script, you need to type `sbatch your_submission_script.slurm`.

### Check for your job completion
To check the status of your current job(s), type `squeue -u $USER`. As a convenience, you can type `watch squeue -u $USER`: it will show the current list of jobs in queue and/or being executed, and refresh this list every 2 seconds by default.

## Graph generators
The graph processed is generated from inside the application. You will see that the code has two graph generators:
- `generate_nice_graph`: generates a graph that will easily offer parallelisation gains.
- `generate_sneaky_graph`: generates a graph that will require some analysis to get parallelisation gains. **This is the graph that will be used in submissions.**
