#!/bin/bash

# Example of running R script with a job array

#SBATCH -J somatic_mutation_simulation
#SBATCH --partition=defq
#SBATCH --array=1-720                    # how many tasks in the array
#SBATCH --cpus-per-task=8                            # one CPU core per task
#SBATCH -t 3:00:00
#SBATCH -o try-%j-%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssrikan2@jhu.edu

# Load software
ml anaconda
conda activate r4-base

# Run R script with a command line argument
Rscript somatic_mut_sim.R $SLURM_ARRAY_TASK_ID