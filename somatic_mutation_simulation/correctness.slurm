#!/bin/bash

# Example of running R script with a job array

#SBATCH -J tree_comparison
#SBATCH --partition=defq
#SBATCH --array=1-240                    # how many tasks in the array
#SBATCH --cpus-per-task=4                            # one CPU core per task
#SBATCH -t 00:02:00
#SBATCH -o try-%j-%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssrikan2@jhu.edu

# Load software
ml anaconda
conda activate r4-base

# Run R script with a command line argument
Rscript tree_comparison.R $SLURM_ARRAY_TASK_ID