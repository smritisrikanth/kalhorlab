#!/bin/bash

#SBATCH -J new_lin_spec
#SBATCH --partition=parallel               # how many tasks in the array
#SBATCH --cpus-per-task=4                          # one CPU core per task
#SBATCH -t 02:00:00
#SBATCH -o a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssrikan2@jhu.edu

# Load software
ml anaconda
conda activate r4-base

# Run R script with a command line argument
Rscript lin_spec_new_rockfish.R