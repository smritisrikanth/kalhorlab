#!/bin/bash

#SBATCH -J tree_reconstruction
#SBATCH --partition=defq
#SBATCH --array=1-200                   # how many tasks in the array
#SBATCH --cpus-per-task=4                          # one CPU core per task
#SBATCH -t 01:00:00
#SBATCH -o a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssrikan2@jhu.edu

# Load software
ml anaconda
conda activate r4-base

# Run R script with a command line argument
Rscript tree_reconstruction_simulation.R