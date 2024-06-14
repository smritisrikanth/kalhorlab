#!/bin/bash

#SBATCH -J somatic_mut_gen
#SBATCH --partition=defq
#SBATCH --array=1-200                    # how many tasks in the array
#SBATCH --cpus-per-task=4                          # one CPU core per task
#SBATCH -t 02:00:00
#SBATCH -o try-%j-%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssrikan2@jhu.edu

# Load software
ml anaconda
conda activate r4-base

# Run R script with a command line argument
Rscript somatic_mut_gen.R
