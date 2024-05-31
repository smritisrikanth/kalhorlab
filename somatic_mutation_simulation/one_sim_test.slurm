#!/bin/bash
  
#SBATCH -J somatic_mutation_one_simulation_test
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8                            # one CPU core per task
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssrikan2@jhu.edu

# Load software
ml anaconda
conda activate r4-base

# Run R script with a command line argument
Rscript one_simulation_test.R