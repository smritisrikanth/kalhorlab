#!/bin/bash
  
#SBATCH -J somatic_mut_sim_create_param_tb
#SBATCH --partition=defq
#SBATCH -c 1                            # one CPU core per task
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssrikan2@jhu.edu

# Load software
ml anaconda
conda activate r4-base

# Run R script with a command line argument
Rscript make_param_tb.R