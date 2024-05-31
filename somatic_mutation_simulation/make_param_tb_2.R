library(tidyverse)

param_tb = as_tibble(expand.grid(sampling = c('fixed', 'proportional'), sample_size = c(50,100,200), num_sim_1 = 1:20))

param_tb$input_file = paste('sample_size_', param_tb$sample_size, '_', param_tb$sampling, '_', as.numeric(rownames(param_tb)), '.rda', sep = "")
#param_tb$input_file_original = paste('sample_size_', param_tb$sample_size, '_', param_tb$sampling, '_', floor((as.numeric(rownames(param_tb))+2)/3), '.rda', sep = "")

write.table(as.data.frame(param_tb),'/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/correctness_output/param_tb_2.txt',col.names = T,row.names = F,quote = F)

results <- data.frame(matrix(nrow=0, ncol = 8))
colnames(results) <- c('sample_size', 'sampling', 'mutation_rate', 'input_file', 'job_id', 'reconstruction_num', 'percent_correct', 'percent_complete')

write.table(results,'/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/correctness_output/results_2.txt',col.names = T,row.names = F,quote = F)