library(tidyverse)

param_tb = as_tibble(expand.grid(mutation_rate = c(1e-8,1e-9,1e-10), num_sim_2 = 1:2, sampling = c('fixed', 'proportional'), sample_size = c(50,100,200), num_sim_1 = 1:20))

param_tb$input_file = paste('sample_size_', param_tb$sample_size, '_', param_tb$sampling, '_', floor((as.numeric(rownames(param_tb))+5)/6), '.rda', sep = "")
#param_tb$input_file_original = paste('sample_size_', param_tb$sample_size, '_', param_tb$sampling, '_', floor((as.numeric(rownames(param_tb))+2)/3), '.rda', sep = "")

write.table(as.data.frame(param_tb),'/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/output/param_tb.txt',col.names = T,row.names = F,quote = F)

results <- data.frame(matrix(nrow=0, ncol = 8))
colnames(results) <- c('sample_size', 'sampling', 'mutation_rate', 'input_file', 'job_id', 'num_mut', 'subsetted_num_mut', 'KC0')

write.table(results,'/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/output/results.txt',col.names = T,row.names = F,quote = F)