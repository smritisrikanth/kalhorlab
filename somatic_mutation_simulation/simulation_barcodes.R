sink(stdout(), type="message")

# simulate lineage barcodes
library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('/home/ssrikan2/data-kreza1/smriti/qfm2')
#devtools::load_all()

#job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
job_id = 49
input_path <- '/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/output/'
param_table <- read.table('/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/output2/param_tb_2.txt',header = T)
load(paste0(input_path,param_table$input_file[job_id]))
global_step_size = 0.01
global_target_time = 9.0
# mut_p = readRDS("./metadata//mut_p_marc1.rds")

# # simulate barcodes
# chr_mat = simulate_all_phylogeny_bb(count_graph, mut_p)

# # reconstruct phylogeny
library(furrr)
# plan(multisession, workers = 12)
# tr2 = phylotime(sc_mat = chr_mat, mut_p = mut_p,
#                  t_total = count_graph$target_time - count_graph$phylo_edges$length[1],
#                  parallel = T)

# x = cophenetic(tr_r)
# y = cophenetic(tr)
# tibble(a = c(x[rownames(y), colnames(y)]),
#        b = c(y)) %>%
#         sample_n(size = 1e4) %>%
#         ggscatter(x = "a", y = "b") + geom_abline()

# reconstruct fate map using reconstructed phylogeny
res3 = ice_fase_mod(tr = tr2,
                    sc_celltypes = get_type_from_id(tr2$tip.label),
                    root_time = count_graph$phylo_edges$out_time[1],
                    total_time = 9.0 - count_graph$phylo_edges$out_time[1])
res4 = collaspe_di(res3, 0.2)
res5 = ice_fase_mutlifurc(tr = tr2,
                          sc_celltypes = get_type_from_id(tr2$tip.label),
                          root_time = count_graph$phylo_edges$out_time[1],
                          total_time = 9.0 - count_graph$phylo_edges$out_time[1],
                          gr = res4$gr)
save(count_graph, tr, res0, res1, res2, res3, res4, res5, 
     file = paste0("/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/output2/",job_id,'-',param_table$input_file[job_id]))








