library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('/home/ssrikan2/data-kreza1/smriti/qfm2')
devtools::load_all()
setwd('/home/ssrikan2/data-kreza1/smriti/MF_Signal_Simulation')

job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#load('param_tb.rda')

load(paste0('./output_ab0.15_ss1000/count_graph_',job_id,'.rda'))

# define a signal curve
#this becomes program value?
signal_func = function(x_vec) {
  val = 1e-9 * ((1/6)* pmin(pmax(0, x_vec), 4) + (-4/15) * pmax(0, x_vec-4))
  out_vec = -pmax(val, 0) + 2*10^(-9)/3
  out_vec[x_vec > 6.5] = 0
  out_vec
}

signal_func_gene_program <- function(edge_tb, x_vec) {
  
}


# somatic with signal
tr_dd = list_dd_and_tips_mod2(tr)$dd
tr_node_time = count_graph$phylo_edges$length
tr_node_start = count_graph$phylo_edges$in_time
tr_node_end = count_graph$phylo_edges$out_time
names(tr_node_time) = names(tr_node_start) = names(tr_node_end) = count_graph$phylo_edges$out_node

tr_tips = list_dd_and_tips_mod2(tr)$tips
make_mut <- function() {
  paste0(sample(c(letters, 1:9), size = 20, replace = T), collapse = "")
}


cell_mut_tb = bind_rows(map(names(tr_tips), function(node_id) {
  # constant rate
  # Lambda = mut_rate * tr_node_time[[node_id]]
  # signal
  Lambda = integrate(signal_func, tr_node_start[node_id], tr_node_end[node_id])$value
  n_mut = rbinom(1, 10^9, 1-exp(-Lambda))
  # this version value is number of occurrences along the edge
  # tibble(cell = tr_tips[[node_id]],
  #        mut = make_mut(),
  #        value = n_mut)
  # this version generates one row for each mut that occurred
  if (n_mut > 0) {
    out_tb = expand.grid(node = node_id,
                         cell = tr_tips[[node_id]],
                         mut = map_chr(1:n_mut, function(i) {
                           make_mut()
                         }))
  } else {
    out_tb = NULL
  }
  out_tb
  #return(NULL)
}))

cell_mut_tb$value = 1
m <- spread(cell_mut_tb[2:4], mut, value, fill = 0)
# m$type <- m$cell
# m <- m[,c(1,16994,2:16993)]
chr_mat = as.matrix(m[-1])
mf_vec = colMeans(chr_mat)

save(mf_vec, file = paste0('./output2_ab0.15_ss1000/mf_vec_', job_id, '.rda'))


