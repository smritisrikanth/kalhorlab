signal_func = function(x_vec) {
  val = 1e-9 * ((1/6)* pmin(pmax(0, x_vec-1), 4) + (-1/6) * pmax(0, x_vec-5))
  out_vec = pmax(val, 0) #+ 2*10^(-9)/3
  #out_vec[x_vec < 1] = 0
  out_vec
}

signal_func_decreased = function(x_vec) {
  val = 1e-10 * ((1/6)* pmin(pmax(0, x_vec-1), 4) + (-1/6) * pmax(0, x_vec-5))
  out_vec = pmax(val, 0) #+ 2*10^(-9)/3
  #out_vec[x_vec < 1] = 0
  out_vec
}

signal_func_flat = function(x_vec) {
  val = rep(5e-10, length(x_vec))
  val
}

make_mut <- function() {
  paste0(sample(c(letters, 1:9), size = 20, replace = T), collapse = "")
}
aveMatFac <- function(mat, fac) {
  # need to be able to handle character or numeric
  if (class(fac) != "factor")
    fac <- factor(fac)
  rown <- length(levels(fac))
  coln <- dim(mat)[2]
  out <- matrix(, rown, coln)
  ind <- as.numeric(fac)
  for (i in 1:rown) {
    out[i,] <- colMeans(mat[ind == i, , drop = F], na.rm = T)
  }
  rownames(out) <- levels(fac)
  return(out)
}

library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('/home/ssrikan2/data-kreza1/smriti/qfm2')
devtools::load_all()
setwd('/home/ssrikan2/data-kreza1/smriti/MF_Signal_Simulation')
job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

input_folder = 'output'
output_folder = 'one_cell_somatic_mut_signal_flat'

print(job_id)
load(paste0('./',input_folder, '/count_graph_', floor(job_id/2),'.rda'))

tr = phylo_edges_to_tr(count_graph$phylo_edges)
tr_dd = list_dd_and_tips_mod2(tr)$dd
tr_node_time = count_graph$phylo_edges$length
tr_node_start = count_graph$phylo_edges$in_time
tr_node_end = count_graph$phylo_edges$out_time
names(tr_node_time) = names(tr_node_start) = names(tr_node_end) = count_graph$phylo_edges$out_node

tr_tips = list_dd_and_tips_mod2(tr)$tips

# node_id = names(tr_tips)[[1]]
cell_mut_tb = bind_rows(map(names(tr_tips), function(node_id) {
  # constant rate
  # Lambda = mut_rate * tr_node_time[[node_id]]
  # signal
  Lambda = integrate(signal_func_flat, tr_node_start[node_id], tr_node_end[node_id])$value
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
rownames(chr_mat) = m[[1]]

tr_cells_append = tr$tip.label[!tr$tip.label %in% rownames(chr_mat)]
chr_mat_append = matrix(0, nrow = length(tr_cells_append), ncol = ncol(chr_mat))
rownames(chr_mat_append) = tr_cells_append
chr_mat = rbind(chr_mat, chr_mat_append)
assertthat::assert_that(nrow(chr_mat) == length(tr$tip.label))

mf_mat_cell_type = aveMatFac(chr_mat, get_type_from_id(rownames(chr_mat)))
colnames(mf_mat_cell_type) = colnames(chr_mat)
mf_mat_cell_type = as.data.frame(mf_mat_cell_type)
mf_mat_cell_type$CellType = rownames(mf_mat_cell_type)


mf_table = gather(mf_mat_cell_type, key = 'sequence', value = 'mosaic_fraction', - CellType)
mf_table = left_join(mf_table, cell_mut_tb %>% mutate(sequence = mut))
mf_table$node_time = ((count_graph$phylo_edges$out_time+count_graph$phylo_edges$in_time)[match(mf_table$node, count_graph$phylo_edges$out_node)])/2


# mf_to_time_tb = nest(mf_table, data = -c(sequence)) %>%
#   mutate(mf = map_dbl(data, function(data) {
#     mean(data$mosaic_fraction)
#   }))
# mf_to_time_tb = left_join(cell_mut_tb, mf_to_time_tb %>% select(mut = sequence, mf)) %>%
#   mutate(log2mf = log2(mf))
# mf_to_time_tb$node_time = ((count_graph$phylo_edges$out_time+count_graph$phylo_edges$in_time)[match(mf_to_time_tb$node, count_graph$phylo_edges$out_node)])/2
# mf_to_time_tb = select(mf_to_time_tb, c(log2mf, node_time))

save(mf_table, file = paste0('./', output_folder, '/mf_table_', job_id, '.rda'))
