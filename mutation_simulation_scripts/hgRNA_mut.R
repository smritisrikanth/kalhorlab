signal_func_hgRNA = function(x_vec) {
  val = 5e-1 * ((1/6)* pmin(pmax(0, x_vec-1), 4) + (-1/6) * pmax(0, x_vec-5))
  out_vec = pmax(val, 0)
  #out_vec[x_vec < 1] = 0
  out_vec
}

signal_func_hgRNA_decreased = function(x_vec) {
  val = 5e-2 * ((1/6)* pmin(pmax(0, x_vec-1), 4) + (-1/6) * pmax(0, x_vec-5))
  out_vec = pmax(val, 0)
  #out_vec[x_vec < 1] = 0
  out_vec
}

# signal_func_hgRNA_flat = function(x_vec) {
#   val = rep(3e-1, length(x_vec))
#   val
# }

chrmat_to_onehot <- function(x, include_unmutated = F) {
  onehot_list = lapply(1:ncol(x), function(i) {
    y = x[, i]
    y[is.na(y)] = "0"
    if (all(y == "0")) {
      return(NULL)
    }
    levels = sort(unique(y))
    out_mat = 1 * do.call(cbind, lapply(levels, function(z) y ==
                                          z))
    colnames(out_mat) = paste0("site", i , "_", levels)
    if (!include_unmutated) {
      out_mat = out_mat[, colnames(out_mat) != paste0("site", i , "_0"), drop = F]
    }
    out_mat
  })
  out = do.call(cbind, onehot_list)
  if (!is.null(out)) {
    rownames(out) = rownames(x)
  }
  out
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
library(purrr)
setwd('/home/ssrikan2/data-kreza1/smriti/qfm2')
devtools::load_all()

mut_p = readRDS("./metadata//mut_p_marc1.rds")

num_hgRNA = 5

mut_p$mut_rate = mut_p$mut_rate[1:num_hgRNA]
mut_p$recur_prob = mut_p$recur_prob[1:num_hgRNA]
mut_p$recur_vec_list = mut_p$recur_vec_list[1:num_hgRNA]

mut_p$mut_rate_func = map(1:length(mut_p$mut_rate), function(i) {
  mr = mut_p$mut_rate[i]
  signal_func_hgRNA_flat = function(x_vec) {
    val = rep(mr, length(x_vec))
    val
  }
  signal_func_hgRNA_flat
})
mut_p$mut_rate = NULL

setwd('/home/ssrikan2/data-kreza1/smriti/MF_Signal_Simulation')
job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

input_folder = 'mouse_gas_cg_ss1000'
output_folder = 'mouse_gas_cg_ss1000_flat_hgRNA'

load(paste0('./',input_folder, '/count_graph_', job_id,'.rda'))

chr_mat = simulate_all_phylogeny_bb_v1(count_graph, mut_p)
chrmat_onehot = chrmat_to_onehot(x = chr_mat$chr_mat, include_unmutated = F)
mut_frac_mat = aveMatFac(chrmat_onehot, fac = get_type_from_id(rownames(chrmat_onehot)))
colnames(mut_frac_mat) = colnames(chrmat_onehot)



mut_frac_mat = as.data.frame(mut_frac_mat)
mut_frac_mat$CellType = rownames(mut_frac_mat)


mut_frac_mat = gather(mut_frac_mat, key = 'sequence', value = 'mosaic_fraction', -CellType)

mut_frac_mat$ID = as.numeric(substr(mut_frac_mat$sequence,
                                    5,
                                    as.numeric(gregexpr('R', mut_frac_mat$sequence))-2))
mut_frac_mat$sequence = substr(mut_frac_mat$sequence,
                               as.numeric(gregexpr('R', mut_frac_mat$sequence)),
                               length(mut_frac_mat$sequence))


# anl = chr_mat$allele_node_list
#
# anl = map(anl, function(id) {
#   setNames(names(id), as.vector(id))
# })

mut_frac_mat$probability = map2_dbl(mut_frac_mat$ID, mut_frac_mat$sequence, function(id, seq) {
  mut_p$recur_vec_list[[id]][[seq]]
})
mut_frac_mat$node = map2(mut_frac_mat$ID, mut_frac_mat$sequence, function(id, seq) {
  node = names(chr_mat$allele_node_list[[id]])[which(seq == chr_mat$allele_node_list[[id]])]
  if (length(node) > 1) {
    node[which.min(((count_graph$phylo_edges$out_time+count_graph$phylo_edges$in_time)[match(node, count_graph$phylo_edges$out_node)])/2)]
  } else {
    node
  }
  #anl[[id]][[seq]]
})
mut_frac_mat$node_time = ((count_graph$phylo_edges$out_time+count_graph$phylo_edges$in_time)[match(mut_frac_mat$node, count_graph$phylo_edges$out_node)])/2


# mf_to_time_tb = nest(mf_to_time_tb, data = -c(ID, sequence)) %>%
#   mutate(mosaic_fraction = map_dbl(data, function(data) {
#     mean(data$mosaic_fraction)
#   }))
# 
# mf_to_time_tb$probability = map2_dbl(mf_to_time_tb$ID, mf_to_time_tb$sequence, function(id, seq) {
#   mut_p$recur_vec_list[[id]][[seq]]
# })
# mf_to_time_tb = mf_to_time_tb[mf_to_time_tb$probability <= 0.01,]
# mf_to_time_tb$node = map2(mf_to_time_tb$ID, mf_to_time_tb$sequence, function(id, seq) {
#   node = names(chr_mat$allele_node_list[[id]])[which(seq == chr_mat$allele_node_list[[id]])]
#   if (length(node) > 1) {
#     node[which.min(((count_graph$phylo_edges$out_time+count_graph$phylo_edges$in_time)[match(node, count_graph$phylo_edges$out_node)])/2)]
#   } else {
#     node
#   }
#   #anl[[id]][[seq]]
# })
# mf_to_time_tb = select(mf_to_time_tb, c(mosaic_fraction, node)) %>% mutate(log2mf = log2(mosaic_fraction))
# mf_to_time_tb$node_time = ((count_graph$phylo_edges$out_time+count_graph$phylo_edges$in_time)[match(mf_to_time_tb$node, count_graph$phylo_edges$out_node)])/2
# mf_to_time_tb = select(mf_to_time_tb, c(log2mf, node_time))

save(mut_frac_mat, file = paste0('./', output_folder, '/mf_tables_', job_id, '.rda'))
