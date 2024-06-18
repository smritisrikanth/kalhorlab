list_dd_and_tips_from_edge_tb <- function(edge_df) {
  tb_tips = edge_df$to[!edge_df$to %in% edge_df$from]
  tr_dd_new = split(edge_df$to, edge_df$from)
  m_seq = make_merge_sequences_mod2(dplyr::rename(edge_df,
                                                  out_node = to, in_node = from), dd_list = tr_dd_new,
                                    tip_id = tb_tips)
  tr_tip_list_new = map(tb_tips, function(x) x)
  names(tr_tip_list_new) = tb_tips
  
  tr_down_list = map(tb_tips, function(x) x)
  names(tr_down_list) = tb_tips
  
  for (x in m_seq) {
    x_tips = map(x, function(dd) {
      purrr::reduce(tr_tip_list_new[dd], c)
    })
    x_downs = map(names(x), function(x_name) {
      dd = x[[x_name]]
      c(purrr::reduce(tr_down_list[dd], c), x_name)
    })
    names(x_downs) = names(x)
    tr_tip_list_new = append(tr_tip_list_new, x_tips)
    tr_down_list = append(tr_down_list, x_downs)
  }
  list(tips = tr_tip_list_new,
       downs = tr_down_list,
       dd = tr_dd_new)
}

tb_to_func_vec <-  function(in_tb) {
  in_tb_cp = in_tb
  out_func <- function(t_vec) {
    map_dbl(t_vec, function(t_val) {
      assertthat::assert_that(t_val <= max(in_tb_cp$time), msg = "time out of range")
      assertthat::assert_that(t_val >= min(in_tb_cp$time), msg = "time out of range")
      t_ind = which.min(abs(in_tb_cp$time - t_val))
      in_tb_cp$value[t_ind]
    })
  }
  out_func(in_tb_cp$time[[1]])
  out_func
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
make_mut <- function() {
  paste0(sample(c(letters, 1:9), size = 20, replace = T), collapse = "")
}

# tb = tibble(time = seq(edges$start_time[1], edges$end_time[1], length = length(edges$program1[[1]])),
#             value = edges$program1[[1]])
#
# gene_func = tb_to_func_vec(tb)

# integrate(gene_func, lower = 0.1, upper = 0.5)$value

# tr_dd = list_dd_and_tips_mod2(tr)$dd
tr_dd = split(edge_tb$to, edge_tb$from)
tr_node_time = edge_tb$length
tr_node_start = edge_tb$from_time
tr_node_end = edge_tb$to_time
names(tr_node_time) = names(tr_node_start) = names(tr_node_end) = edge_tb$to

tr_tips = list_dd_and_tips_from_edge_tb(edge_tb)$tips

program_val = edges %>% select(program3) %>% unnest(cols = program3) %>% pull(program3)
median_program = median(program_val)
base_mean = log(10^(0)) - median_program
base = rnorm(1,base_mean,1)

#V2
program_val = edge_tb %>% select(program1) %>%
  unnest(cols = program1) %>% select(value)
median_program = median(program_val$value)
base_mean = log(10^(0)) - median_program
base = rnorm(1,base_mean,1)

i = 0
library(future)
plan(multisession, workers = 12)
cell_mut_tb = bind_rows(future_map(names(tr_tips), function(node_id) {
  # constant rate
  # Lambda = mut_rate * tr_node_time[[node_id]]
  # signal
  # if (i %% 1000 == 0) {
  #         message(i)
  # }
  i <<- i+1
  if (!(node_id %in% edge_tb$to)) {
    # message('root')
    branch = edges[edges$in_node == 'Root',]
    program_val_tb = tibble(time = seq(0, edge_tb$from_time[1], length = length(edges$program3[[1]])),
                            value = edges$program3[[1]])
  } else {
    program_val_tb = edge_tb$program3[edge_tb$to == node_id][[1]]
  }
  
  program_val_tb$value = exp(base + program_val_tb$value)
  signal_func = tb_to_func_vec(program_val_tb)
  
  
  if (!(node_id %in% edge_tb$to)) {
    Lambda = integrate(signal_func,
                       program_val_tb$time[[1]],
                       program_val_tb$time[[2]],
                       subdivisions = 1000)$value
  } else {
    Lambda = integrate(signal_func,
                       tr_node_start[node_id],
                       tr_node_end[node_id],
                       subdivisions = 1000)$value
  }
  n_mut = rbinom(1, 10^9, 1-exp(-Lambda/10^9))
  # mean(n_mut)
  # this version value is number of occurrences along the edge
  # tibble(cell = tr_tips[[node_id]],
  #        mut = make_mut(),
  #        value = n_mut)
  # this version generates one row for each mut that occurred
  if (n_mut > 0) {
    out_tb = tibble(expand.grid(node = node_id,
                                cell = tr_tips[[node_id]],
                                mut = map_chr(1:n_mut, function(i) {
                                  make_mut()
                                })
                                , stringsAsFactors = F))
  } else {
    out_tb = NULL
  }
  out_tb
  #return(NULL)
}, .progress = T))

# check the number of mutations per edge
# mut_sum = cell_mut_tb %>% group_by(node) %>%
#         summarise(n_mut = length(unique(mut)))

#single cell table
cell_mut_tb = cell_mut_tb %>% select(tip = cell, mut, node_to = node)
cell_mut_tb = left_join(cell_mut_tb, edge_tb %>% select(node_to = to,
                                                        node_from = from,
                                                        node_from_time = from_time,
                                                        node_to_time = to_time),
                        by = 'node_to')
cell_mut_tb = left_join(cell_mut_tb, edge_tb %>% select(node_to = to, block))


cell_mut_tb$value = 1




m <- spread(cell_mut_tb %>% select(tip, mut, value), mut, value, fill = 0)
# m$type <- m$cell
# m <- m[,c(1,16994,2:16993)]
chr_mat = as.matrix(m[-1])
rownames(chr_mat) = m[[1]]

tr_cells_append = edge_tb$to[!(edge_tb$to[!edge_tb$to %in% edge_tb$from] %in% rownames(chr_mat))]
chr_mat_append = matrix(0, nrow = length(tr_cells_append), ncol = ncol(chr_mat))
rownames(chr_mat_append) = tr_cells_append
chr_mat = rbind(chr_mat, chr_mat_append)
#assertthat::assert_that(nrow(chr_mat) == length(tr$tip.label))

tips_type = map_chr(rownames(chr_mat), function(tip) {
  edge_tb$to_type[edge_tb$to == tip]
})
mf_mat_cell_type = aveMatFac(chr_mat, tips_type)
colnames(mf_mat_cell_type) = colnames(chr_mat)
mf_mat_cell_type = as.data.frame(mf_mat_cell_type)
mf_mat_cell_type$CellType = rownames(mf_mat_cell_type)


mf_table = gather(mf_mat_cell_type, key = 'sequence', value = 'mosaic_fraction', -CellType)
mf_table = left_join(mf_table, unique(cell_mut_tb %>% select(sequence = mut, node_to_time, node_from_time)))


#mf table with sparse matrix
library(Matrix)

load('3_13_cell_mut_tb.rda')
cell_mut_tb = left_join(cell_mut_tb, edge_tb %>% select(node_to = to, to_type, from_type))

all_tips = sort(unique(cell_mut_tb$tip))
all_tips_type = edge_tb$to_type[match(all_tips, edge_tb$to)]
all_muts = sort(unique(cell_mut_tb$mut))
mut_mat = sparseMatrix(i = match(cell_mut_tb$tip, all_tips),
                       j = match(cell_mut_tb$mut, all_muts),
                       x = rep(1, nrow(cell_mut_tb)))


# mf_mat_cell_type = t(bind_cols(map(unique(all_tips_type), function(type) {
#   mf_vec = colMeans(mut_mat[all_tips_type == type,])
#   names(mf_vec) = colnames(mut_mat)
#   mf_vec
# })))
mf_mat_cell_type = aveMatFac(mut_mat, all_tips_type)
rownames(mf_mat_cell_type) = unique(all_tips_type)
colnames(mf_mat_cell_type) = all_muts
mf_mat_cell_type = as.data.frame(mf_mat_cell_type)

mf_mat_cell_type$CellType = rownames(mf_mat_cell_type)

mf_table = gather(mf_mat_cell_type, key = 'sequence', value = 'mosaic_fraction', -CellType)
mf_table = left_join(mf_table, unique(cell_mut_tb %>% select(sequence = mut, node_to_time, node_from_time,
                                                             node_from, node_to,
                                                             to_type, from_type)))





