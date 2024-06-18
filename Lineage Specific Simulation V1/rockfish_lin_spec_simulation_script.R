#functions
assign_state_transition <- function(edge_tb, from, root = F) {
  
  i = edge_tb$from_type[edge_tb$from == from][1]
  
  j <<- j+1
  if (j %% 1000 == 0) {
    message(j)
  }
  out_cells = edge_tb$to[edge_tb$from == from]
  for (to in out_cells) {
    if (!root) {
      prev_st = edge_tb$st[edge_tb$to == from]
    } else {
      prev_st = paste0('Root_', i)
    }
    
    o = edge_tb$to_type[edge_tb$to == to]
    
    # i = get_type_from_id(from)
    # o = get_type_from_id(to)
    
    out_types = edge_tb$to_type[edge_tb$from == from]
    
    if (i == o & !any(out_types != o)) {
      edge_tb$st[edge_tb$from == from & edge_tb$to == to] = prev_st
    } else {
      edge_tb$st[edge_tb$from == from & edge_tb$to == to] = paste0(i, '_', o)
    }
    if (to %in% edge_tb$from) {
      edge_tb = assign_state_transition(edge_tb, to)
    }
  }
  edge_tb
}

#edge_tb$block = map_chr(1:nrow(edge_tb), function(n) 'unassigned')
assign_block <- function(edge_tb, from, root = F) {
  
  i = edge_tb$from_type[edge_tb$from == from][1]
  
  j <<- j + 1
  if (j %% 1000 == 0) {
    message(j)
  }
  out_cells = edge_tb$to[edge_tb$from == from]
  for (to in out_cells) {
    if (!root) {
      prev_block = edge_tb$block[edge_tb$to == from]
      prev_st = edge_tb$st[edge_tb$to == from]
    } else {
      prev_block = paste0(sample(c(letters, 1:9), size = 10, replace = T), collapse = "")
      prev_st = paste0('Root_', i)
    }
    
    curr_st = edge_tb$st[edge_tb$to == to]
    
    if (curr_st == prev_st) {
      edge_tb$block[edge_tb$to == to] = prev_block
    } else {
      edge_tb$block[edge_tb$to == to] = paste0(sample(c(letters, 1:9), size = 10, replace = T), collapse = "")
    }
    
    if (to %in% edge_tb$from) {
      edge_tb = assign_block(edge_tb, to)
    }
  }
  edge_tb
}

#edge_tb$from_pseudo = map_dbl(1:nrow(edge_tb), function(n) 0)
#edge_tb$to_pseudo = map_dbl(1:nrow(edge_tb), function(n) 0)
assign_all_pseudotime <- function(edge_tb, edges, tip_list) {
  
  edge_tb$block_copy = edge_tb$block
  edge_tb_nested_blocks = nest(edge_tb, data = -block_copy)
  edge_tb_nested_blocks$data = future_map(1:nrow(edge_tb_nested_blocks), function(n) {
    message(n)
    tb = edge_tb_nested_blocks$data[[n]]
    # j <<- j + 1
    # if (j %% 1 == 0) {
    #   message(j)
    # }
    
    states = strsplit(tb$st[1], '_')[[1]]
    
    if (states[1] == states[2]) {
      tb$from_pseudo = edges$end_time[edges$out_node == states[1]]
      tb$to_pseudo = edges$end_time[edges$out_node == states[2]]
      return(tb)
    }
    
    pseudo_start = edges$start_time[edges$in_node == states[1]][1]
    pseudo_end = edges$end_time[edges$out_node == states[2]]
    
    block_start_time = min(tb$from_time)
    # tb$levels = map_dbl(1:nrow(tb), function(n) 0)
    # tb = assign_levels(tb, block_start_time, 1)
    
    tb = assign_pseudotime(tb, tb$from[tb$from_time == block_start_time], pseudo_start, pseudo_end, tip_list, start = T)
    
    # first_edge = tb[tb$from_time == block_start_time,]
    # while (edges$in_node[edges$out_node == states[2]] != states[1]) {
    #   last_branch = edges[edges$out_node == states[2],]
    #   if (first_edge$to_pseudo >= last_branch$start_time) {
    #     intermediate_time = last_branch$start_time
    #   } else {
    #     intermediate_time = first_edge$to_pseudo - 0.01
    #   }
    #   
    #   new_edge_from_time = first_edge$length*(intermediate_time - first_edge$from_pseudo)/
    #     (first_edge$to_pseudo - first_edge$from_pseudo) + first_edge$from_time
    #   
    #   if (nrow(tb) > 1) {
    #     new_block = first_edge$block
    #   } else {
    #     new_block = paste0(sample(c(letters, 1:9), size = 10, replace = T), collapse = "")
    #   }
    #   
    #   new_edge = list(paste0(first_edge$to, 'prev'),
    #                first_edge$to,
    #                first_edge$to_time - new_edge_from_time,
    #                new_edge_from_time,
    #                first_edge$to_time,
    #                last_branch$in_node,
    #                first_edge$to_type,
    #                paste0(last_branch$in_node, '_', first_edge$to_type),
    #                new_block,
    #                intermediate_time,
    #                first_edge$to_pseudo
    #                )
    #   names(new_edge) = names(tb)
    #   
    #   first_edge = list(first_edge$from,
    #                  new_edge$from,
    #                  new_edge$from_time - first_edge$from_time,
    #                  first_edge$from_time,
    #                  new_edge$from_time,
    #                  first_edge$from_type,
    #                  new_edge$from_type,
    #                  paste0(first_edge$from_type, '_', new_edge$from_type),
    #                  first_edge$block,
    #                  first_edge$from_pseudo,
    #                  new_edge$from_pseudo
    #   )
    #   
    #   names(first_edge) = names(tb)
    #   
    #   tb = rbind(tb, new_edge)
    #   
    #   states = strsplit(first_edge$st, '_')[[1]]
    # }
    # 
    # tb = tb[tb$from_time != block_start_time,]
    # tb$st = paste0(tb$from_type, '_', tb$to_type)
    # first_edge$block = paste0(sample(c(letters, 1:9), size = 10, replace = T), collapse = "")
    # tb = rbind(tb, first_edge)
    
    tb
  }, .progress = T)
  
  edge_tb = unnest(edge_tb_nested_blocks, cols = data)
  edge_tb = edge_tb %>% select(-block_copy)
  edge_tb
}

assign_levels <- function(tb, from_times, level) {
  ind = tb$from_time %in% from_times
  tb$levels[ind] = level
  if (any(tb$levels == 0)) {
    tb = assign_levels(tb, tb$to_time[ind],level+1)
  }
  tb
}

assign_tips <- function(tb, to_list, tips = F) {
  if (tips) {
    tb$tips[tb$to %in% to_list] = tb$to[tb$to %in% to_list]
  }
  for (to in to_list) {
    if (to %in% tb$to) {
      from = tb$from[tb$to == to]
      tb$tips[tb$from == from] = c(tb$tips[tb$from == from], tb$tips[tb$to == to])
      tb = assign_tips(tb, from)
    }
  }
  tb
}

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

assign_pseudotime <- function(tb, from, pseudo_start, pseudo_end, tip_list, start = F) {
  
  # k <<- k+1
  # if (k %% 1 == 0) {
  #   message(k)
  # }
  
  out_cells = tb$to[tb$from == from]
  tips = list_dd_and_tips_from_edge_tb(tb)$tips[[from]]
  
  if (is.null(tips)) {
    tips = tb$to[!(tb$to %in% tb$from)]
  }
  
  avg_time = mean(tb$to_time[tb$to %in% tips] - tb$from_time[tb$from == from][1])
  # avg_time = sum((tb[tb$levels >= level,] %>%
  #                   group_by(levels) %>%
  #                   summarise(avg_time = mean(length)))$avg_time)
  
  if (all(out_cells %in% tip_list)) {
    
    time_so_far = tb$from_time[tb$from == from][1] - min(tb$from_time)
    b = 0.75*(pseudo_end - pseudo_start + 0.5)
    lambda = 4 / (pseudo_end - pseudo_start + 0.5)
    s = pseudo_end - pseudo_start
    
    if (!start) {
      prev_pseudo_end = tb$to_pseudo[tb$to == from]
    } else {
      prev_pseudo_end = pseudo_start
    }
    
    tb$from_pseudo[tb$to %in% out_cells] = prev_pseudo_end
    
    t = rexp(1, 1/lambda) + b
    
    if (time_so_far == 0 && length(out_cells) == 1) {
      to = out_cells[[1]]
      if (t < tb$length[tb$to == to]) {
        tb$to_pseudo[tb$to == to] = pseudo_end
      } else {
        pseudo_length = pseudo_end - prev_pseudo_end
        tb$to_pseudo[tb$to == to] = tb$length[tb$to == to] * pseudo_length / t + prev_pseudo_end
      }
      return(tb)
    }
    
    if (t < time_so_far) {
      tb$to_pseudo[tb$to %in% out_cells] = pseudo_end
    } else {
      for (to in out_cells) {
        t = rexp(1, 1/lambda) + max(b - time_so_far)
        if (t < tb$length[tb$to == to]) {
          tb$to_pseudo[tb$to == to] = pseudo_end
        } else {
          pseudo_length = pseudo_end - prev_pseudo_end
          tb$to_pseudo[tb$to == to] = (time_so_far + tb$length[tb$to == to]) *
            pseudo_length / (time_so_far + t) + prev_pseudo_end
        }
      }
    }
    
    return(tb)
    
  }
  
  for (to in out_cells) {
    
    if (!start) {
      prev_pseudo_end = tb$to_pseudo[tb$to == from]
    } else {
      prev_pseudo_end = pseudo_start
    }
    
    tb$from_pseudo[tb$to == to] = prev_pseudo_end
    
    if (to %in% tb$from) {
      pseudo_length = pseudo_end - prev_pseudo_end
      tb$to_pseudo[tb$to == to] = prev_pseudo_end + tb$length[tb$to == to] * pseudo_length / avg_time
      tb = assign_pseudotime(tb, to, pseudo_start, pseudo_end, tip_list)
    } else {
      tb$to_pseudo[tb$to == to] = pseudo_end
    }
    
  }
  
  tb
}

assign_program_value <- function(edges, edge, program_name) {
  states = strsplit(edge$st, '_')[[1]]
  
  if (states[1] == states[2] || edge$to_pseudo - edge$from_pseudo <= 0.01) {
    branch = edges[edges$out_node == states[2],]
    return(tibble(time = c(edge$from_time, edge$to_time),
                  value = branch[[program_name]][[1]][length(branch[[program_name]][[1]])]))
  }
  
  branch = edges[edges$in_node == states[1] & edges$out_node == states[2],]
  
  # if (nrow(branch) == 0) {
  #   
  #   last_branch = edges[edges$out_node == states[2],]
  #   branches
  #   
  #   
  #   
  #   branches = edges[edges$in_node == states[1] | edges$out_node == states[2],]
  #   branches = branches[branches$out_node %in% branches$in_node |
  #                         branches$in_node %in% branches$out_node,]
  #   branch1 = branches[branches$out_node %in% branches$in_node,]
  #   intermediate_time = branch1$end_time
  #   
  #   steps1 = length(branch1[[program_name]][[1]])
  #   start_val1 = round((edge$from_pseudo - branch1$start_time) / branch1$edge_length * steps1) + 1
  #   end_val1 = round((intermediate_time - branch1$start_time) / branch1$edge_length * steps1)
  #   
  #   branch2 = branches[branches$in_node %in% branches$out_node,]
  #   steps2 = length(branch2[[program_name]][[1]])
  #   start_val2 = round((intermediate_time - branch2$start_time) / branch2$edge_length * steps2) + 1
  #   end_val2 = round((edge$to_pseudo - branch2$start_time) / branch2$edge_length * steps2) + 1
  #   
  #   time = seq(edge$from_time, edge$to_time,
  #              (edge$to_time - edge$from_time)/(end_val2 - start_val2 + end_val1 - start_val1 + 1))
  #   value = c(branch1[[program_name]][[1]][start_val1:end_val1],
  #             branch2[[program_name]][[1]][start_val2:end_val2])
  #   return(tibble(time = time, value = value))
  # }
  
  steps = length(branch[[program_name]][[1]])
  start_val = round((edge$from_pseudo - branch$start_time) / branch$edge_length * steps) + 1
  end_val = round((edge$to_pseudo - branch$start_time) / branch$edge_length * steps)
  
  tibble(time = seq(edge$from_time, edge$to_time, (edge$to_time - edge$from_time)/(end_val - start_val)),
         value = branch[[program_name]][[1]][start_val:end_val])
}

assign_all_program_value <- function(edges, edge_tb, program_name) {
  edge_tb[[program_name]] = map(1:nrow(edge_tb), function(n) {
    if (n %% 10000 == 0) {
      message(n)
    }
    edge = edge_tb[n,]
    assign_program_value(edges, edge, program_name)
  })
  edge_tb
}

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

make_mut <- function() {
  paste0(sample(c(letters, 1:9), size = 20, replace = T), collapse = "")
}

library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('/home/ssrikan2/data-kreza1/smriti/qfm2')
devtools::load_all()

setwd('/home/ssrikan2/data-kreza1/smriti/MF_Signal_Simulation/lineage specific simulation/')

load('./3_13_edge_tb.rda')
load('./3_13_edges.rda')


library(furrr)
plan(multisession, workers = 8)

tip_list = edge_tb$to[!(edge_tb$to %in% edge_tb$from)]

edge_tb$st = map_chr(1:nrow(edge_tb), function(n) 'unassigned')
j = 0
edge_tb = assign_state_transition(edge_tb, edge_tb$from[[1]], T)
#edge_tb$st = paste0(edge_tb$from_type, '_', edge_tb$to_type)
edge_tb$block = map_chr(1:nrow(edge_tb), function(n) 'unassigned')
j = 0
edge_tb = assign_block(edge_tb, edge_tb$from[[1]], T)
edge_tb$from_pseudo = map_dbl(1:nrow(edge_tb), function(n) 0)
edge_tb$to_pseudo = map_dbl(1:nrow(edge_tb), function(n) 0)
j = 0
st = proc.time()
edge_tb = assign_all_pseudotime(edge_tb, edges, tip_list)
proc.time() - st
edge_tb = assign_all_program_value(edges, edge_tb, 'program3')







tr_dd = split(edge_tb$to, edge_tb$from)
tr_node_time = edge_tb$length
tr_node_start = edge_tb$from_time
tr_node_end = edge_tb$to_time
names(tr_node_time) = names(tr_node_start) = names(tr_node_end) = edge_tb$to

tr_tips = list_dd_and_tips_from_edge_tb(edge_tb)$tips

program_val = edges %>% select(program3) %>% unnest(cols = program3) %>% pull(program3)
median_program = median(program_val)
base_mean = log(10^(1)) - median_program
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

save(cell_mut_tb, file = './3_13_cell_mut_tb.rda')
save(edge_tb, file = './3_13_edge_tb_with_pseudotime_and_program.rda')


