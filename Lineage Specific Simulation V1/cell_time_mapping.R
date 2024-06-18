
library(furrr)
plan(multisession, workers = 8)

library(ggplot2)
set.seed(0)

# edge_tb$st = map_chr(1:nrow(edge_tb), function(n) 'unassigned')
i <- 0
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

assign_state_transition_v2 <- function(edge_tb, from, root = F) {
  
  j <<- j+1
  if (j %% 1000 == 0) {
    message(j)
  }
  
  i = edge_tb$from_type[edge_tb$from == from][1]
  
  out_cells = edge_tb$to[edge_tb$from == from]
  
  if (!root) {
    in_type = edge_tb$from_type[edge_tb$to == from]
  } else {
    in_type = 'Root'
  }
  
  for (to in out_cells) {
    
    o = edge_tb$to_type[edge_tb$to == to]
    
    if (!root) {
      prev_st = edge_tb$st[edge_tb$to == from]
    } else {
      prev_st = paste0('Root_', i)
    }
    
    if (in_type == i) {
      edge_tb$st[edge_tb$from == from & edge_tb$to == to] = prev_st
    } else {
      edge_tb$st[edge_tb$from == from & edge_tb$to == to] = paste0(in_type, '_', i)
    }
    
    if (to %in% edge_tb$from) {
      edge_tb = assign_state_transition_v2(edge_tb, to)
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
    
    first_edge = tb[tb$from_time == block_start_time,]
    while (edges$in_node[edges$out_node == states[2]] != states[1]) {
      last_branch = edges[edges$out_node == states[2],]
      if (first_edge$to_pseudo >= last_branch$start_time) {
        intermediate_time = last_branch$start_time
      } else {
        intermediate_time = first_edge$to_pseudo - 0.01
      }

      new_edge_from_time = first_edge$length*(intermediate_time - first_edge$from_pseudo)/
        (first_edge$to_pseudo - first_edge$from_pseudo) + first_edge$from_time

      if (nrow(tb) > 1) {
        new_block = first_edge$block
      } else {
        new_block = paste0(sample(c(letters, 1:9), size = 10, replace = T), collapse = "")
      }
      
      new_edge = first_edge
      
      new_edge$from = paste0(first_edge$to, 'prev')
      new_edge$to = first_edge$to
      new_edge$from_type = last_branch$in_node
      new_edge$to_type = first_edge$to_type
      new_edge$from_time = new_edge_from_time
      new_edge$to_time = first_edge$to_time
      new_edge$length = first_edge$to_time - new_edge_from_time
      new_edge$st = paste0(last_branch$in_node, '_', first_edge$to_type)
      new_edge$block = new_block
      new_edge$from_pseudo = intermediate_time
      new_edge$to_pseudo = first_edge$to_pseudo
      
      first_edge$from = first_edge$from
      first_edge$to = new_edge$from
      first_edge$from_type = first_edge$from_type
      first_edge$to_type = new_edge$from_type
      first_edge$from_time = first_edge$from_time
      first_edge$to_time = new_edge$from_time
      first_edge$length = new_edge$from_time - first_edge$from_time
      first_edge$st = paste0(first_edge$from_type, '_', new_edge$from_type)
      first_edge$block = first_edge$block
      first_edge$from_pseudo = first_edge$from_pseudo
      first_edge$to_pseudo = new_edge$from_pseudo

      tb = rbind(tb, new_edge)

      states = strsplit(first_edge$st, '_')[[1]]
    }

    tb = tb[tb$from_time != block_start_time,]
    tb$st = paste0(tb$from_type, '_', tb$to_type)
    first_edge$block = paste0(sample(c(letters, 1:9), size = 10, replace = T), collapse = "")
    tb = rbind(tb, first_edge)
    
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
  
  if (nrow(edge_df) == 1) {
    tips = rep(edge_df$to[!edge_df$to %in% edge_df$from], 2)
    names(tips) = c(edge_df$from[1], edge_df$to[1])
    downs = list(c(edge_df$from[1], edge_df$to[1]), edge_df$to[1])
    names(downs) = c(edge_df$from[1], edge_df$to[1])
    tr_dd_new = split(edge_df$to, edge_df$from)
    res = list(tips = tips,
         downs = downs,
         dd = tr_dd_new)
    return(res)
  }
  
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
  
  if (nrow(branch) == 0) {

    last_branch = edges[edges$out_node == states[2],]
    branches



    branches = edges[edges$in_node == states[1] | edges$out_node == states[2],]
    branches = branches[branches$out_node %in% branches$in_node |
                          branches$in_node %in% branches$out_node,]
    branch1 = branches[branches$out_node %in% branches$in_node,]
    intermediate_time = branch1$end_time

    steps1 = length(branch1[[program_name]][[1]])
    start_val1 = round((edge$from_pseudo - branch1$start_time) / branch1$edge_length * steps1) + 1
    end_val1 = round((intermediate_time - branch1$start_time) / branch1$edge_length * steps1)

    branch2 = branches[branches$in_node %in% branches$out_node,]
    steps2 = length(branch2[[program_name]][[1]])
    start_val2 = round((intermediate_time - branch2$start_time) / branch2$edge_length * steps2) + 1
    end_val2 = round((edge$to_pseudo - branch2$start_time) / branch2$edge_length * steps2) + 1

    time = seq(edge$from_time, edge$to_time,
               (edge$to_time - edge$from_time)/(end_val2 - start_val2 + end_val1 - start_val1 + 1))
    value = c(branch1[[program_name]][[1]][start_val1:end_val1],
              branch2[[program_name]][[1]][start_val2:end_val2])
    return(tibble(time = time, value = value))
  }
  
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


subsample_tr <- function(tr, ss_val = 1000) {
  all_tip_tb = tibble(tip = tr$tip.label,
                      tip_type = get_type_from_id(tr$tip.label))
  sampled_tip_tb = all_tip_tb %>% nest(tips = -tip_type) %>%
    mutate(tips_sample = map(tips, function(tb) {
      ss_tip = pmin(ss_val, nrow(tb))
      tibble(tip = sample(tb$tip, size = ss_tip, replace = F))
    })) %>%
    select(-tips) %>%
    unnest(cols = "tips_sample")
  tr_samp = drop.tip(tr, tr$tip.label[!tr$tip.label %in% sampled_tip_tb$tip])
  tr_samp
}


load(paste0("~/Documents/R/MF_Signal_Simulation/data_1.rda"))
tr_samp = subsample_tr(tr, ss_val = 200)
edge_tb = make_edge_tb(tr_samp)
edge_tb$from_type = get_type_from_id(edge_tb$from)
edge_tb$to_type = get_type_from_id(edge_tb$to)
phy = readRDS("~/Documents/R/MF_Signal_Simulation/gast_phylo_1.rds")

phy = readRDS('./phy_com_scaled.rds')
edge_tb = readRDS('./tr_edges_com.rds')
colnames(edge_tb) = c('from', 'to', 'from_time', 'to_time', 'length', 'from_type', 'to_type')
ind = map_lgl(edge_tb$from, function(from) {
  length(unique(edge_tb$to_type[edge_tb$from == from])) != 1
})
edge_tb$from_type[ind] = "Foregut1"
edge_tb$to_type[ind] = "Foregut1"

tip_list = edge_tb$to[!(edge_tb$to %in% edge_tb$from)]

edge_tb$st = map_chr(1:nrow(edge_tb), function(n) 'unassigned')
j = 0
edge_tb = assign_state_transition_v2(edge_tb, edge_tb$from[[1]], T)
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

ggplot(data = edge_tb, aes(x = from_time,
                           y = from_pseudo,
                           color = from_type)) +
  geom_point()

#check
sum(map_lgl(edge_tb$to, function(to) {
  if (to %in% edge_tb$from) {
    return(edge_tb$to_time[edge_tb$to == to] != edge_tb$from_time[edge_tb$from == to][1])
  } else {
    return(FALSE)
  }
}))

map_lgl(edge_tb$to, function(to) {
  
})


edge_tb$length_pseudo = edge_tb$to_pseudo - edge_tb$from_pseudo
ggscatter(edge_tb, x = "length", y = "length_pseudo",
          facet.by = "to_type") + 
  geom_smooth(color = "red", se = F, method = "lm") +
  geom_abline(color = "blue", intercept = -0.5)


plot_tb = edge_tb %>% select(to_type, program1) %>% unnest(cols = program1)
ggplot(data = plot_tb,
       aes(x = time,
           y = value,
           color = to_type)) + geom_line()

#visualize
#convert to igraph, assign edge attributes


