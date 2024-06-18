library(uwot)
library(ComplexHeatmap)
library(qfm)
library(pracma)
setwd('~/Documents/qfm2/')
devtools::load_all()
library(furrr)
plan(multisession, workers = 8)

#functions

eta_func <- function(eta, t, delta) {
  eta^t * (1-eta) - delta
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

assign_all_pseudotime <- function(edge_tb, edges, tip_list) {
  
  edge_tb$block_copy = edge_tb$block
  edge_tb_nested_blocks = nest(edge_tb, data = -block_copy)
  edge_tb_nested_blocks$data = future_map(1:nrow(edge_tb_nested_blocks), function(n) {
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
    
    tb = break_up_state_jumps(tb)
    
    tb
  }, .progress = T)
  
  edge_tb = unnest(edge_tb_nested_blocks, cols = data)
  edge_tb = edge_tb %>% select(-block_copy)
  edge_tb
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
      #tb$to_pseudo[tb$to == to] = pseudo_end #MAKE STOCHASTIC
      
      time_so_far = tb$from_time[tb$from == from][1] - min(tb$from_time)
      b = 0.75*(pseudo_end - pseudo_start + 0.5)
      lambda = 4 / (pseudo_end - pseudo_start + 0.5)
      s = pseudo_end - pseudo_start
      
      t = rexp(1, 1/lambda) + b
      
      if (time_so_far == 0 && length(out_cells) == 1) {
        to = out_cells[[1]]
        if (t < tb$length[tb$to == to]) {
          tb$to_pseudo[tb$to == to] = pseudo_end
        } else {
          pseudo_length = pseudo_end - prev_pseudo_end
          tb$to_pseudo[tb$to == to] = tb$length[tb$to == to] * pseudo_length / t + prev_pseudo_end
        }
      } else if (t < time_so_far) {
        tb$to_pseudo[tb$to %in% out_cells] = pseudo_end
      } else {
        t = rexp(1, 1/lambda) + max(b,time_so_far)
        if (t < tb$length[tb$to == to]) {
          tb$to_pseudo[tb$to == to] = pseudo_end
        } else {
          pseudo_length = pseudo_end - prev_pseudo_end
          tb$to_pseudo[tb$to == to] = (time_so_far + tb$length[tb$to == to]) *
            pseudo_length / (time_so_far + t) + prev_pseudo_end
        }
      }
    }
    
  }
  
  tb
}

break_up_state_jumps <- function(tb) {
  first_edge = tb[tb$from_time == min(tb$from_time),]
  
  if (first_edge$from_type == first_edge$to_type) {
    return(tb)
  }
  
  while (edges$in_node[edges$out_node == first_edge$to_type] != first_edge$from_type) {
    print(tb)
    last_branch = edges[edges$out_node == first_edge$to_type,]
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
  }
  
  tb = tb[tb$from_time != min(tb$from_time),]
  tb$st = paste0(tb$from_type, '_', tb$to_type)
  first_edge$block = paste0(sample(c(letters, 1:9), size = 10, replace = T), collapse = "")
  tb = rbind(tb, first_edge)
  
  tb
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

setwd('~/Documents/R/MF_Signal_Simulation/lineage specific simulation/')

#create edge topology tb
#edge_length_tb = readRDS('./phy_com_edges.rds')[,-3]
edge_length_tb = readRDS('./phy_edges_len.rds')
colnames(edge_length_tb) = c('in_node', 'out_node', 'start_time', 'end_time')
edge_length_tb$edge_length = edge_length_tb$end_time - edge_length_tb$start_time
root_edge = list('Root', edge_length_tb$in_node[[1]], 0, edge_length_tb$start_time[[1]],
                 edge_length_tb$start_time[[1]])
edge_length_tb = rbind(root_edge, edge_length_tb)

step_size = 1/100
delta = 0.3^80 * (1-0.3)
edge_proportion = 0.5

edge_length_tb$eta = map_dbl(edge_length_tb$edge_length, function(l) {
  t = l * edge_proportion / step_size
  ep1 = eta_func(0,t,delta)
  ep2 = eta_func(0.99,t,delta)
  # if (ep1 < 0 && ep2 < 0 || ep1 > 0 && ep2 > 0) {
  #   return(0)
  # }
  eta_hat = uniroot(eta_func, c(0,0.99), t, delta)
  eta_hat$root
})
median(edge_length_tb$eta)

edge_length_tb$v0 = rnorm(nrow(edge_length_tb),0,0.2)

edges = edge_length_tb

#load edge_tb, assign st and blocks, and map pseudotime
load(paste0("~/Documents/R/MF_Signal_Simulation/data_1.rda"))
tr_samp = subsample_tr(tr, ss_val = 100)
edge_tb = make_edge_tb(tr_samp)
edge_tb$from_type = get_type_from_id(edge_tb$from)
edge_tb$to_type = get_type_from_id(edge_tb$to)

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
edge_tb$program1 = map(1:nrow(edge_tb), function(n) NULL)
j = 0
edge_tb = assign_edge_program_value(edge_tb, edge_tb$from[[1]])


# assign_all_program_value_v2 <- function(edge_tb, edges) {
#   
#   edge_tb_nested_blocks = nest(edge_tb,data = -block)
#   edge_tb_nested_blocks$data = map(edge_tb_nested_blocks$data, function(tb) {
#     tb = assign_block_program_value(tb)
#   })
#   
#   edge_tb = unnest(edge_tb_nested_blocks, cols = data)
#   edge_tb
#      
# }

# assign_block_program_value <- function(tb) {
#   
#   states = strsplit(tb$st[1], '_')[[1]]
#   eta = edges$eta[edges$in_node == states[1] && edges$out_node == states[2]]
#   v0 = edges$v0[edges$in_node == states[1] && edges$out_node == states[2]]
#   
#   block_root = tb$from[tb$from_time == min(tb$from_time)][1]
#   tips = tb$to[!(tb$to %in% tb$from)]
#   
#   #should program be generated along pseudotime or real time??
#   # avg_pseudo_time = mean(tb$to_pseudo[tb$to %in% tips] - tb$from_time[tb$from == block_root][1])
#   # total_steps = round(avg_time / step_size)
#   tb[['program1']] = map(1:nrow(tb), function(n) NULL)
#   
#   tb = assign_edge_program_value(block_root)
#   
#   
# }

assign_edge_program_value <- function(edge_tb, from, v0 = rnorm(1,0,0.2)) {
  
  j <<- j+1
  if (j %% 1 == 0) {
    message(j)
  }
  
  out_cells = edge_tb$to[edge_tb$from == from]
  
  if (edge_tb$from_time[edge_tb$from == from][1] > 0) {
    parent_edge = edge_tb[edge_tb$to == from,]
    start_val = parent_edge[['program1']][[1]]$value[nrow(parent_edge[['program1']][[1]])]
  } else {
    parent_edge = list(block = 'None')
    start_val = log(runif(1,0,1.5))
  }
  
  
  
  for (to in out_cells) {
    
    if (j == 102) {
      message(to)
    }
    
    edge = edge_tb[edge_tb$to == to,]
    
    if (edge$block != parent_edge$block) {
      states = strsplit(edge$st, '_')[[1]]
      eta = edges$eta[edges$in_node == states[1] & edges$out_node == states[2]]
      v0 = edges$v0[edges$in_node == states[1] & edges$out_node == states[2]]
    }
    
    num_steps = round((edge$to_pseudo - edge$from_pseudo) / step_size)
    
    walk = c(start_val,rep(0,num_steps))
    error = c(0,rnorm(num_steps,0,step_size))
    velocity = c(0,map_dbl(0:(num_steps-1), function(i) eta^i)) * v0 + error
    walk = cumsum(walk + velocity)
    
    edge_tb[['program1']][edge_tb$to == to] = list(tibble(time = linspace(edge$from_time, edge$to_time,num_steps+1),
                                                         value = walk))
    
    if (to %in% edge_tb$from) {
      edge_tb = assign_edge_program_value(edge_tb, to, velocity[length(velocity)])
    }
  }
  
  edge_tb
  
}
