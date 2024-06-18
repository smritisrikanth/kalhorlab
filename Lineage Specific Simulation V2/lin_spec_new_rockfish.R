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

break_up_state_jumps <- function(edge_tb) {
  
  correct_from_type = map_chr(edge_tb$to_type, function(to_type) {
    edges$in_node[edges$out_node == to_type]
  })
  
  jumps = edge_tb[correct_from_type != edge_tb$from_type &
                    edge_tb$from_type != edge_tb$to_type,]
  
  if (nrow(jumps) == 0) {
    return(edge_tb)
  }
  
  new_tb_list = map(1:nrow(jumps), function(n) {
    jump = jumps[n,]
    
    tb = NULL
    
    while(edges$in_node[edges$out_node == jump$to_type] != jump$from_type) {
      last_branch = edges[edges$out_node == jump$to_type,]
      total_pseudotime = edges$end_time[edges$out_node == jump$to_type] - 
        edges$start_time[edges$in_node == jump$from_type][1]
      frac_pseudotime = (last_branch$end_time - last_branch$start_time)/total_pseudotime
      new_from_time = jump$to_time - frac_pseudotime * jump$length
      
      new_edge = jump
      
      new_edge$from = paste0(jump$to, 'prev')
      new_edge$to = jump$to
      new_edge$from_type = last_branch$in_node
      new_edge$to_type = jump$to_type
      new_edge$from_time = new_from_time
      new_edge$to_time = jump$to_time
      new_edge$length = new_edge$to_time - new_edge$from_time
      
      jump$from = jump$from
      jump$to = new_edge$from
      jump$from_type = jump$from_type
      jump$to_type = new_edge$from_type
      jump$from_time = jump$from_time
      jump$to_time = new_edge$from_time
      jump$length = jump$to_time - jump$from_time
      
      tb = rbind(tb, new_edge)
      
    }
    
    tb = rbind(jump,tb)
  })
  
  new_tb = bind_rows(new_tb_list)
  
  edge_tb = edge_tb[!(edge_tb$to %in% jumps$to),]
  
  return(rbind(edge_tb, new_tb))
  
}


assign_edge_program_value <- function(edge_tb, from, v0 = rnorm(1,0,0.2), program_name, root = F) {
  
  j <<- j+1
  if (j %% 1 == 0) {
    message(j)
  }
  
  out_cells = edge_tb$to[edge_tb$from == from]
  
  if (!root) {
    parent_edge = edge_tb[edge_tb$to == from,]
    start_val = parent_edge[[program_name]][[1]]$value[nrow(parent_edge[[program_name]][[1]])]
  } else {
    parent_edge = list(block = 'None')
    start_val = log(runif(1,0,1.5))
  }
  
  
  
  for (to in out_cells) {
    
    edge = edge_tb[edge_tb$to == to,]
    
    if (root) {
      
      t_stable = edge$length * edge_proportion / step_size
      eta_hat = uniroot(eta_func, c(0,0.99), t_stable, delta)
      eta = eta_hat$root
      v0 = rnorm(1,0,0.2)
    } else if (edge$block != parent_edge$block) {
      states = strsplit(edge$st, '_')[[1]]
      eta = edges$eta[edges$in_node == states[1] & edges$out_node == states[2]]
      v0 = edges$v0[edges$in_node == states[1] & edges$out_node == states[2]][[1]][[program_name]]
    }
    
    num_steps = round((edge$to_time - edge$from_time) / step_size)
    
    # walk = c(start_val,rep(0,num_steps))
    # error = c(0,rnorm(num_steps,0,step_size))
    # velocity = c(0,map_dbl(0:(num_steps-1), function(i) eta^i)) * v0 + error
    # walk = cumsum(walk + velocity)
    
    error = rnorm(num_steps-1,0,step_size)
    coefs = cumsum(map_dbl(0:(num_steps-1), function(i) eta^i))
    error_terms = c(0,error * coefs[1:length(error)])
    walk = c(start_val,coefs*v0 + error_terms + start_val)
    v_end = walk[length(walk)] - walk[length(walk) - 1]
    
    edge_tb[[program_name]][edge_tb$to == to] = list(tibble(time = linspace(edge$from_time, edge$to_time,num_steps+1),
                                                            value = walk))
    
    if (to %in% edge_tb$from) {
      edge_tb = assign_edge_program_value(edge_tb, to, v_end,program_name)
    }
  }
  
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

library(uwot)
library(ComplexHeatmap)
library(qfm)
library(pracma)
library(readxl)
setwd('/home/ssrikan2/data-kreza1/smriti/qfm2')
devtools::load_all()
library(furrr)
plan(multisession, workers = 8)

setwd('/home/ssrikan2/data-kreza1/smriti/MF_Signal_Simulation/lineage specific simulation/')
load('gene_expr_data_v0.rda')
edge_length_tb = phy_edges[-3]
colnames(edge_length_tb) = c('in_node', 'out_node', 'start_time', 'end_time')
edges = edge_length_tb
# edge_tb = complete_edges
# colnames(edge_tb) = c('from', 'to', 'from_type', 'to_type', 'from_time', 'to_time', 'length')
# 
# # #create edge topology tb
# # edge_length_tb = readRDS('./phy_com_edges.rds')[,-3]
# # #edge_length_tb = readRDS('./phy_edges_len.rds')
# # colnames(edge_length_tb) = c('in_node', 'out_node', 'start_time', 'end_time')
# # edge_length_tb$edge_length = edge_length_tb$end_time - edge_length_tb$start_time
# # root_edge = list('Root', edge_length_tb$in_node[[1]], 0, edge_length_tb$start_time[[1]],
# #                  edge_length_tb$start_time[[1]])
# # edge_length_tb = rbind(root_edge, edge_length_tb)
# # 
# # edges = edge_length_tb
# # 
# # #load in edge phylogeny tb
# # edge_tb = readRDS('./tr_edges_com.rds')
# # colnames(edge_tb) = c('from', 'to', 'from_time', 'to_time', 'length', 'from_type', 'to_type')
# # ind = map_lgl(edge_tb$from, function(from) {
# #   length(unique(edge_tb$to_type[edge_tb$from == from])) != 1
# # })
# # edge_tb$from_type[ind] = "Foregut1"
# # edge_tb$to_type[ind] = "Foregut1"
# 
# #check
# sum(map_lgl(edge_tb$to_type, function(to_type) {
#   !(to_type %in% edges$out_node)
# }))
# sum(map_lgl(edge_tb$from_type, function(from_type) {
#   !(from_type %in% edges$in_node | from_type %in% edges$out_node)
# }))
# 
# #break up state jumps
# edge_tb = break_up_state_jumps(edge_tb)
# 
# #check
# sum(map2_lgl(edge_tb$to_type, edge_tb$from_type, function(to_type, from_type) {
#   edges$in_node[edges$out_node == to_type] != from_type &
#     to_type != from_type
# }))
# 
# #assign blocks
# tip_list = edge_tb$to[!(edge_tb$to %in% edge_tb$from)]
# 
# edge_tb$st = map_chr(1:nrow(edge_tb), function(n) 'unassigned')
# j = 0
# edge_tb = assign_state_transition_v2(edge_tb, edge_tb$from[[1]], T)
# edge_tb$block = map_chr(1:nrow(edge_tb), function(n) 'unassigned')
# j = 0
# edge_tb = assign_block(edge_tb, edge_tb$from[[1]], T)
# 
# 

load('./6_10_edge_tb_before_programs.rda')


#set eta
t_total = edge_tb %>% nest(data = -c(block,st)) %>% mutate(block_length = map_dbl(data, function(tb) {
  tips = list_dd_and_tips_from_edge_tb(tb)$tips[[tb$from[1]]]

  if (is.null(tips)) {
    tips = tb$to[!(tb$to %in% tb$from)]
  }

  avg_time = mean(tb$to_time[tb$to %in% tips] - tb$from_time[1])
  avg_time
})) %>% group_by(st) %>% summarise(avg_block_length = mean(block_length))

step_size = 1/100
delta = 0.3^60 * (1-0.3)
edge_proportion = 0.5

edge_length_tb$eta = map_dbl(1:nrow(edge_length_tb), function(n) {
  st = paste0(edge_length_tb$in_node[n], '_', edge_length_tb$out_node[n])
  states = strsplit(st, '_')[[1]]
  l = t_total$avg_block_length[t_total$st == st]
  t_stable = l * edge_proportion / step_size
  eta_hat = uniroot(eta_func, c(0,0.99), t_stable, delta)
  eta_hat$root
})

#set v0
num_program = 2
edge_length_tb$v0 = map(1:nrow(edge_length_tb), function(n) {
  v0_vec = rnorm(num_program,0,0.2)
  names(v0_vec) = map_chr(1:num_program, function(i) paste0('program', i))
  v0_vec
})

edges = edge_length_tb
# 
# save(edge_tb, file = './6_10_edge_tb_before_programs.rda')

#simulate gene programs
for (i in 1:num_program) {
  
  j = 0
  program_name = paste0('program', i)
  edge_tb[[program_name]] = map(1:nrow(edge_tb), function(n) NULL)
  edge_tb = assign_edge_program_value(edge_tb, edge_tb$from[[1]],
                                      program_name = program_name,
                                      root = T)
}

save(edge_tb, file = './6_10_edge_tb.rda')



