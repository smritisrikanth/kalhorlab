library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
#library(randomcoloR)
setwd('~/Documents/qfm2')
devtools::load_all()

cell_state_meta = read_csv("./metadata/gas_states.csv")

set.seed(0)

setwd('~/Documents/R/MF_Signal_Simulation')
phy = readRDS('./gast_phylo.rds')


tree = phy
n_rep = 5
max_n_sample = 3

#functions
compare_lists <- function (l1, l2) {
  if ((length(intersect(l1,l2)) == length(l1)) && (length(intersect(l1,l2)) == length(l2))) {
    return(TRUE)
  }
  return(FALSE)
}

compare_two_nodes <- function(node1, node2) {
  if (length(node1) != length(node2)) {
    FALSE
  }
  count = 0
  for (i in 1:length(node1)) {
    for (j in 1:length(node2)) {
      l1 <- node1[[i]]
      l2 <- node2[[j]]
      
      if (compare_lists(l1,l2)) {
        count <- count + 1
      }
    }
  }
  if (count == length(node1)) {
    TRUE
  } else {
    FALSE
  }
}


compare_two_partitions <- function(p, p_ref) {
  temp_list = c()
  final_list = c()
  for (node1 in p_ref) {
    for (node2 in p) {
      temp_list = c(temp_list, compare_two_nodes(node1, node2))
    }
    final_list = c(final_list, any(temp_list))
    temp_list = c()
  }
  final_list
}

get_pooled_partition <- function(node) {
  pooled = as.vector(c())
  for (l in node) {
    pooled = c(pooled, l)
  }
  pooled
}

get_split <- function(pooled_partition) {
  complement = setdiff(tree$tip.label, pooled_partition)
  split = list(pooled_partition, complement)
}

get_unique_splits <- function(splits) {
  repeat_list = as.vector(c())
  for (i in 1:(length(splits)-1)) {
    for (j in (i+1):length(splits)) {
      if (compare_two_nodes(splits[[i]], splits[[j]])) {
        repeat_list = c(repeat_list, j)
      }
    }
  }
  if (is.null(repeat_list)) {
    return(splits)
  }
  splits[-repeat_list]
}

param_tb = as_tibble(expand.grid(distance_method = c('manhattan', 'jaccard'),
                                 reconstruction_method = c('upgma', 'nj', 'ward2'),
                                 filtering = c(-10:0),
                                 dropout = c(1:10),
                                 rep = c(1:n_rep),
                                 n_sample = max_n_sample))


output_folder = 'mouse_ss1000_hgRNA_flat_dist_mats_v3'


param_tb$matrix = map(1:nrow(param_tb), function(n) {
  print(n)
  distance_method = param_tb$distance_method[n]
  reconstruction_method = param_tb$reconstruction_method[n]
  filtering = param_tb$filtering[n]
  dropout = param_tb$dropout[n]
  n_sample = param_tb$n_sample[n]
  
  dist_mat_list = map(1:n_sample, function(s) {
    # num = (n+cumsum(param_tb$n_sample[1:n])[n])%%100 +s-1
    num = floor(runif(1, min=1, max=101))
    file = paste0('./', output_folder, '/dist_mat_', distance_method, '_', filtering, '_', dropout, '_', num, '.rda')
    if (file.exists(file)) {
      load(file)
    }
    while (!file.exists(file) || dim(as.matrix(dist_mat)) != 14) {
      num = floor(runif(1, min=1, max=101))
      file = paste0('./', output_folder, '/dist_mat_', distance_method, '_', filtering, '_', dropout, '_', num, '.rda')
      if (file.exists(file)) {
        load(file)
      }
    }
    dist_mat
  })
  names = rownames(as.matrix(dist_mat_list[[1]]))
  .dist_mat_list = map(dist_mat_list, function(mat) {
    mat = as.matrix(mat)
    ordered_mat = mat[names, names]
    ordered_mat
  })
  cum_mat = Reduce('+', .dist_mat_list)
  
  avg_mat = cum_mat/length(dist_mat_list)
  avg_mat
})

param_tb$dist_list = map(param_tb$matrix, function(mat) {
  df = as.data.frame(mat) %>%
    rownames_to_column(var = "rows") %>%
    pivot_longer(-rows, names_to = "columns", values_to = "distance")
  df_unique <- df %>%
    mutate(min_val = pmin(rows, columns),
           max_val = pmax(rows, columns)) %>%
    distinct(min_val, max_val, distance, .keep_all = TRUE) %>%
    select(-min_val, -max_val)
  df_unique$truth = 9 - map2_dbl(df_unique$rows, df_unique$columns, function(r,c) {
    mrca = phy$node.label[getMRCA(phy, c(r,c)) - Ntip(phy)]
    (cell_state_meta$commit_start[cell_state_meta$name == mrca] +
      cell_state_meta$commit_end[cell_state_meta$name == mrca])/2
  })
  # df_unique$truth[df_unique$rows == df_unique$columns] = 0
  df_unique = df_unique[df_unique$rows != df_unique$columns, ]
  df_unique
})

param_tb$group = paste0(param_tb$reconstruction_method, param_tb$filtering)
plot_tb = param_tb %>% nest(data= -c(group, n_sample)) %>%
  mutate(comb_df = map(data, function(df) {bind_rows(df$dist_list)})) %>%
  mutate(corr = map_dbl(comb_df, function(df) {cor(df$distance, df$truth, method = 'spearman')}))
  
ggplot(data = plot_tb, aes(x = n_sample,
                            y = corr,
                            group = group,
                            color = group)) +
  geom_line()

param_tb$splits = map(1:nrow(param_tb), function(n) {
  print(n)
  avg_mat = param_tb$matrix[[n]]
  reconstruction_method = param_tb$reconstruction_method[n]
  if (reconstruction_method == 'upgma') {
    average_tree =  ape::as.phylo(hclust(as.dist(avg_mat), method = 'average'))
  } else if (reconstruction_method == 'ward2') {
    average_tree =  ape::as.phylo(hclust(as.dist(avg_mat), method = 'ward.D2'))
  } else {
    average_tree = phytools::midpoint.root(nj(as.dist(avg_mat)))
  }
  average_tree <- reorder(average_tree, 'postorder')
  average_tree = name_nodes(average_tree)
  p1 = get_partitions(average_tree)
  # pooled_partitions1 = map(p1, get_pooled_partition)
  pooled_partitions1 = map(p1, function(x) do.call(c, x))
  splits1 = map(pooled_partitions1, get_split)
  unique_splits1 = get_unique_splits(splits1)
  unique_splits1
})

p2 = get_partitions(phy)
pooled_partitions2 = map(p2, get_pooled_partition)
splits2 = map(pooled_partitions2, get_split)
unique_splits2 = get_unique_splits(splits2)

param_tb$complete = map(param_tb$splits, function(split) {
  tibble(split_name = names(unique_splits2), reconstructed = compare_two_partitions(split, unique_splits2))
})
param_tb$correct = map(param_tb$splits, function(split) {
  tibble(split_name = names(split), in_original = compare_two_partitions(unique_splits2, split))
})

param_tb$percent_complete = map_dbl(param_tb$splits, function(split) {
  mean(compare_two_partitions(split, unique_splits2))
})
param_tb$percent_correct = map_dbl(param_tb$splits, function(split) {
  mean(compare_two_partitions(unique_splits2, split))
})

p_comp_per_node_overall = param_tb %>%
  unnest(cols = c(complete)) %>% select(split_name, reconstructed, reconstruction_method, filtering) %>%
  group_by(split_name, reconstruction_method, filtering) %>%
  summarise(percent_complete = mean(reconstructed))
ggplot(data = p_comp_per_node_overall, aes(x = split_name,
                                           y = percent_complete)) +
  geom_bar(stat = 'identity') + facet_wrap(~ reconstruction_method + filtering)

plot_tb_complete = param_tb %>% group_by(reconstruction_method, filtering, n_sample) %>%
  summarise(avg_percent_complete = mean(percent_complete), se = sd(percent_complete)/sqrt(n_sample))
plot_tb_complete$group = paste0(plot_tb_complete$reconstruction_method, plot_tb_complete$filtering)

complete_plot = ggplot(data = plot_tb_complete) +
  geom_errorbar(aes(x = n_sample, ymin = avg_percent_complete-se,
                                                  ymax = avg_percent_complete+se, group = group, color = group)) +
  geom_line(aes(x = n_sample, y = avg_percent_complete, group = group, color = group))

plot_tb_complete = param_tb %>% nest(data = -c(distance_method,
                                         reconstruction_method,
                                         dropout)) %>%
  mutate(opt_filter = map_dbl(data, function(tb) {
    tb$filtering[which(tb$percent_complete == max(tb$percent_complete))][1]
  })) %>% unnest(cols = data) %>% filter(filtering == opt_filter) %>%
  group_by(distance_method,
                   reconstruction_method,
                   dropout,
              opt_filter) %>%
  summarise(avg_percent_complete = mean(percent_complete))
  
  
plot_tb_complete$group = paste0(plot_tb_complete$reconstruction_method,
                                plot_tb_complete$distance_method)
complete_plot = ggplot(data = plot_tb_complete) +
  geom_line(aes(x = dropout/1000, y = avg_percent_complete, group = group, color = group)) +
  geom_point(aes(x = dropout/1000, y = avg_percent_complete, group = group, color = group)) +
  xlab('Dropout Threshold') + ylab('Average Percent Complete')

plot_tb_complete = param_tb %>% filter(dropout == 10) %>%
  group_by(distance_method, reconstruction_method, filtering) %>%
  summarise(avg_percent_complete = mean(percent_complete))
plot_tb_complete$group = paste0(plot_tb_complete$reconstruction_method,
                                plot_tb_complete$distance_method)
complete_plot = ggplot(data = plot_tb_complete) +
  geom_line(aes(x = filtering, y = avg_percent_complete, group = group, color = group)) +
  geom_point(aes(x = filtering, y = avg_percent_complete, group = group, color = group)) +
  xlab('Log2 Filtering Threshold') + ylab('Average Percent Complete') +
  ggtitle('Dropout Threshold = 1/100')




#correct
plot_tb_correct = param_tb %>% nest(data = -c(distance_method,
                                               reconstruction_method,
                                               dropout)) %>%
  mutate(opt_filter = map_dbl(data, function(tb) {
    tb$filtering[which(tb$percent_correct == max(tb$percent_correct))][1]
  })) %>% unnest(cols = data) %>% filter(filtering == opt_filter) %>%
  group_by(distance_method,
           reconstruction_method,
           dropout,
           opt_filter) %>%
  summarise(avg_percent_correct = mean(percent_correct))


plot_tb_correct$group = paste0(plot_tb_correct$reconstruction_method,
                                plot_tb_correct$distance_method)
correct_plot = ggplot(data = plot_tb_correct) +
  geom_line(aes(x = dropout/1000, y = opt_filter, group = group, color = group)) +
  geom_point(aes(x = dropout/1000, y = opt_filter, group = group, color = group)) +
  xlab('Dropout Threshold') + ylab('Log2 Optimum Filtering Threshold')

plot_tb_correct = param_tb %>% filter(dropout == 10) %>%
  group_by(distance_method, reconstruction_method, filtering) %>%
  summarise(avg_percent_correct = mean(percent_correct))
plot_tb_correct$group = paste0(plot_tb_correct$reconstruction_method,
                                plot_tb_correct$distance_method)
correct_plot = ggplot(data = plot_tb_correct) +
  geom_line(aes(x = filtering, y = avg_percent_correct, group = group, color = group)) +
  geom_point(aes(x = filtering, y = avg_percent_correct, group = group, color = group)) +
  xlab('Log2 Filtering Threshold') + ylab('Average Percent Correct') +
  ggtitle('Dropout Threshold = 1/100')




plot_tb_complete = param_tb %>% group_by(distance_method,
                                         reconstruction_method,
                                         dropout,
                                         filtering) %>%
  summarise(avg_percent_complete = mean(percent_complete), se = sd(percent_complete)/sqrt(n_sample))
plot_tb_complete$group = paste0(plot_tb_complete$reconstruction_method,
                                plot_tb_complete$distance_method,
                                plot_tb_complete$filtering)

complete_plot = ggplot(data = plot_tb_complete) +
 # geom_ribbon(aes(x = filtering, ymin = avg_percent_complete-se,
 #                  ymax = avg_percent_complete+se, group = group)) +
  geom_line(aes(x = dropout, y = avg_percent_complete, group = group, color = group)) +
  geom_point(aes(x = dropout, y = avg_percent_complete, group = group, color = group)) +
  xlab('Dropout Threshold') + ylab('Percent Complete')

plot_tb_correct = param_tb %>% group_by(reconstruction_method, filtering, n_sample) %>%
  summarise(avg_percent_correct = mean(percent_correct), se = sd(percent_correct)/sqrt(n_sample))
plot_tb_correct$group = paste0(plot_tb_correct$reconstruction_method, plot_tb_correct$filtering)

correct_plot = ggplot(data = plot_tb_correct) +
  geom_errorbar(aes(x = n_sample, ymin = avg_percent_correct-se,
                    ymax = avg_percent_correct+se, group = group, color = group)) +
  geom_line(aes(x = n_sample, y = avg_percent_correct, group = group, color = group))

complete_plot + correct_plot
