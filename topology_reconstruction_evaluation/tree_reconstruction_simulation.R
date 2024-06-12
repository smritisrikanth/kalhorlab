#functions

dist_mat <- function(data,dist=c('jaccard','manhattan')){
  if(dist == 'jaccard'){ # calculate the jaccard index distance 
    ji_df = data.frame(expand.grid(file1=unique(data$sample),file2=unique(data$sample),d=NA))
    if (nrow(ji_df) < 3) {
      return(NULL)
    }
    for(i in 1:nrow(ji_df)){
      if (i %% 1000 == 0) message(i)
      if(ji_df$file1[i] != ji_df$file2[i]){
        file1 = data %>% filter(sample == ji_df$file1[i]) %>% 
          transmute(ID = paste0(ID,'_',sequence)) %>% pull(ID) %>% unique(.)
        file2 = data %>% filter(sample == ji_df$file2[i]) %>% 
          transmute(ID = paste0(ID,'_',sequence)) %>% pull(ID) %>% unique(.)
        s = length(intersect(file1,file2)) / length(union(file1,file2))
        ji_df$d[i]=1-s
      }else{
        ji_df$d[i]=0
      }
    }
    D = ji_df %>% pivot_wider(names_from = file1, values_from = d) %>%
      column_to_rownames(var='file2')
    return(D)
  } else{ # add manhattan distance here 
    log_binary_matrix = normalized_matrix(data, value_choice = 'logMF')
    return(dist(log_binary_matrix,method = "manhattan"))
    
  }
  phy = name_nodes(phy)
  phy = phytools::midpoint.root(phy)
  return(phy)
}


build_tree <- function(data,dist=c('jaccard','manhattan')){
  if(dist == 'jaccard'){ # calculate the jaccard index distance 
    ji_df = data.frame(expand.grid(file1=unique(data$sample),file2=unique(data$sample),d=NA))
    if (nrow(ji_df) < 3) {
      return(NULL)
    }
    for(i in 1:nrow(ji_df)){
      if (i %% 1000 == 0) message(i)
      if(ji_df$file1[i] != ji_df$file2[i]){
        file1 = data %>% filter(sample == ji_df$file1[i]) %>% 
          transmute(ID = paste0(ID,'_',sequence)) %>% pull(ID) %>% unique(.)
        file2 = data %>% filter(sample == ji_df$file2[i]) %>% 
          transmute(ID = paste0(ID,'_',sequence)) %>% pull(ID) %>% unique(.)
        s = length(intersect(file1,file2)) / length(union(file1,file2))
        ji_df$d[i]=1-s
      }else{
        ji_df$d[i]=0
      }
    }
    D = ji_df %>% pivot_wider(names_from = file1, values_from = d) %>%
      column_to_rownames(var='file2')
    phy=nj(as.dist(D))
  } else{ # add manhattan distance here 
    log_binary_matrix = normalized_matrix(data)
    phy= nj(dist(log_binary_matrix,method = "manhattan"))
    
  }
  phy = name_nodes(phy)
  phy = phytools::midpoint.root(phy)
  return(phy)
}

plot_tree <- function(phy, tip_label, node_col, type=c('igraph','rooted')){
  if(type == 'igraph'){
    phy_ig = as.igraph(phy)
    names <- c()
    colors <- c()
    for(i in 1:length(V(phy_ig)$name)){
      if(grepl('txt',V(phy_ig)$name[i])){
        name = map_chr(str_split(V(phy_ig)$name[i],'_'),1)
        color = map_chr(str_split(name,'-'),1)
        # color = map_chr(str_split(name,'-'),3)
        # color = map_chr(str_split(name,'-'), function(x_str) paste0(x_str[1:2], collapse = "_"))
      }else{
        name = V(phy_ig)$name[i]
        color = 'node'
      }
      names <- c(names,name)
      colors <- c(colors,color)
    }
    
    V(phy_ig)$name <- names
    V(phy_ig)$random_col = colors
    p <- ggraph(phy_ig, layout = "igraph", algorithm = 'kk') +
      # geom_node_label(aes(label = name, fill = random_col)) +
      # scale_fill_manual(values = node_col) +
      geom_edge_diagonal() +
      geom_node_point(aes(label = name, color = random_col))
  }else{
    p <- plot_gr_clean(phy,
                       type_col = node_col,
                       # define color 
                       # type_col = c("Node-1" = "black",
                       #              "Node-2" = "red",
                       #              "Root" = "yellow",
                       #              "yo" = "blue",
                       #              'ta' = "red",
                       #              'he'='orange',
                       #              'pl'='purple'
                       # ),
                       show_node_label = T)
  }
  return(p)
}


normalized_matrix <- function(data, value_choice = c("logMF", "binary")){
  data_2 = data %>% transmute(ID = paste0(ID,'_',sequence),reads=mosaic_fraction,sample=sample) %>% 
    pivot_wider(names_from = ID,values_from = reads) %>%
    column_to_rownames(.,var='sample') %>%
    mutate_all(~ifelse(is.na(.), 0, .))
  normalized_matrix <- t(apply(data_2, 1, function(x) x/sum(x)))
  if(value_choice == 'logMF'){
    log_matrix <- log2(normalized_matrix + 1)
  }else{
    log_matrix <- 1*(log(normalized_matrix + 1) > 0)
  }
  return(log_matrix)
}


load_simulated_data <- function(path, num_samples, mat_name = 'mut_frac_mat') {
  param_tb = tibble(sim = num_samples)
  param_tb$filename = paste0('mf_tables_', param_tb$sim, '.rda')
  param_tb[[mat_name]] = map(param_tb$filename, function(fn) {
    if (!file.exists(paste0(path, fn))) {
      return(NULL)
    }
    load(paste0(path, fn))
    mut_frac_mat %>% filter(mosaic_fraction > 0) %>% nest(data = -c(ID, sequence)) %>%
      mutate(avg_log2mf = map_dbl(data, function(data) {log2(mean(data$mosaic_fraction))})) %>% unnest(cols = data)
  })
  param_tb = param_tb[!is.null(param_tb[[mat_name]]),]
  param_tb_unnest = unnest(param_tb, cols = mat_name) %>%
    select(sim, ID, tissue = CellType, sequence, mosaic_fraction, avg_log2mf, node_time, probability)
  # param_tb_unnest = unnest(param_tb, cols = mat_name) %>%
  #   select(sim, ID, tissue = CellType, sequence, mosaic_fraction, avg_log2mf, node_time)
  param_tb_unnest$sample = paste0(param_tb_unnest$sim, param_tb_unnest$tissue)
  param_tb_unnest = unique(param_tb_unnest)
  raw_data = structure(list(data = param_tb_unnest), class = 'raw')
  raw_data
}


library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
#library(randomcoloR)
setwd('/home/ssrikan2/data-kreza1/smriti/qfm2')
devtools::load_all()

setwd('/home/ssrikan2/data-kreza1/smriti/MF_Signal_Simulation')
job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

input_folder = 'mouse_gas_cg_ss1000_flat_hgRNA'
output_folder = 'mouse_ss1000_hgRNA_flat_dist_mats_v3'



phy = readRDS("./gast_phylo.rds")
raw_data = load_simulated_data('./mouse_gas_cg_ss1000_flat_hgRNA/', job_id)

load('param_tb.rda')
# param_tb = as_tibble(expand.grid(distance_method = c('manhattan',
#                                                      'jaccard'),
#                                  filtering = c(-10:0),
#                                  dropout = c(1:10)))

param_tb$mat = map(1:nrow(param_tb), function(n) {
  print(n)
  distance_method = param_tb$distance_method[n]
  filtering = param_tb$filtering[n]
  dropout = param_tb$dropout[n]
  p = 2^(param_tb$filtering[n])
  d = param_tb$dropout[n]/1000
  
  raw_data_filtered = raw_data
  raw_data_filtered$data = raw_data$data[raw_data$data$probability < p,]
  raw_data_filtered$data = raw_data_filtered$data[raw_data_filtered$data$mosaic_fraction > d,]
  data = raw_data_filtered$data
  
  dist_mat = dist_mat(data = data %>% mutate(sample = tissue), dist = distance_method)
  save(dist_mat, file = paste0('./', output_folder, '/dist_mat_', distance_method, '_', filtering, '_', dropout, '_', job_id, '.rda'))
  dist_mat
})



