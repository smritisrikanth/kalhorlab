library(ggpubr)
library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('~/Documents/qfm2')
devtools::load_all()

setwd("~/Documents/R/MosaicFractionAnalysis/R")
devtools::load_all()

setwd('~/Documents/R/Jasmine/')
indelphi_table <- readr::read_tsv('MARC1_indelphi_result.txt', col_names = T)

raw_data = load_samples_true_pairs('~/Documents/R/Jessica/240111_Bulk_Brain_Files/')
sample_metadata = as_tibble(stringr::str_match(unique(raw_data$data$sample), "(.*?)_(.*?)_(\\d+)_(\\d+_\\d+)-(.*?)_truepairs.txt"))[1:5]
colnames(sample_metadata) = c('sample', 'side', 'region', 'rep', 'tech_rep')
raw_data = append_metadata(raw_data, sample_metadata)
raw_data$data = raw_data$data %>% group_by(ID, sequence, side, region, rep) %>%
  summarise(reads = sum(reads))
raw_data$data = raw_data$data %>% mutate(sample = paste0(region, side))
parent = load_parent('./../../David/PB21-founder_filteredpairs.txt')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)
raw_data = filter_unmutated(raw_data)
raw_data$data = unique(raw_data$data)
raw_data$data = raw_data$data[!is.na(raw_data$data$sample),]

nested_data = nest_raw(raw_data, c('ID', 'sample', 'side', 'region', 'rep'))
breaks <- c(-Inf, seq(-10, 0, by = 0.5))
bin_sizes = rev(diff(breaks))
bin_sizes[is.infinite(bin_sizes)] = min(bin_sizes)
nested_data = append_density(nested_data, breaks)
colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
raw_data = filter_recurring_spacers(raw_data, indelphi_table, 0.001, filter_missing = FALSE)
nested_data_filtered = nest_raw(raw_data, c('ID', 'sample', 'side', 'region', 'rep'))
nested_data_filtered = append_density(nested_data_filtered, breaks)
nested_data = adjust_density(nested_data_filtered, nested_data, c('ID', 'sample', 'side', 'region', 'rep'))

nested_data = normalize_density_list(nested_data)
nested_density = average_density_by(nested_data, c('side', 'region'), normalized = TRUE)

plot_tb_unnormalized = nested_density$nested_density %>% select(side, region, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    density$table
  })) %>%
  filter(!is.na(side)) %>%
  unnest(cols = table)

plot_unnormalized = ggplot(data = plot_tb_unnormalized) +
  geom_line(aes(x = -bin, y = count, group = region, color = region)) +
  facet_wrap(~ side) +
  xlab('Negative Log2 Mosaic Fraction') +
  ylab('Adjusted Frequency') +
  theme_pubr() +
  theme(legend.position = "right") + labs(colour = 'Side')



#same pipeline with mutations assigned to nodes in reconstructed topology
raw_data = load_samples_true_pairs('~/Documents/R/Jessica/240111_Bulk_Brain_Files/')
sample_metadata = as_tibble(stringr::str_match(unique(raw_data$data$sample), "(.*?)_(.*?)_(\\d+)_(\\d+_\\d+)-(.*?)_truepairs.txt"))[1:5]
colnames(sample_metadata) = c('sample', 'side', 'region', 'rep', 'tech_rep')
raw_data = append_metadata(raw_data, sample_metadata)
raw_data$data = raw_data$data %>% group_by(ID, sequence, side, region) %>%
  summarise(reads = sum(reads))
raw_data$data = raw_data$data[!is.na(raw_data$data$region),]
raw_data$data = raw_data$data %>% mutate(sample = paste0(region, '_', side))
parent = load_parent('./../../David/PB21-founder_filteredpairs.txt')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)
raw_data = filter_unmutated(raw_data)
raw_data$data = unique(raw_data$data)
raw_data$data = raw_data$data[!is.na(raw_data$data$sample),]

#reconstruct topology
mat = dist_mat(data = raw_data$data, dist = 'jaccard')
phy = nj(as.dist(mat))
#phy = ape::as.phylo(hclust(as.dist(mat), method = 'average'))
#phy = ape::as.phylo(hclust(as.dist(mat), method = 'ward.D2'))
phy = name_nodes(phy)
phy = phytools::midpoint.root(phy)

# mat = dist_mat(data = raw_data$data, dist = 'manhattan')
# phy = nj(as.dist(mat))
# phy = name_nodes(phy)
# phy = phytools::midpoint.root(phy)

#assign mutations to nodes
raw_data$data = raw_data$data %>% mutate(tissue = sample)
raw_data = annotate_mut_node(raw_data, phy)

nested_data = nest_raw(raw_data, c('ID', 'sample', 'gr_node'))
breaks <- c(-Inf, seq(-15, 0, by = 0.5))
bin_sizes = rev(diff(breaks))
bin_sizes[is.infinite(bin_sizes)] = min(bin_sizes)
nested_data = append_density(nested_data, breaks)
colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
raw_data_filtered = filter_recurring_spacers(raw_data, indelphi_table, 0.0001, filter_missing = FALSE)
nested_data_filtered = nest_raw(raw_data_filtered, c('ID', 'sample', 'gr_node'))
nested_data_filtered = append_density(nested_data_filtered, breaks)
nested_data = adjust_density(nested_data_filtered, nested_data, c('ID', 'sample', 'gr_node'))

nested_data = normalize_density_list(nested_data, apply_cum_adjustment = T)
nested_density = average_density_by(nested_data, c('gr_node'))

plot_tb = nested_density$nested_density %>% select(gr_node, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    density$table
  })) %>% filter(substr(gr_node, 1, 4) != 'Node') %>%
  unnest(cols = table)
plot_tb$region = map_chr(strsplit(plot_tb$gr_node, '_'), function(x) {
  x[1]
})
plot_tb$side = map_chr(strsplit(plot_tb$gr_node, '_'), function(x) {
  x[2]
})

plot_tb_smooth = nested_density$nested_density %>% select(gr_node, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    #density$table$count = density$table$count/sum(density$table$count)/bin_sizes
    loess = loess(formula = log(count+0.0001) ~ bin,
                  data = density$table,
                  span = 0.3)
    table = tibble(log2mf = seq(-15, 0, by = 0.01),
                   smooth_density = exp(predict(loess, newdata = tibble(bin = seq(-15, 0, by = 0.01))))-0.0001)
    table$max_mf = table$log2mf[which.max(table$smooth_density)]
    # density$table
    table
  })) %>% filter(substr(gr_node, 1, 4) != 'Node') %>%
  unnest(cols = table) %>%
  mutate(region = unlist(strsplit(gr_node, '_'))[1]) %>%
  mutate(side = unlist(strsplit(gr_node, '_'))[2])
plot_tb_smooth$region = map_chr(strsplit(plot_tb_smooth$gr_node, '_'), function(x) {
  x[1]
})
plot_tb_smooth$side = map_chr(strsplit(plot_tb_smooth$gr_node, '_'), function(x) {
  x[2]
})

plot = ggplot(data = plot_tb) +
  geom_line(aes(x = -bin, y = count, group = gr_node, color = region)) +
  facet_wrap(~ side) +
  xlab('Negative Log2 Mosaic Fraction') +
  ylab('Adjusted Frequency') +
  theme_pubr() +
  theme(legend.position = "right") + labs(colour = 'Node')

plot_smooth = ggplot(data = plot_tb_smooth) +
  geom_line(aes(x = -log2mf, y = smooth_density, group = gr_node, color = region)) +
  facet_wrap(~ side) +
  xlab('Negative Log2 Mosaic Fraction') +
  ylab('Adjusted Frequency') +
  theme_pubr() +
  theme(legend.position = "right") + labs(colour = 'Node') 

plot / plot_smooth

#functions
dist_mat <- function(data,dist=c('jaccard','manhattan')){
  if(dist == 'jaccard'){ # calculate the jaccard index distance 
    ji_df = data.frame(expand.grid(file1=unique(data$sample),file2=unique(data$sample),d=NA))
    if (nrow(ji_df) < 3) {
      return(NULL)
    }
    for(i in 1:nrow(ji_df)){
      if (i %% 10 == 0) message(i)
      message(i)
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

