library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
library(randomcoloR)
setwd('~/Documents/qfm2')
devtools::load_all()

setwd("~/Documents/R/MosaicFractionAnalysis/R")
devtools::load_all()

setwd('~/Documents/R/Jasmine/')
indelphi_table <- readr::read_tsv('MARC1_indelphi_result.txt', col_names = T)
mut_rate_data <- read_csv('SupplementaryTable3 - TableS3-hgRNA Activities.csv')

setwd('~/Documents/R/David/')
raw_data = load_samples('~/Documents/R/David/')

parent = load_parent('./PB21-founder_filteredpairs.txt')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)

#raw_data = filter_unmutated(raw_data)

sample_metadata = as_tibble(stringr::str_match(unique(raw_data$data$sample), "(\\d+)-(.*?)_filteredpairs\\.txt"))
colnames(sample_metadata) = c('sample', 'mouse', 'tissue')
raw_data = append_metadata(raw_data, sample_metadata)
colnames(mut_rate_data)[3] = 'ID'
raw_data = append_metadata(raw_data, mut_rate_data, cols = c('Class'), join_by = 'ID')
#raw_data$data = unique(raw_data$data)
animal_age_list = c('5728' = 1.7, '5729' = 1.7, '6554' = 0.58, '6555' = 0.58, '6351' = 1.36,
                    '6352' = 1.36, '6559' = 1.1, '65576558' = 1.1, '5662' = 2.4, '65606561' = 1.1)
raw_data$data$age = map_dbl(raw_data$data$mouse, function(mouse) {
  if (mouse %in% names(animal_age_list)) {
    animal_age_list[[mouse]]
  } else {
    return(0)
  }
})

nested_data_mutated_fraction = nest_raw(raw_data, c('ID', 'sample', 'mouse', 'tissue', 'age', 'Class'))
nested_data_mutated_fraction = compute_mutation_information(nested_data_mutated_fraction)
ID_list = unique(nested_data_mutated_fraction$nested$ID)

data = select(nested_data_mutated_fraction$nested, ID, mouse, tissue, age, mutated_fraction, Class)
  
ggplot(data = data %>% filter(mouse == '65576558'), aes(x = tissue, y = mutated_fraction)) + geom_boxplot() +
  facet_wrap(~ Class)

raw_data = filter_unmutated(raw_data)
breaks <- c(-Inf, seq(-10, 0, by = 0.5))
bin_sizes = rev(diff(breaks))
bin_sizes[is.infinite(bin_sizes)] = min(bin_sizes)
nested_data = nest_raw(raw_data, c('ID', 'sample', 'mouse', 'tissue', 'age'))
nested_data = append_density(nested_data, breaks)

colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
raw_data_filtered = filter_recurring_spacers(raw_data, indelphi_table, 0.001, filter_missing = FALSE)
nested_data_filtered = nest_raw(raw_data_filtered, c('ID', 'sample', 'mouse', 'tissue', 'age'))
nested_data_filtered = append_density(nested_data_filtered, breaks)
nested_data = adjust_density(nested_data_filtered, nested_data, c('ID', 'sample', 'mouse', 'tissue', 'age'))


nested_data = normalize_density_list(nested_data)
nested_data$nested$tissue_group = substr(nested_data$nested$tissue, 1, 1)
nested_density = average_density_by(nested_data, c('ID', 'mouse', 'tissue_group', 'tissue', 'age'), normalized = TRUE)

plot_tb = nested_density$nested_density %>% select(mouse, tissue_group, ID, tissue, age, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    loess = loess(formula = log(count+0.0001) ~ bin,
                  data = density$table)
    table = tibble(log2mf = seq(-10, 0, by = 0.01),
                   smooth_density = exp(predict(loess, newdata = tibble(bin = seq(-10, 0, by = 0.01))))-0.0001)
    table
  })) %>% unnest(cols = table) %>%
  filter(ID == ID_list[[25]]) %>% filter(mouse == '6351')
plot_tb$group = paste0(plot_tb$mouse,
                       plot_tb$tissue,
                       plot_tb$age)
plot_tb$age = factor(plot_tb$age)


ggplot(data = plot_tb, aes(x = -log2mf, y = smooth_density, group = group, color = tissue_group)) + geom_line()


#with node assignment
phy = build_tree(data = raw_data$data %>% mutate(sample = tissue)
                 %>% filter(mouse == '65576558'), dist = 'jaccard')
node_colors = map(phy$node.label, function(node) {
  randomColor()
})
names(node_colors) = phy$node.label
plot_tree(phy = phy, type = 'igraph', tip_label = phy$tip.label, node_col = node_colors)
plot_gr_clean(phy, show_node_label = T)

pdf(file = './david_data_each_mouse_tree.pdf', onefile = T)
for (m in unique(raw_data$data$mouse)) {
  data = raw_data$data %>% mutate(sample = tissue) %>% filter(mouse == m)
  phy = build_tree(data = data, dist = 'jaccard')
  if (!is.null(phy)) {
    print(plot(phy, show.node.label = T))
  } else {
    print(m)
  }
}
dev.off()
data = raw_data$data %>% mutate(sample = tissue) %>% filter(mouse == '65576558')
raw_data_test = raw_data
raw_data_test$data = data
raw_data_test = annotate_mut_node(raw_data_test, phy)
nested_data_test = nest_raw(raw_data_test, c('ID', 'sample', 'gr_node'))
nested_data_test = append_density(nested_data_test, breaks)

colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
raw_data_filtered_test = filter_recurring_spacers(raw_data_test, indelphi_table, 1, filter_missing = FALSE)
nested_data_filtered_test = nest_raw(raw_data_filtered_test, c('ID', 'sample', 'gr_node'))
nested_data_filtered_test = append_density(nested_data_filtered_test, breaks)
nested_data_test = adjust_density(nested_data_filtered_test, nested_data_test, c('ID', 'sample', 'gr_node'))

node_order_list = c(1:39)
names(node_order_list) = c(phy$node.label, phy$tip.label)
nested_data_test = increment_cum_mf(nested_data_test, node_order_list)
nested_data_test = normalize_density_list(nested_data_test)
nested_data_test$nested$condition = substr(nested_data_test$nested$gr_node, 1, 1)
#nested_data_test$nested$condition = nested_data_test$nested$gr_node
nested_density_test = average_density_by(nested_data_test, c('condition'), normalized = T)
log2mf_df = align_counts(nested_density_test)
plot_tb = gather(data = log2mf_df,
                           key = 'condition', value = 'count', -c(log2mf))
ggplot(data = plot_tb, aes(x = -log2mf, y = count, group = condition, color = condition)) +
  geom_line() + xlab('Negative Log2MF') + ylab('Frequency') +
  ggtitle('Unfiltered') + labs(colour = 'Node')

violin_plot_data = nested_data_test$nested
violin_plot_data$raw_mf = map(violin_plot_data$data, function(data) {
  data$log2MF
})
violin_plot_tb = violin_plot_data %>% ungroup() %>% select(condition, raw_mf) %>%
  unnest(cols = raw_mf)

ggviolin(data = violin_plot_tb, x = 'condition', y = 'raw_mf')

fsenames = !!(names(log2mf_df))

fig <- plot_ly(log2mf_df, x = ~(-log2mf), y = , name = names[2], type = 'scatter', mode = 'lines') 
  
#functions
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
    phy=nj(dist(log_binary_matrix,method = "manhattan"))
    
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

