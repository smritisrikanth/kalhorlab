library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
library(randomcoloR)
setwd('~/Documents/qfm2')
devtools::load_all()

#functions
dist_mat <- function(data,dist=c('jaccard','manhattan')){
  if(dist == 'jaccard'){ # calculate the jaccard index distance 
    ji_df = data.frame(expand.grid(file1=unique(data$sample),file2=unique(data$sample),d=NA))
    if (nrow(ji_df) < 3) {
      return(NULL)
    }
    for(i in 1:nrow(ji_df)){
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



setwd("~/Documents/R/MosaicFractionAnalysis/R")
devtools::load_all()

setwd('~/Documents/R/Jasmine/')
indelphi_table <- readr::read_tsv('MARC1_indelphi_result.txt', col_names = T)
mut_rate_data <- read_csv('SupplementaryTable3 - TableS3-hgRNA Activities.csv')

setwd('~/Documents/R/David/')
raw_data = load_samples('~/Documents/R/David/5-pair_filtering/')

parent = load_parent('./PB21-founder_filteredpairs.txt')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)

raw_data = filter_unmutated(raw_data)

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
    return(0.5)
  }
})

mice = unique(raw_data$data$mouse)

nest_var_by = c('ID', 'sample', 'mouse', 'tissue', 'age')
nest_density_by = c('mouse', 'tissue', 'age')

nested_data = nest_raw(raw_data, c('ID', 'sample', 'mouse', 'tissue', 'age'))
breaks <- c(-Inf, seq(-10, 0, by = 0.5))
bin_sizes = rev(diff(breaks))
bin_sizes[is.infinite(bin_sizes)] = min(bin_sizes)
nested_data = append_density(nested_data, breaks)
colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
raw_data_filtered = filter_recurring_spacers(raw_data, indelphi_table, 0.0001, filter_missing = FALSE)
nested_data_filtered = nest_raw(raw_data_filtered, c('ID', 'sample', 'mouse', 'tissue', 'age'))
nested_data_filtered = append_density(nested_data_filtered, breaks)
nested_data = adjust_density(nested_data_filtered, nested_data, c('ID', 'sample', 'mouse', 'tissue', 'age'))

nested_data = normalize_density_list(nested_data, apply_cum_adjustment = T)
nested_density = average_density_by(nested_data, nest_density_by, raw = T)

plot_tb_unnormalized = nested_density$nested_density %>% select(mouse, age, tissue, avg_density) %>%
  filter(mouse == mice[[15]]) %>%
  mutate(table = map(avg_density, function(density) {
    density$table
  })) %>%
  unnest(cols = table)
plot_tb_unnormalized$group = paste0(plot_tb_unnormalized$mouse, plot_tb_unnormalized$tissue)

plot_unnormalized = ggplot(data = plot_tb_unnormalized) +
  geom_line(aes(x = -bin, y = count, group = group, color = tissue)) +
  facet_wrap(~ mouse + age) +
  xlab('Negative Log2 Mosaic Fraction') +
  ylab('Adjusted Frequency') +
  theme_pubr() +
  theme(legend.position = "right") + labs(colour = 'Tissue')

plot_unnormalized


#with topology reconstruction
setwd('~/Documents/R/David/')
raw_data = load_samples('~/Documents/R/David/5-pair_filtering/')

parent = load_parent('./PB21-founder_filteredpairs.txt')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)

raw_data = filter_unmutated(raw_data)

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
    return(0.5)
  }
})

mice = unique(raw_data$data$mouse)
mice = mice[mice != 6553]

raw_data$data$tissue_group = substr(raw_data$data$tissue, 1, 1)
dist_mat_list = map(mice, function(m) {
  message(m)
  mat = dist_mat(data = raw_data$data %>% filter(mouse == m) %>%
                   mutate(sample = tissue_group), dist = 'jaccard')
  mat
})
names_6 = rownames(as.matrix(dist_mat_list[[14]]))
names_7 = rownames(as.matrix(dist_mat_list[[2]]))
names_5 = rownames(as.matrix(dist_mat_list[[10]]))
dist_mat_list_to_avg = map(dist_mat_list, function(mat) {
  orig = mat
  mat = as.matrix(mat)
  if (length(intersect(names(orig), names_7)) == 7) {
    ordered_mat = mat[names_7, names_7]
    return(ordered_mat)
  } else {
    return(NULL)
  }
})
null_dist_mats = map_lgl(dist_mat_list_to_avg, function(mat) is.null(mat))
.dist_mat_list = dist_mat_list_to_avg[!null_dist_mats]
cum_mat = Reduce('+', .dist_mat_list)


avg_mat = cum_mat/length(.dist_mat_list)
avg_phy = nj(as.dist(avg_mat))
avg_phy = name_nodes(avg_phy)
avg_phy = phytools::midpoint.root(avg_phy)

raw_data$data$tissue = substr(raw_data$data$tissue, 1, 1)
raw_data = annotate_mut_node(raw_data, avg_phy)

nested_data = nest_raw(raw_data, c('ID', 'sample', 'mouse', 'tissue', 'age', 'gr_node'))
breaks <- c(-Inf, seq(-10, 0, by = 0.5))
bin_sizes = rev(diff(breaks))
bin_sizes[is.infinite(bin_sizes)] = min(bin_sizes)
nested_data = append_density(nested_data, breaks)
colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
raw_data_filtered = filter_recurring_spacers(raw_data, indelphi_table, 0.0001, filter_missing = FALSE)
nested_data_filtered = nest_raw(raw_data_filtered, c('ID', 'sample', 'mouse', 'tissue', 'age', 'gr_node'))
nested_data_filtered = append_density(nested_data_filtered, breaks)
nested_data = adjust_density(nested_data_filtered, nested_data, c('ID', 'sample', 'mouse', 'tissue', 'age', 'gr_node'))

nested_data = normalize_density_list(nested_data, apply_cum_adjustment = T)
nested_density = average_density_by(nested_data, c('ID', 'mouse', 'gr_node'), normalized = T)

mice = unique(nested_density$nested_density$mouse)
IDs = unique(nested_density$nested_density$ID)

pdf('~/Documents/R/David/MF_plots.pdf')
for (m in mice) {
  message(m)
  for (i in IDs) {
    plot_tb = nested_density$nested_density %>% select(mouse, ID, gr_node, avg_density) %>%
      filter(mouse == m) %>% filter(ID == i) %>%
      mutate(table = map(avg_density, function(density) {
        density$table
      })) %>%
      unnest(cols = table)
    plot_tb$group = paste0(plot_tb$mouse, plot_tb$ID, plot_tb$gr_node)
    if(nrow(plot_tb) == 0) {
      next
    }
    plot = ggplot(data = plot_tb) +
      geom_line(aes(x = -bin, y = count, group = group, color = gr_node)) +
      facet_wrap(~ mouse + ID)
    print(plot)
  }
}
dev.off()




#tree reconstruction
#do percent correct/complete on reconstructed trees per mouse in comparison to average tree

pdf(file = './../2_16_mouse_tree_pdf.pdf', onefile = T)
for (m in mice) {
  message(m)
  mat = dist_mat(data = raw_data$data %>% filter(mouse == m) %>%
                   mutate(sample = tissue), dist = 'jaccard')
  if (is.null(mat)) {
    next
  }
  phy = nj(as.dist(mat))
  phy = name_nodes(phy)
  phy = phytools::midpoint.root(phy)
  plot(phy)
  title(m)
}
dev.off()

#dist_mat_list
raw_data$data$tissue_group = substr(raw_data$data$tissue, 1, 1)
dist_mat_list = map(mice, function(m) {
  message(m)
  mat = dist_mat(data = raw_data$data %>% filter(mouse == m) %>%
                   mutate(sample = tissue_group), dist = 'jaccard')
  mat
})
names_6 = rownames(as.matrix(dist_mat_list[[14]]))
names_7 = rownames(as.matrix(dist_mat_list[[2]]))
names_5 = rownames(as.matrix(dist_mat_list[[10]]))

num_tissue = map_dbl(dist_mat_list, function(mat) {
  length(names(mat))
})
names_list = map(dist_mat_list, function(mat) {
  names(mat)
})

dist_mat_df = tibble(mouse = mice, dist_mat = dist_mat_list,
                     num_tissue_group = num_tissue)

dist_mat_list_to_avg = map(dist_mat_list, function(mat) {
  orig = mat
  mat = as.matrix(mat)
  if (length(intersect(names(orig), names_6)) == 6) {
    ordered_mat = mat[names_6, names_6]
    return(ordered_mat)
  } else {
    return(NULL)
  }
})
null_dist_mats = map_lgl(dist_mat_list_to_avg, function(mat) is.null(mat))
.dist_mat_list = dist_mat_list_to_avg[!null_dist_mats]
cum_mat = Reduce('+', .dist_mat_list)

avg_mat = cum_mat/length(dist_mat_list)
avg_phy = nj(as.dist(avg_mat))
avg_phy = name_nodes(avg_phy)
avg_phy = phytools::midpoint.root(avg_phy)

#average tree for combined tissue groups
pdf(file = './../mouse_trees_tissue_groups.pdf', onefile = T)
raw_data$data$tissue_group = substr(raw_data$data$tissue, 1, 1)
for (m in mice) {
  message(m)
  mat = dist_mat(data = raw_data$data %>% filter(mouse == m) %>%
                   mutate(sample = tissue_group), dist = 'jaccard')
  if (is.null(mat)) {
    next
  }
  phy = nj(as.dist(mat))
  phy = name_nodes(phy)
  phy = phytools::midpoint.root(phy)
  plot(phy)
  title(m)
}
dev.off()

dist_mat_df$trees = map(dist_mat_df$mouse, function(m) {
  message(m)
  mat = dist_mat(data = raw_data$data %>% filter(mouse == m) %>%
                   mutate(sample = tissue_group), dist = 'jaccard')
  if (is.null(mat)) {
    next
  }
  phy = nj(as.dist(mat))
  phy = name_nodes(phy)
  phy = phytools::midpoint.root(phy)
  phy
})

dist_mat_df$splits = map(dist_mat_df$trees, function(t) {
  p1 = get_partitions(t)
  # pooled_partitions1 = map(p1, get_pooled_partition)
  pooled_partitions1 = map(p1, function(x) do.call(c, x))
  splits1 = map(pooled_partitions1, get_split)
  unique_splits1 = get_unique_splits(splits1)
  unique_splits1
})

p2 = get_partitions(avg_phy)
pooled_partitions2 = map(p2, get_pooled_partition)
splits2 = map(pooled_partitions2, get_split)
unique_splits2 = get_unique_splits(splits2)

param_tb = dist_mat_df

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

plot_tb = param_tb %>% select(mouse, percent_complete, percent_correct)
plot_tb = plot_tb[!null_dist_mats,]

complete_plot = ggplot(data = plot_tb,
                       aes(x = mouse, y = percent_complete)) +
  geom_col()


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

get_split <- function(pooled_partition, true_tree) {
  complement = setdiff(avg_phy$tip.label, pooled_partition)
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

#percent complete/correct for each mouse
#one plot per mouse per ID

#mutated fraction per mouse per tissue

setwd('~/Documents/R/David/')
raw_data = load_samples('~/Documents/R/David/5-pair_filtering/')

parent = load_parent('./PB21-founder_filteredpairs.txt')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)

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
    return(0.5)
  }
})

nested_data = nest_raw(raw_data, c('mouse', 'tissue', 'ID'))
nested_data = compute_mutation_information(nested_data)

mut_frac_table = nested_data$nested %>% select(mouse, tissue, mutated_fraction, ID) %>%
  filter(!is.na(mouse))

mut_frac_nested = nest(mut_frac_table,
                       data = -c(ID))
logit <- function(p) {
  log((p + 10^-10)/(1 - p + 10^-10))
}
mut_frac_nested$model = map(mut_frac_nested$data, function(tb) {
  lm(formula = logit(mutated_fraction) ~ mouse + tissue, data = tb)
})

pdf('./../mut_frac_mouse_tissue.pdf')
for (i in unique(mut_frac_table$ID)) {
  message(i)
  plot = ggplot(data = mut_frac_table %>% filter(ID == i),
         aes(x = mouse, y = mutated_fraction, color = tissue)) +
    geom_point() + ggtitle(id)
  print(plot)
}
dev.off()

ggplot(data = mut_frac_table,
       aes(x = mouse, y = mutated_fraction, color = tissue)) +
  geom_point()
plot
dev.off()
