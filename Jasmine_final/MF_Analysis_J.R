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
mut_rate_data <- read_csv('SupplementaryTable3 - TableS3-hgRNA Activities.csv')

raw_data = load_samples('~/Downloads/Jasmine Data (1)/MARC1-Pipeline-MARCC-master/5-pair_filtering/')
parent = load_parent('./PB3-founder_filteredpairs.txt')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)
raw_data = filter_unmutated(raw_data)
colnames(mut_rate_data)[3] = 'ID'
raw_data = append_metadata(raw_data, mut_rate_data, cols = c('Class'), join_by = 'ID')
sample_metadata = as_tibble(stringr::str_match(unique(raw_data$data$sample), "(t\\d+)-(\\d+)-(.*?)_filteredpairs\\.txt"))
colnames(sample_metadata) = c('sample', 'mother', 'mouse', 'tissue')
sample_metadata$condition = map_chr(sample_metadata$mother, function(m) {
  if (grepl('t67', m)) {
    'MARC1'
  } else {
    'Cas9'
  }
})
raw_data = append_metadata(raw_data, sample_metadata)
raw_data$data$tissue = substr(raw_data$data$tissue, 1, 2)
raw_data$data = raw_data$data[raw_data$data$tissue != 'MA' &
                                raw_data$data$tissue != 'pl',]
raw_data$data = unique(raw_data$data)
#phy = name_nodes(ape::read.tree(text = "(yo,(he,ta));"))
#raw_data = annotate_mut_node(raw_data, phy)
nested_data = nest_raw(raw_data, c('ID', 'sample', 'mouse', 'tissue', 'Class', 'condition'))
colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
raw_data = filter_recurring_spacers(raw_data, indelphi_table, 0.001, filter_missing = FALSE)
nested_data_filtered = nest_raw(raw_data, c('ID', 'sample', 'mouse', 'tissue', 'Class', 'condition'))
breaks <- c(-Inf, seq(-10, 0, by = 0.5))
bin_sizes = rev(diff(breaks))
bin_sizes[is.infinite(bin_sizes)] = min(bin_sizes)
nested_data = append_density(nested_data, breaks)
nested_data_filtered = append_density(nested_data_filtered, breaks)
nested_data = adjust_density(nested_data_filtered, nested_data, c('ID', 'sample', 'mouse', 'tissue', 'Class', 'condition'))

# node_order_list = c(1:5)
# names(node_order_list) = c(phy$node.label, phy$tip.label)
# nested_data = increment_cum_mf(nested_data, node_order_list)

nested_data = normalize_density_list(nested_data)

nested_density = average_density_by(nested_data, c('tissue', 'Class', 'condition'), normalized = TRUE)
#names(nested_density$nested_density)[which(names(nested_density$nested_density) == 'gr_node')] = 'condition'

# log2mf_df = align_counts(nested_density)
# log2mf_df_norm = log2mf_df %>% select(-c(log2mf)) %>%
#   mutate(across(everything(), ~ ./sum(.)/bin_sizes)) %>%
#   mutate(log2mf = log2mf_df$log2mf)

plot_tb_unnormalized = nested_density$nested_density %>% select(Class, tissue, condition, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    #density$table$count = density$table$count/sum(density$table$count)/bin_sizes
    density$table
  })) %>%
  unnest(cols = table) %>%
  filter(Class != 'inactive')
plot_tb_unnormalized$group = paste0(plot_tb_unnormalized$Class, plot_tb_unnormalized$tissue, plot_tb_unnormalized$condition)

nested_density$nested_density$num_mice = map_dbl(nested_density$nested_density$density_list, function(density) {
  length(unique(density$mouse))
})
plot_tb_smooth_unnormalized = nested_density$nested_density %>% select(num_mice, Class, tissue, condition, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    #density$table$count = density$table$count/sum(density$table$count)/bin_sizes
    loess = loess(formula = log(count+0.0001) ~ bin,
                  data = density$table)
    table = tibble(log2mf = seq(-10, 0, by = 0.01),
                   smooth_density = exp(predict(loess, newdata = tibble(bin = seq(-10, 0, by = 0.01))))-0.0001)
    table$max_mf = table$log2mf[which.max(table$smooth_density)]
    table
  })) %>%
  unnest(cols = table) %>%
  filter(Class != 'inactive')
plot_tb_smooth_unnormalized$group = paste0(plot_tb_smooth_unnormalized$Class,
                                           plot_tb_smooth_unnormalized$tissue,
                                           plot_tb_smooth_unnormalized$condition)

plot_smooth_unnormalized = ggplot(data = plot_tb_smooth_unnormalized) +
  geom_line(aes(x = -log2mf, y = smooth_density, group = group, color = condition)) +
  facet_wrap(~ Class + tissue) +
  xlab('Negative Log2 Mosaic Fraction') +
  ylab('Adjusted Frequency') +
  theme_pubr() +
  theme(legend.position = "right") + labs(colour = 'Mother') +
  ylim(0,0.05)

plot_tb = nested_density$nested_density %>% select(Class, tissue, condition, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    density$table$count = density$table$count/sum(density$table$count)/bin_sizes
    density$table
  })) %>%
  unnest(cols = table) %>%
  filter(Class != 'inactive')
plot_tb$group = paste0(plot_tb$Class, plot_tb$tissue, plot_tb$condition)

plot_tb_smooth = nested_density$nested_density %>% select(Class, tissue, condition, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    #density$table$count = density$table$count/sum(density$table$count)/bin_sizes
    loess = loess(formula = log(count+0.0001) ~ bin,
                  data = density$table)
    table = tibble(log2mf = seq(-10, 0, by = 0.01),
                   smooth_density = exp(predict(loess, newdata = tibble(bin = seq(-10, 0, by = 0.01))))-0.0001)
    table$max_mf = table$log2mf[which.max(table$smooth_density)]
    table
  })) %>%
  unnest(cols = table) %>%
  filter(Class != 'inactive')
plot_tb_smooth$group = paste0(plot_tb_smooth$Class, plot_tb_smooth$tissue, plot_tb_smooth$condition)

plot_unnormalized = ggplot(data = plot_tb_unnormalized) +
  geom_line(aes(x = -bin, y = count, group = group, color = Class)) +
  facet_wrap(~ condition + tissue) +
  xlab('Negative Log2 Mosaic Fraction') +
  ylab('Adjusted Frequency') +
  theme_pubr() +
  theme(legend.position = "right") + labs(colour = 'Class')

num_mice = nested_density$nested_density %>% filter(Class != 'inactive') %>%
  group_by(Class, tissue) %>%
  summarize(num = sum(num_mice))
plot_smooth_unnormalized = ggplot(data = plot_tb_smooth_unnormalized) +
  geom_line(aes(x = -log2mf, y = smooth_density, group = group, color = condition)) +
  facet_wrap(~ tissue + Class) +
  xlab('Negative Log2 Mosaic Fraction') +
  ylab('Adjusted Frequency') +
  theme_pubr() +
  labs(colour = 'Mother Background') +
  geom_text(data = num_mice, size = 2,
            aes(x = 2, y = 0.15, label = paste0('n = ', num))) +
  theme(legend.position = 'right',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 9))

png(filename = '~/Downloads/color_mother.png', width = 4, height = 4, units = "in", res = 500)
print(plot_smooth_unnormalized)
dev.off()
  

plot = ggplot(data = plot_tb) +
  geom_line(aes(x = -bin, y = count, group = group, color = Class)) +
  facet_wrap(~ condition + tissue) +
  xlab('Negative Log2 Mosaic Fraction') +
  ylab('Adjusted Density') +
  theme_pubr() +
  theme(legend.position = "right") + labs(colour = 'Class')

plot_smooth = ggplot(data = plot_tb_smooth) +
  geom_line(aes(x = -log2mf, y = smooth_density, group = group, color = condition)) +
  xlab('Negative Log2 Mosaic Fraction') +
  ylab('Adjusted Density') +
  theme_pubr() +
  theme(legend.position = "right") + labs(colour = 'Class')

(plot_unnormalized + plot_smooth_unnormalized) / (plot + plot_smooth)

#with nodes
plot_tb = gather(data = log2mf_df_norm %>% select(-c(average, average_nonzero)),
                 key = 'condition', value = 'count', -c(log2mf))
plot_tb_avg = gather(data = log2mf_df_norm %>% select(c(average, average_nonzero, log2mf)),
                 key = 'condition', value = 'count', -c(log2mf))

plot_nodes = ggplot(data = plot_tb) +
  geom_line(aes(x = -log2mf, y = count, group = condition, color = condition)) +
  ylab('Mosaic Fraction Density') +
  xlab('Negative Log2 Mosaic Fraction') + theme_pubr() +
  theme(legend.position = 'right') + labs(colour = 'Cell Type') +
  ggtitle('Stratified by Cell Type')

plot_avg = ggplot(data = plot_tb_avg) +
  geom_line(aes(x = -log2mf, y = count, group = condition, color = condition)) +
  scale_colour_manual(values = c('average' = 'black', 'average_nonzero' = 'gray')) +
  ylab('Mosaic Fraction Density') +
  xlab('Negative Log2 Mosaic Fraction') + theme_pubr() +
  theme(legend.position = 'right') + labs(colour = 'Signal Stiching Method') +
  ggtitle('Average Recovered Signal')

plot_nodes / plot_avg + theme(text = element_text(size = 10))
