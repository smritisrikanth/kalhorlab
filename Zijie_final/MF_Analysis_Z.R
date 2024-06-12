source("MF_Functions.R")
source("MF_util.R")

library(tidyverse)
library(ggpubr)
library(tidyr)
library(purrr)
library("readxl")
library(patchwork)

setwd('~/Documents/R/MosaicFractionAnalysis/')
devtools::load_all()

# setwd('~/Documents/R/filteredpairs_Z')
# #parent <- read.csv('PB7-founder_filteredpairs.txt', sep = '\t', col.names = c('ID', 'sequence', 'reads', '%ID', '%total'))
# parent = load_parent('PB7-founder_filteredpairs.txt')
# breaks <- c(-Inf, seq(-10, 0, by = 1))
# by_var = c('ID', 'sample')

# compute average normalized density
raw_data = load_samples('~/Documents/R/Zijie/filteredpairs_Z/')
parent = load_parent('~/Documents/R/Zijie/filteredpairs_Z/PB7-founder_filteredpairs.txt')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)
raw_data = filter_unmutated(raw_data)
raw_data$data$condition = substr(raw_data$data$sample, 10, 10)
raw_data$data = raw_data$data[raw_data$data$condition != 'A',]
raw_data$data$condition[raw_data$data$condition == 'G'] = 'A'
nested_data = nest_raw(raw_data, c('ID', 'sample', 'condition'))
id_sample_condition = select(nested_data$nested, ID, sample, condition)
sparse_samples = id_sample_condition %>% group_by(sample) %>%
  summarise(num_id = length(unique(ID))) %>% filter(num_id < 5)
sparse_ids = id_sample_condition %>% group_by(ID) %>%
  summarise(num_sample = length(unique(sample))) %>% filter(num_sample < 10)
#colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
#raw_data = filter_recurring_spacers(raw_data, indelphi_table, 0.01, filter_missing = FALSE)
recur_mut_one <- read_excel("~/Documents/R/Zijie/recurring mutation/one_percent_recurring_mutation.xlsx")
recur_mut_pointone <- read_excel("~/Documents/R/Zijie/recurring mutation/pointone_percent_recurring_mutation.xlsx")
raw_data$data = raw_data$data[!(raw_data$data$sequence %in% recur_mut_pointone$allele),]
raw_data$data = raw_data$data[raw_data$data$log2MF < 0,]
raw_data$data = raw_data$data[!(raw_data$data$ID %in% sparse_ids$ID),]
raw_data$data = raw_data$data[!(raw_data$data$sample %in% sparse_samples$sample),]
nested_data_filtered = nest_raw(raw_data, c('ID', 'sample', 'condition'))
breaks <- c(-Inf, seq(-10, 0, by = 0.5))
bin_sizes = rev(diff(breaks))
bin_sizes[is.infinite(bin_sizes)] = min(bin_sizes)
nested_data = append_density(nested_data, breaks)
nested_data_filtered = append_density(nested_data_filtered, breaks)
nested_data = adjust_density(nested_data_filtered, nested_data, c('ID', 'sample', 'condition'))
nested_data = normalize_density_list(nested_data)
# nested_data$nested$sum_mf = map_dbl(nested_data$nested$data, function(data) {
#   sum(data$mosaic_fraction)
# })

nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)
nested_density_id_cum_mf = average_density_by(nested_data, c('ID'), cum_mf = TRUE)

plot_tb_smooth = nested_density$nested_density %>% select(condition, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    loess = loess(formula = log(count+0.0001) ~ bin,
                  data = density$table)
    table = tibble(log2mf = seq(-10, 0, by = 0.01),
           smooth_density = exp(predict(loess, newdata = tibble(bin = seq(-10, 0, by = 0.01))))-0.0001)
    table$max_mf = table$log2mf[which.max(table$smooth_density)]
    table
  })) %>%
  select(condition, table) %>%
  unnest(cols = table)

plot_tb = nested_density$nested_density %>% select(condition, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    #density$table$count = density$table$count/sum(density$table$count)/bin_sizes
    tibble(log2mf = density$table$bin, density = density$table$count)
  })) %>%
  select(condition, table) %>%
  unnest(cols = table)


plot = ggplot(plot_tb) +
  geom_line(aes(x = log2mf, y = density, group = condition, color = condition)) +
  scale_colour_manual(values = c("A" = "#f0f0f0",
                                 "B" = "#7bccc4",
                                 "C" = "#e34a33",
                                 "D" = "#43a2ca",
                                 "E" = "#b30000",
                                 "F" = "#0868ac",
                                 "G" = "#fc8d59",
                                 "H" = "#636363")) +
  xlab('Log2 Mosaic Fraction') +
  ylab('Density (Adjusted and Normalized)') +
  theme_pubr() +
  theme(legend.position = "right")

plot_smooth = ggplot(plot_tb) +
  geom_line(aes(x = -log2mf, y = density, group = condition, color = condition)) +
  scale_colour_manual(values = c("A" = "#fc8d59",
                                 "B" = "#7bccc4",
                                 "C" = "#e34a33",
                                 "D" = "#43a2ca",
                                 "E" = "#b30000",
                                 "F" = "#0868ac",
                                 "G" = "#fc8d59",
                                 "H" = "#636363")) +
  xlab('Time (Days)') +
  ylab('Mosaic Fraction Frequency') +
  ggtitle('Time Scaled Recovered Signal') +
  theme_pubr() +
  theme(legend.position = "right")

plot + plot_smooth

plot_tb_reg_id = nest(plot_tb, data = -ID)
plot_tb_id = nest(plot_tb_smooth, data = -ID)
plot_tb_id = left_join(plot_tb_id, nested_density_id_cum_mf$nested_density %>% select(ID, avg_density))
plot_tb_id$table = map(plot_tb_id$avg_density, function(density) {
  density$table
})
plot_tb_id$condition_max_mf = map(1:nrow(plot_tb_id), function(x) NULL)


#pdf(file = '~/Documents/R/MF_Signal_Simulation/plot_per_id_with_max_mf_sparse_IDs_samples_filtered.pdf', onefile = T)
for (i in 1:nrow(plot_tb_id)) {
  message(i)
  plot_tb = plot_tb_reg_id$data[[i]]
  plot_tb_smooth = plot_tb_id$data[[i]]
  # plot = ggplot(plot_tb) +
  #   geom_line(aes(x = log2mf, y = density, group = condition, color = condition)) +
  #   scale_colour_manual(values = c("A" = "#f0f0f0",
  #                                  "B" = "#7bccc4",
  #                                  "C" = "#e34a33",
  #                                  "D" = "#43a2ca",
  #                                  "E" = "#b30000",
  #                                  "F" = "#0868ac",
  #                                  "G" = "#fc8d59",
  #                                  "H" = "#636363")) +
  #   xlab('Log2 Mosaic Fraction') +
  #   ylab('Density (Adjusted and Normalized)') +
  #   theme_pubr() +
  #   theme(legend.position = "right") +
  #   ggtitle(plot_tb_reg_id$ID[[i]])
  #
  # plot_smooth = ggplot(plot_tb_smooth) +
  #   geom_line(aes(x = log2mf, y = smooth_density, group = condition, color = condition)) +
  #   scale_colour_manual(values = c("A" = "#f0f0f0",
  #                                  "B" = "#7bccc4",
  #                                  "C" = "#e34a33",
  #                                  "D" = "#43a2ca",
  #                                  "E" = "#b30000",
  #                                  "F" = "#0868ac",
  #                                  "G" = "#fc8d59",
  #                                  "H" = "#636363")) +
  #   xlab('Log2 Mosaic Fraction') +
  #   ylab('Density (Adjusted and Normalized)') +
  #   theme_pubr() +
  #   theme(legend.position = "right")
  #
  # plot_cum_mf = ggplot(data = plot_tb_id$table[[i]], aes(x = bin, y = count, color = 'cum_mf')) +
  #   geom_line()

  plot_tb_id$condition_max_mf[[i]] = plot_tb_smooth %>% group_by(condition) %>% summarise(max_mf = unique(max_mf))
  max_mf_plot = ggplot(data = plot_tb_id$condition_max_mf[[i]], aes(x = condition, y = max_mf)) +
    geom_point()

  print((plot + plot_smooth) / max_mf_plot)
}
dev.off()

plot_tb_id_max_mf = plot_tb_id %>% select(ID, condition_max_mf) %>% unnest(condition_max_mf)
ggplot(data = plot_tb_id_max_mf, aes(x = condition, y = -max_mf, fill = condition)) +
  scale_fill_manual(values = c("A" = "#fc8d59",
                                 "B" = "#7bccc4",
                                 "C" = "#e34a33",
                                 "D" = "#43a2ca",
                                 "E" = "#b30000",
                                 "F" = "#0868ac",
                                 "G" = "#fc8d59",
                                 "NC" = "#636363")) +
  geom_boxplot() + theme_pubr() + xlab('Treatment Group') +
  ylab('Negative Log2 Mosaic Fraction') + labs(fill = 'Treatment Group') +
  theme(legend.position = 'right') +
  ggtitle('Distribution of Most Frequently Occuring Mosaic Fraction Across Treatment Groups')

late = plot_tb_id_max_mf$max_mf[plot_tb_id_max_mf$condition %in% c('A', 'C', 'E')]
early = plot_tb_id_max_mf$max_mf[plot_tb_id_max_mf$condition %in% c('B', 'D', 'F')]
t.test(early, late, 'less')

#plotting
plot_id_sample_heatmap(raw_data, c('condition'))
plot_averaged_curves(nested_id_sample_mutated_filtered_id_condition,
                     group_str = 'condition', color_str = 'condition', facet_str = 'ID',
                     color_values = c("A" = "#f0f0f0",
                                      "B" = "#7bccc4",
                                      "C" = "#e34a33",
                                      "D" = "#43a2ca",
                                      "E" = "#b30000",
                                      "F" = "#0868ac",
                                      "G" = "#fc8d59",
                                      "H" = "#636363"))
plot_tb = nested_id_sample_mutated_filtered_id_condition$nested_density %>% select(ID, condition, avg_density)
plot_tb$table = map(plot_tb$avg_density, function(density) {
  density$table
})
plot_tb = unnest(plot_tb, cols = table)
ggplot(data = plot_tb, aes(x = bin,
                                  y = count,
                                  group = condition,
                                  color = condition)) +
  scale_colour_manual(values = c("A" = "#f0f0f0",
                                 "B" = "#7bccc4",
                                 "C" = "#e34a33",
                                 "D" = "#43a2ca",
                                 "E" = "#b30000",
                                 "F" = "#0868ac",
                                 "G" = "#fc8d59",
                                 "H" = "#636363")) +
  geom_line() +
  facet_wrap(~ID) +
  theme(legend.position = "right")

#testing
raw <- raw_data
parent_table <- parent
#all_samples_sequences_one_mouse <- all_samples_sequences[all_samples_sequences$mouse_number == '04',]
all_samples_sequences <- all_samples_sequences %>% select(sample, mouse_number, tissue, ID, sequence, reads)
all_samples_sequences$mutation_annotation <- 'observed'
save(all_samples_sequences, file = './all_mice_raw.rda')

colnames(indelphi_table)[2] <- 'sequence'
colnames(indelphi_table)[4] <- 'probability'

mf_vec = raw_data$data$mosaic_fraction[1:10]
mf_vec2 = raw_data$data$mosaic_fraction[11:20]
bin_size <- 1
breaks <- c(-Inf, seq(-10, 0, by = bin_size))
density1 = calculate_density(mf_vec, breaks)
density2 = calculate_density(mf_vec2, breaks)
density_list = list(density1, density2)
avg_density = average_density(density_list)
by_var = c('ID', 'sample')
nested_data = nested_id_sample_mutated
nested_data_filtered = nested_id_sample_mutated_filtered
