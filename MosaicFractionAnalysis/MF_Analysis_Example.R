source("MF_Functions.R")
source("MF_util.R")

library(tidyverse)
library(ggpubr)
library(tidyr)
library(purrr)

setwd('~/Documents/R/Zijie/filteredpairs_Z')
parent = load_parent('PB7-founder_filteredpairs.txt')
breaks <- c(-Inf, seq(-10, 0, by = 1))

raw_data = load_samples('~/Documents/R/Zijie/filteredpairs_Z/')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)

raw_data = compute_mosaic_fraction(raw_data)
raw_data$data$condition = substr(raw_data$data$sample, 10, 10)
nested_id_sample = nest_raw(raw_data, c('ID', 'sample', 'condition'))
nested_id_sample = compute_mutation_information(nested_id_sample)
raw_data = filter_unmutated(raw_data)
nested_id_sample_mutated = nest_raw(raw_data, c('ID', 'sample', 'condition'))
nested_id_sample_mutated = append_density(nested_id_sample_mutated, breaks)
nested_id_sample_density = average_density_by(nested_id_sample_mutated, c('condition'))

raw_data_filtered = structure(list(data = raw_data$data[!(raw_data$data$sequence %in% recur_mut_one$allele), ]), class = 'raw')
nested_id_sample_mutated_filtered = nest_raw(raw_data_filtered, c('ID', 'sample', 'condition'))
nested_id_sample_mutated_filtered = append_density(nested_id_sample_mutated_filtered, breaks)
nested_id_sample_mutated_filtered = adjust_density(nested_id_sample_mutated_filtered, nested_id_sample_mutated, c('ID', 'sample', 'condition'))
nested_id_sample_mutated_filtered = normalize_density_list(nested_id_sample_mutated_filtered)
nested_id_sample_mutated_filtered_id_condition = average_density_by(nested_id_sample_mutated_filtered, c('ID', 'condition'))

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
