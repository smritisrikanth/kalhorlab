library(tidyverse)
library(ggpubr)
library(tidyr)
library(purrr)
library("readxl")
library(patchwork)

setwd('~/Documents/R/MosaicFractionAnalysis/')
devtools::load_all()

raw_data = load_samples('~/Documents/R/Zijie/filteredpairs_Z/')
parent = load_parent('~/Documents/R/Zijie/filteredpairs_Z/PB7-founder_filteredpairs.txt')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)
raw_data = filter_unmutated(raw_data)
raw_data$data$condition = substr(raw_data$data$sample, 10, 10)
raw_data$data = raw_data$data[raw_data$data$condition != 'A',]
raw_data$data$condition[raw_data$data$condition == 'G'] = 'A'
raw_data$data$condition[raw_data$data$condition == 'H'] = 'NC'
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
recur_mut_point001 <- read_excel("~/Documents/R/Zijie/recurring mutation/point001_percent_recurring_mutation.xlsx")
raw_data$data = raw_data$data[!(raw_data$data$sequence %in% recur_mut_point001$`0`),]
raw_data$data = raw_data$data[raw_data$data$log2MF < 0,]
raw_data$data = raw_data$data[!(raw_data$data$ID %in% sparse_ids$ID),]
raw_data$data = raw_data$data[!(raw_data$data$sample %in% sparse_samples$sample),]
nested_data_filtered = nest_raw(raw_data, c('ID', 'sample', 'condition'))
breaks <- c(-Inf, seq(-10, 0, by = 0.5))
bin_sizes = rev(diff(breaks))
bin_sizes[is.infinite(bin_sizes)] = min(bin_sizes)
nested_data = append_density(nested_data, breaks)
nested_data_filtered = append_density(nested_data_filtered, breaks)
nested_data_filtered = adjust_density(nested_data_filtered, nested_data, c('ID', 'sample', 'condition'))
nested_data = normalize_density_list(nested_data)
nested_data_filtered = normalize_density_list(nested_data_filtered)
nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)
nested_density_filtered = average_density_by(nested_data_filtered, c('condition'), normalized = TRUE)

plot_tb = nested_density$nested_density %>% select(condition, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    #density$table$count = density$table$count/sum(density$table$count)/bin_sizes
    tibble(log2mf = density$table$bin, density = density$table$count)
  })) %>%
  select(condition, table) %>%
  unnest(cols = table)

plot_tb_filtered = nested_density_filtered$nested_density %>% select(condition, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    #density$table$count = density$table$count/sum(density$table$count)/bin_sizes
    tibble(log2mf = density$table$bin, density = density$table$count)
  })) %>%
  select(condition, table) %>%
  unnest(cols = table)

plot_tb_smooth = nested_density_filtered$nested_density %>% select(condition, avg_density) %>%
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


plot = ggplot(plot_tb) +
  geom_line(aes(x = -log2mf, y = density, group = condition, color = condition)) +
  scale_colour_manual(values = c("A" = "#fc8d59",
                                 "B" = "#7bccc4",
                                 "C" = "#e34a33",
                                 "D" = "#43a2ca",
                                 "E" = "#b30000",
                                 "F" = "#0868ac",
                                 "G" = "#fc8d59",
                                 "NC" = "#636363")) +
  xlab('') +
  ylab('') +
  ggtitle('Time Scaled Recovered Signal') +
  theme_pubr() +
  theme(legend.position = "right") + labs(colour = 'Treatment Group')

plot_filtered = ggplot(plot_tb_filtered) +
  geom_line(aes(x = -log2mf, y = density, group = condition, color = condition)) +
  scale_colour_manual(values = c("A" = "#fc8d59",
                                 "B" = "#7bccc4",
                                 "C" = "#e34a33",
                                 "D" = "#43a2ca",
                                 "E" = "#b30000",
                                 "F" = "#0868ac",
                                 "G" = "#fc8d59",
                                 "NC" = "#636363")) +
  xlab('Time (Days)') +
  ylab('Mosaic Fraction Frequency') +
  ggtitle('Filtered') +
  theme_pubr() +
  theme(legend.position = "right") + labs(colour = 'Treatment Group')

plot_smooth = ggplot(plot_tb_smooth) +
  geom_line(aes(x = -log2mf, y = smooth_density, group = condition, color = condition)) +
  # geom_vline(aes(xintercept = -max_mf, group = condition, color = condition)) +
  scale_colour_manual(values = c("A" = "#fc8d59",
                                 "B" = "#7bccc4",
                                 "C" = "#e34a33",
                                 "D" = "#43a2ca",
                                 "E" = "#b30000",
                                 "F" = "#0868ac",
                                 "G" = "#fc8d59",
                                 "NC" = "#636363")) +
  xlab('Time (Days)') +
  ylab('') +
  ggtitle('Smoothed') +
  theme_pubr() +
  theme(legend.position = "right") + labs(colour = 'Treatment Group')

((plot + ylab('MF Frequency') + ggtitle('Unfiltered') + theme(legend.position = 'right',
                                                             axis.text = element_text(size = 6),
                                                             axis.title = element_text(size = 8),
                                                             legend.title = element_text(size = 8),
                                                             legend.text = element_text(size = 6),
                                                             plot.title = element_text(size = 9))) /
  (plot_filtered + ylab('MF Frequency') + theme(legend.position = 'right',
                                                  axis.text = element_text(size = 6),
                                                  axis.title = element_text(size = 8),
                                                  legend.title = element_text(size = 8),
                                                  legend.text = element_text(size = 6),
                                                  plot.title = element_text(size = 9)))) +
  (plot_smooth + ylab('MF Frequency') + theme(legend.position = 'right',
                                              axis.text = element_text(size = 6),
                                              axis.title = element_text(size = 8),
                                              legend.title = element_text(size = 8),
                                              legend.text = element_text(size = 6),
                                              plot.title = element_text(size = 9))) +
  plot_layout(guides = 'collect')

plot_smooth + ylab('MF Frequency') + theme(legend.position = 'right',
                                           axis.text = element_text(size = 6),
                                           axis.title = element_text(size = 8),
                                           legend.title = element_text(size = 8),
                                           legend.text = element_text(size = 6),
                                           plot.title = element_text(size = 9))

expected_mf = tibble(condition = unique(plot_tb_smooth$condition),
                     time = list(c(8:10), c(0:4), c(6:10), c(0:2), c(8:10), c(0:2), c(0:10)),
                     value = c(1, 1.1, 1.1, 2, 2, 1, 0))
expected_mf$condition = c('Late Low', 'Early Long', 'Late Long', 'Early High', 'Late High', 'Early Low', 'NC')
plot_tb_truth = unnest(expected_mf, cols = time)
truth_plot = ggplot(data = plot_tb_truth,
                    aes(x = time, y = value, group = condition, color = condition)) +
  geom_line() +
  scale_colour_manual(values = c("Early Low" = "#fc8d59",
                                 "Late Low" = "#7bccc4",
                                 "Early Long" = "#e34a33",
                                 "Late Long" = "#43a2ca",
                                 "Early High" = "#b30000",
                                 "Late High" = "#0868ac",
                                 "Early Low" = "#fc8d59",
                                 "NC" = "#636363")) + theme_pubr() +
  xlab('Time (Days)') + ylab('Relative Expected MF Frequency') + labs(colour = 'Treatment Group') +
  theme(legend.position = 'right',
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 0),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 9))

plot / plot_filtered / plot_smooth + plot_layout(guides = 'collect')

(plot_smooth + ylab('MF Frequency') + theme(legend.position = 'right',
                                           axis.text = element_text(size = 6),
                                           axis.title = element_text(size = 8),
                                           legend.title = element_text(size = 8),
                                           legend.text = element_text(size = 6),
                                           plot.title = element_text(size = 9))) /
  truth_plot + plot_layout(guides = 'collect')

png(filename = './z_smooth_truth.png', width = 3.5, height = 4, units = "in", res = 500)
print((plot_smooth + ylab('MF Frequency') + theme(legend.position = 'none',
                                                  axis.text = element_text(size = 6),
                                                  axis.title = element_text(size = 8),
                                                  legend.title = element_text(size = 8),
                                                  legend.text = element_text(size = 6),
                                                  plot.title = element_text(size = 9)) +
         ggtitle('Smoothed Reconstructed')) /
        truth_plot + ggtitle('Input') + plot_layout(guides = 'collect'))
dev.off()

png(filename = './z_unfiltered_filtered.png', width = 3.5, height = 4, units = "in", res = 500)
print((plot + ylab('MF Frequency') + ggtitle('Recurrent Mutations Not Filtered') + theme(legend.position = 'none',
                                                                   axis.text = element_text(size = 6),
                                                                   axis.title = element_text(size = 8),
                                                                   legend.title = element_text(size = 8),
                                                                   legend.text = element_text(size = 6),
                                                                   plot.title = element_text(size = 9))) /
        (plot_filtered + ylab('MF Frequency') + theme(legend.position = 'none',
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 8),
                                                      legend.title = element_text(size = 8),
                                                      legend.text = element_text(size = 6),
                                                      plot.title = element_text(size = 9)) +
           ggtitle('0.001% Recurrent Mutations Filtered')) +
        plot_layout(guides = 'collect'))
dev.off()

pdf(file = 'group_D_filtering.pdf', onefile = T)
dev.off()
print(plot_smooth)





