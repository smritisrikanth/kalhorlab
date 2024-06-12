library(purrr)
library(ggpubr)
library(ggplot2)
library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('~/Documents/qfm2')
devtools::load_all()

setwd('~/Documents/R/MF_Signal_Simulation/')
#load('param_tb.rda')

param_tb = tibble(sample = 1:100)

#param_tb$filename = paste0('count_graph_', param_tb$sample, '.rda_1.rda')
param_tb$filename = paste0('mf_vec_', param_tb$sample, '.rda')

#fn = param_tb$filename[[1]]
param_tb$sim_data = map(param_tb$filename, function(fn) {
  load(paste0('./output2_ab0.15_ss1000/', fn))
  tibble(sequence = names(mf_vec),
         mosaic_fraction = mf_vec)
})


# param_tb$sim_data = map2(param_tb$filename, param_tb$num_sim, function(fn, num) {
#   filename = paste0(fn, '_', num, '.rda')
#   load(paste0('./output2/', filename))
#   tibble(sequence = names(mf_vec),
#          mosaic_fraction = mf_vec)
# })
param_tb$ID = 'ID'

#param_tb$sample = paste0(param_tb$filename, '_', param_tb$num_sim)

param_tb = unnest(param_tb, cols = "sim_data")

param_tb = param_tb[param_tb$mosaic_fraction < 1,]


raw_data = structure(list(data = param_tb), class = 'raw')
nested_data = nest_raw(raw_data, c('ID', 'sample'))
#nested_data = structure(list(nested = nest(param_tb, data = -c(ID, sample, filename)), class = 'nested_data'))
breaks <- c(-Inf, seq(-13.5, -0.5, by = 1), 1)
nested_data = append_density(nested_data, breaks)
nested_data = normalize_density_list(nested_data, apply_cum_adjustment = FALSE)
#nested_data$nested$condition = substr(nested_data$nested$filename, 1, 13)
#nested_data$nested$condition = nested_data$nested$ID
nested_data$nested$condition = nested_data$nested$sample
nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)

#plotting
plot_averaged_curves(nested_density, group_str = 'condition', color_str = 'condition')



#second set of simulations

param_tb_2 = tibble(sample = 1:100)
param_tb_2$filename = paste0('mf_vec_', param_tb_2$sample, '.rda')
param_tb_2$sim_data = map(param_tb_2$filename, function(fn) {
  load(paste0('./output2_ab0.65_ss1000/', fn))
  tibble(sequence = names(mf_vec),
         mosaic_fraction = mf_vec)
})
param_tb_2$ID = 'ID'
param_tb_2 = unnest(param_tb_2, cols = "sim_data")
#param_tb_2 = param_tb_2[param_tb_2$mosaic_fraction < 1,]

raw_data = structure(list(data = param_tb_2), class = 'raw')
nested_data = nest_raw(raw_data, c('ID', 'sample'))
#nested_data = structure(list(nested = nest(param_tb, data = -c(ID, sample, filename)), class = 'nested_data'))
#breaks <- c(-Inf, seq(-13, 0, by = 1))
nested_data = append_density(nested_data, breaks)
nested_data = normalize_density_list(nested_data, apply_cum_adjustment = FALSE)
#nested_data$nested$condition = substr(nested_data$nested$filename, 1, 13)
#nested_data$nested$condition = nested_data$nested$ID
nested_data$nested$condition = nested_data$nested$sample %% 10
nested_density2 = average_density_by(nested_data, c('condition'), normalized = FALSE)

plot_averaged_curves(nested_density2, group_str = 'condition', color_str = 'condition')

#plotting two sets together
plot_tb = nested_density$nested_density %>% select(condition, avg_density1 = avg_density) %>% left_join(nested_density2$nested_density[,c('condition', 'avg_density')], by = 'condition')
plot_tb$table = map(plot_tb$avg_density1, function(density) {
  colnames(density$table) = c('bin', 'async')
  density$table
})
plot_tb$table2 = map(plot_tb$avg_density, function(density) {
  colnames(density$table) = c('bin2', 'synchronized')
  density$table
})
plot_tb = unnest(plot_tb, cols = c(table, table2))
plot_tb = plot_tb %>% select(condition, bin, async, synchronized)
#plot_tb = plot_tb %>% gather(key = 'division_rate', value = 'MF_density', -c('condition', 'bin'))
#plot_tb$group = paste0(plot_tb$condition, '_', plot_tb$division_rate)
# ggplot(data = plot_tb, aes(x = bin,
#                            y = count,
#                            group = group,
#                            color = simulation)) +
#   #geom_point() +
#   geom_line() +
#   theme(legend.position = 'right')

#true time vs log2mf scatter plot
#ggscatter(cell_mut_tb_sync, y = "node_time", x = "log2mf") + geom_smooth(method = "lm")
#ggscatter(cell_mut_tb_async, y = "node_time", x = "log2mf") + geom_smooth(method = "lm")
lm_sync = lm(data = cell_mut_tb_sync, formula = node_time ~ log2mf)
plot_tb$time_sync = predict(lm_sync, newdata = list(log2mf = plot_tb$bin))
lm_async = lm(data = cell_mut_tb_async, formula = node_time ~ log2mf)
plot_tb$time_async = predict(lm_async, newdata = list(log2mf = plot_tb$bin))
# ggplot(data = plot_tb, aes(x = time,
#                            y = MF_density,
#                            group = group,
#                            color = division_rate)) +
#   #geom_point() +
#   geom_line() +
#   theme(legend.position = 'right') +
#   geom_line(aes(x = time, y = signal_func(time)*10^9, color = 'truth'))

ggplot(data = plot_tb) +
  geom_line(aes(x = time_sync, y = synchronized/sum(synchronized), group = condition, color = 'synchronized')) +
  geom_line(aes(x = time_sync, y = signal_func(time_sync)/sum(signal_func(time_sync)), color = 'truth')) +
  geom_line(aes(x = time_async, y = async/sum(async), group = condition, color = 'async')) +
  theme(legend.position = 'right') +
  ylab('Mosaic Fraction Density') +
  xlab('time (days)') +
  ggtitle('Unajusted')

# ggplot(data = plot_tb) +
#   geom_line(aes(x = bin, y = async, group = condition, color = 'async')) +
#   geom_line(aes(x = bin, y = synchronized, group = condition, color = 'synchronized')) +
#   geom_line(aes(x = bin, y = signal_func(bin), color = 'truth')) +
#   theme(legend.position = 'right') +
#   ylab('Mosaic Fraction Density') +
#   xlab('bin') +
#   ggtitle('Adjusted') +
#   geom_vline(xintercept = -log2(2500))


# 
