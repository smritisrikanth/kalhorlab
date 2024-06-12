library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('~/Documents/qfm2')
devtools::load_all()
setwd('~/Documents/R/MF_Signal_Simulation/')

load('phy.rda')
#somatic mutations



param_tb = tibble(sim = 1:100)
param_tb$filename = paste0('mf_tables_', param_tb$sim, '.rda')
param_tb$mf_table = map(param_tb$filename, function(fn) {
  load(paste0('./three_cell_somatic_mut_ss5000/', fn))
  mf_table
})
# param_tb$mf_to_time_tb = map(param_tb$filename, function(fn) {
#   load(paste0('./three_cell_somatic_mut_ss5000/', fn))
#   mf_to_time_tb
# })

param_tb_unnest = unnest(param_tb, cols = mf_table) %>%
  select(sim, tissue = CellType, sequence, mosaic_fraction)
param_tb_unnest$sample = paste0(param_tb_unnest$sim, param_tb_unnest$tissue)
param_tb_unnest$ID = 'ID'

raw_data = structure(list(data = param_tb_unnest), class = 'raw')
raw_data$data = filter(raw_data$data, mosaic_fraction > 0)
raw_data = annotate_mut_node(raw_data, phy)
#raw_data$data = raw_data$data[raw_data$data$mosaic_fraction < 1,]
nested_data = nest_raw(raw_data, c('sim', 'ID', 'sample', 'tissue'))
breaks <- c(-Inf, seq(-13.5, -0.5, by = 1))
nested_data = append_density(nested_data, breaks)
nested_data = normalize_density_list(nested_data, apply_cum_adjustment = FALSE)
nested_data$nested$condition = nested_data$nested$tissue
param_tb$lm = map(param_tb$mf_to_time_tb, function(tb) {
  lm(data = tb, formula = node_time ~ log2mf)
})
nested_data$nested = left_join(nested_data$nested, select(param_tb, c(sim, lm)), by = 'sim')
nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)
#plot_averaged_curves(nested_density, group_str = 'condition', color_str = 'condition')

#fit a separate linear model for each node assignment, convert to time separately, and then merge again
lm = lm(data = param_tb$mf_to_time_tb[[1]], formula = node_time ~ log2mf)
plot_tb = nested_density$nested_density %>% select(condition, avg_density)
plot_tb$table = map(plot_tb$avg_density, function(density) {
  density$table
})
plot_tb = unnest(plot_tb, cols = table)
plot_tb$time = predict(lm, newdata = list(log2mf = plot_tb$bin))
ggplot(data = plot_tb) +
  geom_line(aes(x = time, y = count/sum(count), group = condition, color = condition)) +
  geom_line(aes(x = time, y = signal_func(time)/sum(signal_func(time)), color = 'truth')) +
  ylab('Mosaic Fraction Density') +
  xlab('time (days)')

plot_tb = nested_density$nested_density %>% select(condition, avg_density, density_list) %>%
  unnest(cols = density_list) %>% select(condition, avg_density, lm)
plot_tb$avg_density = map2(plot_tb$avg_density, plot_tb$lm, function(density, lm) {
  density$table$time = predict(lm, newdata = list(log2mf = density$table$bin))
  density
})
time_vec = average_density(plot_tb$avg_density, colname = 'time')$table$count
plot_tb$table = map(plot_tb$avg_density, function(density) {
  density$table$time = time_vec
  density$table
})
plot_tb = unnest(plot_tb, cols = table)
#plot_tb$time = predict(lm, newdata = list(log2mf = plot_tb$bin))
plot_tb$truth = signal_func_hgRNA(plot_tb$time)
ggplot(data = plot_tb) +
  #geom_line(aes(x = time, y = truth/sum(truth), color = 'truth')) +
  geom_line(aes(x = time, y = count/sum(count), group = condition, color = condition)) +
  theme(legend.position = "right")


#hgRNA mutations
param_tb = tibble(sample = 1:100)
param_tb$filename = paste0('mf_tables_', param_tb$sample, '.rda')
param_tb$mut_frac_mat = map(param_tb$filename, function(fn) {
  if (file.exists(paste0('./one_cell_hgRNA_mut_ss5000/', fn))) {
    load(paste0('./one_cell_hgRNA_mut_ss5000/', fn))
    mut_frac_mat
  } else {
    NULL
  }
})
param_tb$mf_to_time_tb = map(param_tb$filename, function(fn) {
  if (file.exists(paste0('./one_cell_hgRNA_mut_ss5000/', fn))) {
    load(paste0('./one_cell_hgRNA_mut_ss5000/', fn))
    mf_to_time_tb
  } else {
    NULL
  }
})
#param_tb = param_tb[-c(32,74,86),]

param_tb_unnest = unnest(param_tb, cols = mut_frac_mat) %>%
  select(ID, sample, sequence, mosaic_fraction)
raw_data = structure(list(data = param_tb_unnest), class = 'raw')
raw_data$data$probability = map2_dbl(raw_data$data$ID, raw_data$data$sequence, function(id, seq) {
  mut_p$recur_vec_list[[id]][[seq]]
})
pdf("./prob_cutoffs_total_avg_se_no_cum_adj.pdf", onefile = T)
i = -13
for (i in -14:0) {
  message(i)
  raw_data = structure(list(data = param_tb_unnest), class = 'raw')
  raw_data$data$probability = map2_dbl(raw_data$data$ID, raw_data$data$sequence, function(id, seq) {
    mut_p$recur_vec_list[[id]][[seq]]
  })
  nested_data = nest_raw(raw_data, c('ID', 'sample'))
  raw_data$data = raw_data$data[raw_data$data$probability < 2^i,]
  nested_data_filtered = nest_raw(raw_data, c('ID', 'sample'))
  breaks <- c(-Inf, seq(-13.5, -0.5, by = 1), 1.0)
  nested_data = append_density(nested_data, breaks)
  nested_data_filtered = append_density(nested_data_filtered, breaks)
  nested_data_filtered = adjust_density(nested_data_filtered, nested_data, c('ID', 'sample'))
  nested_data = nested_data_filtered
  nested_data = normalize_density_list(nested_data, apply_cum_adjustment = FALSE)
  nested_data_adjusted = normalize_density_list(nested_data, apply_cum_adjustment = TRUE)
  #nested_data$nested$condition = paste0(nested_data$nested$sample, nested_data$nested$ID)
  #nested_data$nested$condition = as.factor(nested_data$nested$sample %% 10)
  nested_data$nested$condition = 'all'
  nested_data_adjusted$nested$condition = 'all'
  param_tb$lm = map(param_tb$mf_to_time_tb, function(tb) {
    lm(data = tb, formula = node_time ~ log2mf)
  })
  nested_data$nested = left_join(nested_data$nested, select(param_tb, c(sample, lm)), by = 'sample')
  nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)
  nested_density_adjusted = average_density_by(nested_data_adjusted, c('condition'), normalized = TRUE)
  adjusted = nested_density_adjusted$nested_density$avg_density[[1]]$table$count
  #nested_density2 = average_density_by(nested_data, c('condition2'), normalized = TRUE)
  
  plot_tb = nested_density$nested_density %>% select(condition, avg_density, density_list) %>%
    unnest(cols = density_list) %>% select(condition, avg_density, 'lm')
  s = plot_tb$avg_density[[1]]$se
  plot_tb$avg_density = map2(plot_tb$avg_density, plot_tb$lm, function(density, lm) {
    density$table$time = predict(lm, newdata = list(log2mf = density$table$bin))
    density
  })
  time_vec = average_density(plot_tb$avg_density, colname = 'time')$table$count
  plot_tb$table = map(plot_tb$avg_density, function(density) {
    density$table$time = time_vec
    density$table$se = s
    density$table$adjusted = adjusted
    #density$table$total = nested_density2$nested_density$avg_density[[1]]$table$count
    density$table
  })
  plot_tb = unnest(plot_tb, cols = table)
  #plot_tb$time = predict(lm, newdata = list(log2mf = plot_tb$bin))
  plot_tb$truth = signal_func_hgRNA(plot_tb$time)
  print(ggplot(data = plot_tb) +
    #geom_ribbon(aes(x = time, ymin = count/sum(count) - se/sum(count), ymax = count/sum(count) + se/sum(count)), fill = 'gray') +
    #geom_line(aes(x = time, y = truth/sum(truth), color = 'truth')) +
    geom_line(aes(x = time, y = adjusted, group = condition, color = 'adjusted')) +
    geom_line(aes(x = time, y = count, group = condition, color = 'unadjusted')) +
    #geom_line(aes(x = time, y = adjusted, group = condition, color = 'total average')) +
    #geom_line(aes(x = time, y = total/sum(total), color = 'average across all samples')) +
    theme(legend.position = "right") +
    ylab('Mosaic Fraction Density') +
    xlab('time (days)') +
    ggtitle('Mosaic Fraction Density With and Without Cumulative Adjustment'))
}
dev.off()
raw_data$data = raw_data$data[raw_data$data$probability < 0.0001,]
nested_data = nest_raw(raw_data, c('ID', 'sample'))
breaks <- c(-Inf, seq(-13.5, -0.5, by = 1), 1.0)
nested_data = append_density(nested_data, breaks)
nested_data = normalize_density_list(nested_data, apply_cum_adjustment = TRUE)
#nested_data$nested$condition = paste0(nested_data$nested$sample, nested_data$nested$ID)
#nested_data$nested$condition = as.factor(nested_data$nested$sample %% 10)
nested_data$nested$condition = 'all'
param_tb$lm = map(param_tb$mf_to_time_tb, function(tb) {
  lm(data = tb, formula = node_time ~ log2mf)
})
nested_data$nested = left_join(nested_data$nested, select(param_tb, c(sample, lm)), by = 'sample')
nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)


plot_tb = nested_density$nested_density %>% select(condition, avg_density, density_list) %>%
  unnest(cols = density_list) %>% select(condition, avg_density, 'lm')
plot_tb$avg_density = map2(plot_tb$avg_density, plot_tb$lm, function(density, lm) {
  density$table$time = predict(lm, newdata = list(log2mf = density$table$bin))
  density
})
time_vec = average_density(plot_tb$avg_density, colname = 'time')$table$count
plot_tb$table = map(plot_tb$avg_density, function(density) {
  density$table$time = time_vec
  density$table
})
plot_tb = unnest(plot_tb, cols = table)
#plot_tb$time = predict(lm, newdata = list(log2mf = plot_tb$bin))
plot_tb$truth = signal_func_hgRNA(plot_tb$time)
ggplot(data = plot_tb) +
  geom_line(aes(x = time, y = truth/sum(truth), color = 'truth')) +
  geom_line(aes(x = time, y = count/sum(count), group = condition, color = condition)) +
  theme(legend.position = "right") +
  ylab('Mosaic Fraction Density') +
  xlab('time (days)')

#three cell hgRNA
param_tb = tibble(sim = 1:100)
param_tb$filename = paste0('mf_tables_', param_tb$sim, '.rda')
param_tb$mut_frac_mat = map(param_tb$filename, function(fn) {
  if (file.exists(paste0('./three_cell_hgRNA_mut_ss5000/', fn))) {
    load(paste0('./three_cell_hgRNA_mut_ss5000/', fn))
    mut_frac_mat
  } else {
    NULL
  }
})
param_tb$mf_to_time_tb = map(param_tb$filename, function(fn) {
  if (file.exists(paste0('./three_cell_hgRNA_mut_ss5000/', fn))) {
    load(paste0('./three_cell_hgRNA_mut_ss5000/', fn))
    mf_to_time_tb
  } else {
    NULL
  }
})
#param_tb = param_tb[-c(43,45),]
param_tb_unnest = unnest(param_tb, cols = mut_frac_mat) %>%
  select(ID, sim, tissue = CellType, sequence, mosaic_fraction)
param_tb_unnest$sample = paste0(param_tb_unnest$sim, param_tb_unnest$tissue)

raw_data = structure(list(data = param_tb_unnest), class = 'raw')
raw_data$data = filter(raw_data$data, mosaic_fraction > 0)
raw_data$data$probability = map2_dbl(raw_data$data$ID, raw_data$data$sequence, function(id, seq) {
  mut_p$recur_vec_list[[id]][[seq]]
})
raw_data$data = raw_data$data[raw_data$data$probability < 0.001,]
raw_data = annotate_mut_node(raw_data, phy)
nested_data = nest_raw(raw_data, c('ID', 'sample', 'tissue'))
breaks <- c(-Inf, seq(-13.5, -0.5, by = 1))
nested_data = append_density(nested_data, breaks)
nested_data = normalize_density_list(nested_data, apply_cum_adjustment = TRUE)
nested_data$nested$condition = nested_data$nested$tissue
nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)
#plot_averaged_curves(nested_density, group_str = 'condition', color_str = 'condition')

#fit a separate linear model for each node assignment, convert to time separately, and then merge again
lm = lm(data = param_tb$mf_to_time_tb[[1]], formula = node_time ~ log2mf)
plot_tb = nested_density$nested_density %>% select(condition, avg_density)
plot_tb$table = map(plot_tb$avg_density, function(density) {
  density$table
})
plot_tb = unnest(plot_tb, cols = table)
plot_tb$time = predict(lm, newdata = list(log2mf = plot_tb$bin))
ggplot(data = plot_tb) +
  geom_line(aes(x = time, y = count/sum(count), group = condition, color = condition)) +
  geom_line(aes(x = time, y = signal_func_hgRNA(time)/sum(signal_func_hgRNA(time)), color = 'truth')) +
  ylab('Mosaic Fraction Density') +
  xlab('time (days)')

