library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('~/Documents/qfm2')
devtools::load_all()

mut_p = readRDS("./metadata//mut_p_marc1.rds")
mut_p$mut_rate_func = map(1:length(mut_p$mut_rate), function(i) {
  signal_func_hgRNA
})
mut_p$mut_rate = NULL

setwd("~/Documents/R/MosaicFractionAnalysis/R")
devtools::load_all()


setwd('~/Documents/R/MF_Signal_Simulation/')

load('phy.rda')

signal_func_hgRNA = function(x_vec) {
  val = 5e-1 * ((1/6)* pmin(pmax(0, x_vec-1), 4) + (-1/6) * pmax(0, x_vec-5))
  out_vec = pmax(val, 0)
  #out_vec[x_vec < 1] = 0
  out_vec
}
signal_func = function(x_vec) {
  val = 1e-9 * ((1/6)* pmin(pmax(0, x_vec-1), 4) + (-1/6) * pmax(0, x_vec-5))
  out_vec = pmax(val, 0) #+ 2*10^(-9)/3
  #out_vec[x_vec < 1] = 0
  out_vec
}

#somatic mutations

param_tb = tibble(sim = 1:100)
param_tb$filename = paste0('mf_table_', param_tb$sim, '.rda')
param_tb$mf_table = map(param_tb$filename, function(fn) {
  load(paste0('./three_cell_somatic_mut_ss5000_signal4/', fn))
  mf_table %>% select(-c(value,mut,cell)) %>% nest(data = -sequence) %>%
    mutate(avg_log2mf = map_dbl(data, function(data) {log2(mean(data$mosaic_fraction))})) %>% unnest(cols = data)
})
param_tb_unnest = unnest(param_tb, cols = mf_table) %>%
  select(sim, tissue = CellType, sequence, mosaic_fraction, avg_log2mf, node_time)
param_tb_unnest$sample = paste0(param_tb_unnest$sim, param_tb_unnest$tissue)
param_tb_unnest$ID = 'ID'
raw_data = structure(list(data = param_tb_unnest), class = 'raw')
raw_data$data = filter(raw_data$data, mosaic_fraction > 0)
raw_data = annotate_mut_node(raw_data, phy)
nested_data = nest_raw(raw_data, c('ID', 'sample', 'gr_node'))
breaks <- c(-Inf, seq(-13.5, -0.5, by = 1), 1.0)
nested_data = append_density(nested_data, breaks)
nested_data = normalize_density_list(nested_data, apply_cum_adjustment = FALSE)
nested_data$nested$condition = nested_data$nested$gr_node
nested_density = average_density_by(nested_data, c('condition'), raw = TRUE)

nested_density$nested_density$lm = map(nested_density$nested_density$density_list, function (list) {
  log2mf = as.vector(c())
  node_time = as.vector(c())
  for (tb in list$data) {
    log2mf = c(log2mf, tb$avg_log2mf)
    node_time = c(node_time, tb$node_time)
  }
  lm(formula = node_time ~ log2mf)
})
nested_density$nested_density$avg_density = map(nested_density$nested_density$avg_density, function(density) {
  normalize_density(density, apply_cum_adjustment = FALSE)
})
plot_tb = nested_density$nested_density %>% select(condition, avg_density, lm)
plot_tb$table = map2(plot_tb$avg_density, plot_tb$lm, function(density, lm) {
  density$table$time = predict(lm, newdata = list(log2mf = density$table$bin))
  density$table
})
matrix_mf = as.vector(c())
matrix_time = as.vector(c())
for (i in plot_tb$table) {
  matrix_mf = rbind(matrix_mf, i$count)
  matrix_time = rbind(matrix_time, i$time)
}
matrix_time_zero = matrix_time
matrix_time_zero[matrix_mf == 0] = 0
time_vec = colSums(matrix_time_zero)/colSums(matrix_mf != 0)
overall_mf_vec = colSums(matrix_mf)/colSums(matrix_mf != 0)
time_vec[is.na(time_vec)] = colMeans(matrix_time)[is.na(time_vec)]
overall_mf_vec[is.na(overall_mf_vec)] = 0
overall_mf_vec = colMeans(matrix_mf)
plot_tb$table = map(plot_tb$table, function(table) {
  table$time = time_vec
  table$overall = overall_mf_vec
  table
})
plot_tb = unnest(plot_tb, cols = table)
overall_density = structure(list(table = plot_tb %>% select(bin, count = overall, cum_mf = cum_mf)), class = 'mf_density')
overall_density = normalize_density(overall_density, apply_cum_adjustment = FALSE)
plot_tb$norm_count = overall_density$table$norm_count
ggplot(data = plot_tb) +
  geom_line(aes(x = time, y = count/sum(count), group = condition, color = condition)) +
  geom_line(aes(x = time, y = signal_func(time)/sum(signal_func(time)), color = 'truth')) +
  geom_line(aes(x = time, y = overall/sum(overall), color = 'raw averaged')) +
  geom_line(aes(x = time, y = norm_count/sum(norm_count), color = 'normalized after')) +
  ylab('Mosaic Fraction Density') +
  xlab('time (days)') +
  ggtitle('')

#one cell type hgRNA
param_tb = tibble(sample = 1:100)
param_tb$filename = paste0('mf_tables_', param_tb$sample, '.rda')
param_tb$mut_frac_mat = map(param_tb$filename, function(fn) {
  if (file.exists(paste0('./one_cell_hgRNA_mut_ss5000_signal3/', fn))) {
    load(paste0('./one_cell_hgRNA_mut_ss5000_signal3/', fn))
    mut_frac_mat  %>% nest(data = -sequence) %>%
      mutate(avg_log2mf = map_dbl(data, function(data) {log2(mean(data$mosaic_fraction))})) %>% unnest(cols = data)
  } else {
    NULL
  }
})

param_tb_unnest = unnest(param_tb, cols = mut_frac_mat) %>%
  select(ID, sample, sequence, mosaic_fraction, tissue = CellType, probability, avg_log2mf, node_time)
pdf('one_cell_hgRNA_prob_cutoffs_dist_normalized.pdf', onefile = T)
i = 0
dist = as.vector(c())
dist_unadj = as.vector(c())
cutoffs = as.vector(c())
corr = as.vector(c())
corr_unadj = as.vector(c())
for (i in -14:0) {
  message(i)
  raw_data = structure(list(data = param_tb_unnest), class = 'raw')
  lm = lm(data = raw_data$data, formula = node_time ~ avg_log2mf)
  nested_data = nest_raw(raw_data, c('ID', 'sample'))
  #i = -13
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
  nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)
  nested_density_adjusted = average_density_by(nested_data_adjusted, c('condition'), normalized = TRUE)
  adjusted = nested_density_adjusted$nested_density$avg_density[[1]]$table$count
  #nested_density2 = average_density_by(nested_data, c('condition2'), normalized = TRUE)
  
  plot_tb = nested_density$nested_density %>% select(condition, avg_density)
  #s = plot_tb$avg_density[[1]]$se
  plot_tb$table = map(plot_tb$avg_density, function(density) {
    density$table
  })
  plot_tb$table = map(plot_tb$avg_density, function(density) {
    #density$table$time = time_vec
    #density$table$se = s
    density$table$adjusted = adjusted
    #density$table$total = nested_density2$nested_density$avg_density[[1]]$table$count
    density$table
  })
  plot_tb = unnest(plot_tb, cols = table)
  plot_tb$time = predict(lm, newdata = list(avg_log2mf = plot_tb$bin))
  plot_tb$truth = signal_func_hgRNA(plot_tb$time)
  cutoffs = c(cutoffs, 2^i)
  dist = c(dist, sqrt(sum((plot_tb$truth/sum(plot_tb$truth) - plot_tb$adjusted/sum(plot_tb$adjusted))^2)))
  dist_unadj = c(dist_unadj, sqrt(sum((plot_tb$truth/sum(plot_tb$truth) - plot_tb$count/sum(plot_tb$count))^2)))
  corr = c(corr, cor(plot_tb$truth/sum(plot_tb$truth), plot_tb$adjusted/sum(plot_tb$adjusted)))
  corr_unadj = c(corr_unadj, cor(plot_tb$truth/sum(plot_tb$truth), plot_tb$count/sum(plot_tb$count)))
  print(ggplot(data = plot_tb) +
          #geom_ribbon(aes(x = time, ymin = count/sum(count) - se/sum(count), ymax = count/sum(count) + se/sum(count)), fill = 'gray') +
          geom_line(aes(x = time, y = truth/sum(truth), color = 'truth')) +
          geom_line(aes(x = time, y = adjusted/sum(adjusted), group = condition, color = 'adjusted')) +
          #geom_line(aes(x = time, y = count/sum(count), group = condition, color = 'unadjusted')) +
          #geom_line(aes(x = time, y = adjusted, group = condition, color = 'adjusted')) +
          #geom_line(aes(x = time, y = total/sum(total), color = 'average across all samples')) +
          theme(legend.position = "right") +
          ylab('Mosaic Fraction Density') +
          xlab('time (days)') +
          ggtitle(paste0('p < ', 2^i, '; dist = ', dist)))
}
dev.off()

plot_tb_dist = tibble(cutoffs, dist, dist_unadj)
ggplot(data = plot_tb_dist) +
  geom_line(aes(x = log2(cutoffs), y = dist, color = 'adjusted')) +
  geom_line(aes(x = log2(cutoffs), y = dist_unadj, color = 'unadjusted')) +
  theme(legend.position = "right") +
  ylab('Euclidean Distance') +
  xlab('Log2 Mutation Occurrence Probability Threshold') +
  ggtitle('Euclidean Distance Between True and Predicted Signal')

plot_tb_corr = tibble(cutoffs, corr, corr_unadj)
ggplot(data = plot_tb_dist) +
  geom_line(aes(x = log2(cutoffs), y = corr, color = 'adjusted')) +
  geom_line(aes(x = log2(cutoffs), y = corr_unadj, color = 'unadjusted')) +
  theme(legend.position = "right") +
  ylab('Pearson Correlation Coefficiet') +
  xlab('Log2 Mutation Occurrence Probability Threshold') +
  ggtitle('Pearson Correlation Between True and Predicted Signal')
  

raw_data = structure(list(data = param_tb_unnest), class = 'raw')
lm = lm(data = raw_data$data, formula = node_time ~ avg_log2mf)
nested_data = nest_raw(raw_data, c('ID', 'sample'))
#i = -13
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
nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)
nested_density_adjusted = average_density_by(nested_data_adjusted, c('condition'), normalized = TRUE)
adjusted = nested_density_adjusted$nested_density$avg_density[[1]]$table$count
#nested_density2 = average_density_by(nested_data, c('condition2'), normalized = TRUE)

plot_tb = nested_density$nested_density %>% select(condition, avg_density)
#s = plot_tb$avg_density[[1]]$se
plot_tb$table = map(plot_tb$avg_density, function(density) {
  density$table
})
plot_tb$table = map(plot_tb$avg_density, function(density) {
  #density$table$time = time_vec
  #density$table$se = s
  density$table$adjusted = adjusted
  #density$table$total = nested_density2$nested_density$avg_density[[1]]$table$count
  density$table
})
plot_tb = unnest(plot_tb, cols = table)
plot_tb$time = predict(lm, newdata = list(avg_log2mf = plot_tb$bin))
plot_tb$truth = signal_func_hgRNA(plot_tb$time)
print(ggplot(data = plot_tb) +
        #geom_ribbon(aes(x = time, ymin = count/sum(count) - se/sum(count), ymax = count/sum(count) + se/sum(count)), fill = 'gray') +
        #geom_line(aes(x = time, y = truth/sum(truth), color = 'truth')) +
        geom_line(aes(x = time, y = adjusted, group = condition, color = 'adjusted')) +
        #geom_line(aes(x = time, y = count, group = condition, color = 'unadjusted')) +
        #geom_line(aes(x = time, y = adjusted, group = condition, color = 'total average')) +
        #geom_line(aes(x = time, y = total/sum(total), color = 'average across all samples')) +
        theme(legend.position = "right") +
        ylab('Mosaic Fraction Density') +
        xlab('time (days)') +
        ggtitle(paste0('p < ', 2^i)))


#three cell hgRNA
param_tb = tibble(sim = 1:100)
param_tb$filename = paste0('mf_tables_', param_tb$sim, '.rda')
param_tb$mut_frac_mat = map(param_tb$filename, function(fn) {
  load(paste0('./three_cell_hgRNA_mut_ss5000_signal4/', fn))
  mut_frac_mat  %>% nest(data = -sequence) %>%
    mutate(avg_log2mf = map_dbl(data, function(data) {log2(mean(data$mosaic_fraction))})) %>% unnest(cols = data)
})
param_tb_unnest = unnest(param_tb, cols = mut_frac_mat) %>%
  select(sim, ID, tissue = CellType, sequence, mosaic_fraction, probability, avg_log2mf, node_time)
param_tb_unnest$sample = paste0(param_tb_unnest$sim, param_tb_unnest$tissue)
raw_data = structure(list(data = param_tb_unnest), class = 'raw')
raw_data$data = filter(raw_data$data, mosaic_fraction > 0)
raw_data = annotate_mut_node(raw_data, phy)
nested_data = nest_raw(raw_data, c('ID', 'sample', 'gr_node'))
raw_data_original = raw_data
pdf("./three_cell_hgRNA_prob_cutoffs_dist_normalized.pdf", onefile = T)
dist = as.vector(c())
cutoffs = as.vector(c())
i = -15
for (i in -10:-16) {
  j = i
  message(i)
  # raw_data = structure(list(data = param_tb_unnest), class = 'raw')
  # raw_data$data = filter(raw_data$data, mosaic_fraction > 0)
  # raw_data = annotate_mut_node(raw_data, phy)
  # nested_data = nest_raw(raw_data, c('ID', 'sample', 'gr_node'))
  #i = -10
  raw_data$data = raw_data$data[raw_data$data$probability < 2^i,]
  nested_data_filtered = nest_raw(raw_data, c('ID', 'sample', 'gr_node'))
  breaks <- c(-Inf, seq(-13.5, -0.5, by = 1), 1.0)
  nested_data = append_density(nested_data, breaks)
  nested_data_filtered = append_density(nested_data_filtered, breaks)
  nested_data_filtered = adjust_density(nested_data_filtered, nested_data, c('ID', 'sample', 'gr_node'))
  nested_data = nested_data_filtered
  nested_data = normalize_density_list(nested_data, apply_cum_adjustment = FALSE)
  nested_data$nested$condition = nested_data$nested$gr_node
  nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)
  
  nested_density$nested_density$lm = map(nested_density$nested_density$density_list, function (list) {
    log2mf = as.vector(c())
    node_time = as.vector(c())
    for (tb in list$data) {
      log2mf = c(log2mf, tb$avg_log2mf)
      node_time = c(node_time, tb$node_time)
    }
    lm(formula = node_time ~ log2mf)
  })
  
  plot_tb = nested_density$nested_density %>% select(condition, avg_density, lm)
  plot_tb$table = map2(plot_tb$avg_density, plot_tb$lm, function(density, lm) {
    density$table$time = predict(lm, newdata = list(log2mf = density$table$bin))
    density$table
  })
  matrix_mf = as.vector(c())
  matrix_time = as.vector(c())
  for (i in plot_tb$table) {
    matrix_mf = rbind(matrix_mf, i$count)
    matrix_time = rbind(matrix_time, i$time)
  }
  matrix_time_zero = matrix_time
  matrix_time_zero[matrix_mf == 0] = 0
  time_vec = colSums(matrix_time_zero)/colSums(matrix_mf != 0)
  overall_mf_vec = colSums(matrix_mf)/colSums(matrix_mf != 0)
  time_vec[is.na(time_vec)] = colMeans(matrix_time)[is.na(time_vec)]
  overall_mf_vec[is.na(overall_mf_vec)] = 0
  plot_tb$table = map(plot_tb$table, function(table) {
    table$time = time_vec
    table$overall = overall_mf_vec
    table
  })
  plot_tb = unnest(plot_tb, cols = table)
  #i = -10
  cutoffs = c(cutoffs, 2^j)
  dist = c(dist, sqrt(sum((signal_func_hgRNA(time_vec)/sum(signal_func_hgRNA(time_vec)) - overall_mf_vec/sum(signal_func_hgRNA(time_vec)))^2)))
  print(ggplot(data = plot_tb) +
    geom_line(aes(x = time, y = count/sum(count), group = condition, color = condition)) +
    geom_line(aes(x = time, y = signal_func_hgRNA(time)/sum(signal_func_hgRNA(time)), color = 'truth')) +
    geom_line(aes(x = time, y = overall/sum(overall), color = 'overall')) +
    ylab('Mosaic Fraction Density') +
    xlab('time (days)'))
    #ggtitle(paste0('p < ', 2^j, '; distance = ', dist)))
}
dev.off()
raw_data = structure(list(data = param_tb_unnest), class = 'raw')
raw_data$data = filter(raw_data$data, mosaic_fraction > 0)
raw_data = annotate_mut_node(raw_data, phy)
nested_data = nest_raw(raw_data, c('ID', 'sample', 'gr_node'))
i = -15
raw_data$data = raw_data$data[raw_data$data$probability < 2^i,]
nested_data_filtered = nest_raw(raw_data, c('ID', 'sample', 'gr_node'))
breaks <- c(-Inf, seq(-13.5, -0.5, by = 1), 1.0)
nested_data = append_density(nested_data, breaks)
nested_data_filtered = append_density(nested_data_filtered, breaks)
nested_data_filtered = adjust_density(nested_data_filtered, nested_data, c('ID', 'sample', 'gr_node'))
nested_data = nested_data_filtered
nested_data = normalize_density_list(nested_data, apply_cum_adjustment = TRUE)
nested_data$nested$condition = nested_data$nested$gr_node
nested_density = average_density_by(nested_data, c('condition'), raw = TRUE)

nested_density$nested_density$lm = map(nested_density$nested_density$density_list, function (list) {
  log2mf = as.vector(c())
  node_time = as.vector(c())
  for (tb in list$data) {
    log2mf = c(log2mf, tb$avg_log2mf)
    node_time = c(node_time, tb$node_time)
  }
  lm(formula = node_time ~ log2mf)
})

plot_tb = nested_density$nested_density %>% select(condition, avg_density, lm)
plot_tb$table = map2(plot_tb$avg_density, plot_tb$lm, function(density, lm) {
  density$table$time = predict(lm, newdata = list(log2mf = density$table$bin))
  density$table
})
matrix_mf = as.vector(c())
matrix_time = as.vector(c())
for (i in plot_tb$table) {
  matrix_mf = rbind(matrix_mf, i$count)
  matrix_time = rbind(matrix_time, i$time)
}
matrix_time_zero = matrix_time
matrix_time_zero[matrix_mf == 0] = 0
time_vec = colSums(matrix_time_zero)/colSums(matrix_mf != 0)
overall_mf_vec = colSums(matrix_mf)/colSums(matrix_mf != 0)
time_vec[is.na(time_vec)] = colMeans(matrix_time)[is.na(time_vec)]
overall_mf_vec[is.na(overall_mf_vec)] = 0
plot_tb$table = map(plot_tb$table, function(table) {
  table$time = time_vec
  table$overall = overall_mf_vec
  table
})
plot_tb = unnest(plot_tb, cols = table)
i = -15
dist = sqrt(sum((signal_func_hgRNA(time_vec) - overall_mf_vec)^2))
ggplot(data = plot_tb) +
  geom_line(aes(x = time, y = count, group = condition, color = condition)) +
  geom_line(aes(x = time, y = signal_func_hgRNA(time), color = 'truth')) +
  geom_line(aes(x = time, y = overall, color = 'overall')) +
  ylab('Mosaic Fraction Density') +
  xlab('time (days)') +
  ggtitle(paste0('p < ', 2^i, '; distance = ', dist))

