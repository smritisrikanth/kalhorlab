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


raw_data = load_simulated_data('./three_cell_somatic_mut_ss5000_signal4/', 100, mat_name = 'mf_table')

#raw_data$data = filter(raw_data$data, mosaic_fraction > 0)
raw_data = annotate_mut_node(raw_data, phy)
nested_data = nest_raw(raw_data, c('ID', 'sample', 'gr_node'))
breaks <- c(-Inf, seq(-13.5, -0.5, by = 1), 1.0)
nested_data = append_density(nested_data, breaks)


#correct cum_mf values
node_order_list = list('TypeS'=1, 'TypeT'=2, 'TypeC'=3, 'TypeA'=4, 'TypeB'=5)
nested_data$nested$gr_node_num = map_dbl(nested_data$nested$gr_node, function(node) {node_order_list[[node]]})
#nested_data$nested = nested_data$nested[order(nested$gr_node_num),]
nested_sample = nest(nested_data$nested, outer_data = -sample)
nested_sample$outer_data = map(nested_sample$outer_data, function(outer_data) {
  outer_data = outer_data[order(outer_data$gr_node_num),]
  if (nrow(outer_data) == 1) {
    outer_data
  }
  cum_mf_total = outer_data$density[[1]]$table$cum_mf[length(outer_data$density[[1]]$table$cum_mf)]
  for (i in 2:nrow(outer_data)) {
    outer_data$density[[i]]$table$cum_mf = outer_data$density[[i]]$table$cum_mf + cum_mf_total
    cum_mf_total = outer_data$density[[i]]$table$cum_mf[length(outer_data$density[[1]]$table$cum_mf)]
  }
  outer_data
})
nested_data$nested = unnest(nested_sample, cols = outer_data)


#correct cum_mf
node_order_list = list('TypeS'=1, 'TypeT'=2, 'TypeC'=3, 'TypeA'=4, 'TypeB'=5)
nested_data = increment_cum_mf(nested_data, node_order_list)

#cumulative adjustment, 1/MF normalization, and averaging across node type
nested_data = normalize_density_list(nested_data, apply_cum_adjustment = FALSE)
nested_data$nested$condition = nested_data$nested$gr_node
nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)

#time conversion
nested_density$nested_density$lm = map(nested_density$nested_density$density_list, function (list) {
  log2mf = as.vector(c())
  node_time = as.vector(c())
  for (tb in list$data) {
    log2mf = c(log2mf, tb$avg_log2mf)
    node_time = c(node_time, tb$node_time)
  }
  lm(formula = node_time ~ log2mf)
})
plot_info = nested_density$nested_density %>% select(condition, avg_density, lm)
plot_info$table = map2(plot_info$avg_density, plot_info$lm, function(density, lm) {
  density$table$time = predict(lm, newdata = list(log2mf = density$table$bin))
  density$table
})

nested_density = convert_mf_to_time(nested_density)

#bind matrix with counts for all conditions along time_vec
doubling_time = 0.4
time_vec = seq(0,length(breaks)*doubling_time, doubling_time)
time_df = as.data.frame(time_vec)
for (i in 1:nrow(plot_info)) {
  table = plot_info$table[[i]]
  time_df[[plot_info$condition[[i]]]] = map_dbl(time_vec, function(t) {
    t1 = max(table$time[table$time <= t])
    t2 = min(table$time[table$time >= t])
    if (is.infinite(t1)) {
      t1 = t
      c1 = 0
    } else {
      c1 = table$count[table$time == t1]
    }
    if (is.infinite(t2)) {
      t2 = t
      c2 = 0
    } else {
      c2 = table$count[table$time == t2]
    }
    c1 + (c2-c1)*(t-t1)/(t2-t1)
  })
}



#calculating overall average
all_counts = time_df %>% select(-time_vec)
time_df$average_nonzero = rowSums(all_counts)/rowSums(all_counts != 0)
time_df$average = rowMeans(all_counts)
time_df$average_nonzero[is.na(time_df$average_nonzero)] = 0
#time_df$average[is.na(time_df$average)] = 0

doubling_time = 0.4
time_vec = seq(0,length(breaks)*doubling_time, doubling_time)

time_df = align_time_scale(nested_density, time_vec)



#plotting
plot_tb = gather(data = time_df, key = 'condition', value = 'count', -time_vec)
plot_tb$truth = signal_func(plot_tb$time_vec)
ggplot(data = plot_tb) +
  geom_line(aes(x = time_vec, y = count, group = condition, color = condition)) +
  geom_line(aes(x = time_vec, y = truth, color = 'truth'))


plot_tb = gather(data = time_df %>% select(time_vec, average_nonzero, average), key = 'condition', value = 'count', -time_vec)
plot_tb$truth = signal_func_hgRNA(plot_tb$time_vec)
ggplot(data = plot_tb) +
  geom_line(aes(x = time_vec, y = count/sum(count), group = condition, color = condition)) +
  geom_line(aes(x = time_vec, y = truth/sum(truth), color = 'truth'))
