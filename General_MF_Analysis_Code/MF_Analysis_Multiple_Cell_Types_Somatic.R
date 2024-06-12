library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('~/Documents/qfm2')
devtools::load_all()

# mut_p = readRDS("./metadata//mut_p_marc1.rds")
# mut_p$mut_rate_func = map(1:length(mut_p$mut_rate), function(i) {
#   signal_func_hgRNA
# })
# mut_p$mut_rate = NULL

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

raw_data = load_simulated_data('./mutation simulation output/three_cell_somatic_mut_ss5000_015_2_signal4/', 1, mat_name = 'mf_table')

#load('./mouse_gas_cg_ss2000_signal4/mf_table_1.rda')


raw_data = annotate_mut_node(raw_data, phy)
nested_data = nest_raw(raw_data, c('ID', 'sample', 'gr_node'))
breaks <- c(-Inf, seq(-13.5, -0.5, by = 1), 1)
bin_sizes = rev(diff(breaks))
bin_sizes[is.infinite(bin_sizes)] = min(bin_sizes)
nested_data = append_density(nested_data, breaks)

#correct cum_mf
node_order_list = list('TypeS'=1, 'TypeT'=2, 'TypeC'=3, 'TypeA'=4, 'TypeB'=5)
nested_data = increment_cum_mf(nested_data, node_order_list)

#cumulative adjustment, 1/MF normalization, and averaging across node type
nested_data = normalize_density_list(nested_data, apply_cum_adjustment = FALSE)
nested_data$nested$condition = nested_data$nested$gr_node
nested_density_raw = average_density_by(nested_data, c('condition'), raw = TRUE)
nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)

doubling_time = c(0.6, 0.55, 0.55, 0.55, rep(0.425, 6), 0.55, rep(0.35, 10), rep(0.55, 4))
names(doubling_time) = names(node_order_list)
nested_density$nested_density$doubling_time = map_dbl(nested_density$nested_density$condition,
                                                      function(condition) {
  doubling_time[[condition]]
})
nested_density = convert_mf_to_time(nested_density)
nested_density_raw = convert_mf_to_time(nested_density_raw)


doubling_time = 0.65
time_vec = seq(0,length(breaks)*doubling_time, doubling_time)
time_df = align_time_scale(nested_density, time_vec)
time_df_norm = time_df %>% select(-time_vec) %>%
  mutate(truth = signal_func(time_vec)) %>%
  mutate(across(everything(), ~ ./sum(.)/(doubling_time))) %>%
  mutate(time_vec = time_vec)

error = sqrt(sum((time_df_norm$truth - time_df_norm$average)^2))
plot_tb = gather(data = time_df_norm %>% select(-c(average, average_nonzero, truth)),
                 key = 'condition', value = 'count', -c(time_vec))
plot_tb_avg = gather(data = time_df_norm %>% select(c(time_vec, average, average_nonzero, truth)),
                key = 'condition', value = 'count', -c(time_vec))
plot_tb_avg$condition[plot_tb_avg$condition == 'truth'] = 'Input'
plot_nodes = ggplot(data = plot_tb) +
  geom_line(aes(x = time_vec, y = count, group = condition, color = condition)) +
  ylab('Mosaic Fraction Density') +
  xlab('Time (days)') + theme_pubr() + 
  theme(legend.position = 'right') + labs(colour = 'Cell Type') +
  ggtitle('Mosaic Fraction Density Stratified By Cell Type')
plot_avg = ggplot(data = plot_tb_avg) +
  geom_line(aes(x = time_vec, y = count, group = condition, color = condition)) +
  ylab('Mosaic Fraction Density') +
  xlab('Time (days)') + theme_pubr() + 
  theme(legend.position = 'right') + labs(colour = 'Signal') +
  scale_colour_manual(values = c('Input' = 'black',
                                 'average' = 'green',
                                 'average_nonzero' = 'red')) +
  ggtitle('Average Mosaic Fraction Density')

plot_nodes / plot_avg

log2mf_df = align_counts(nested_density)
log2mf_df$log2mf = log2mf_df$log2mf + log2(3/4)
log2mf_df$time = -log2mf_df$log2mf*0.6
log2mf_df_norm = log2mf_df %>% select(-c(log2mf,time)) %>%
  mutate(truth = signal_func(log2mf_df$time)) %>%
  mutate(across(everything(), ~ ./sum(.))) %>%
  mutate(time = log2mf_df$time) %>%
  mutate(log2mf = log2mf_df$log2mf)

log2mf_df_norm$average = rowMeans(log2mf_df_norm %>% select(ICM, Epi, PS, APS, Node, DE))
plot_tb = gather(data = log2mf_df_norm %>% select(time, truth, ICM, Epi, PS, APS, Node, DE, average),
                 key = 'condition', value = 'count', -c(time))
error = error = sqrt(sum((log2mf_df_norm$truth - log2mf_df_norm$average)^2))
ggplot(data = plot_tb) +
  geom_line(aes(x = time, y = count, group = condition, color = condition)) +
  ylab('Mosaic Fraction Density') +
  xlab('Time (days)') +
  ggtitle(paste0('Error: ', error))


log2mf_df_55 = log2mf_df
log2mf_df_425 = log2mf_df
log2mf_df_35 = log2mf_df
log2mf_df_55$time = -log2mf_df_55$log2mf*0.55
log2mf_df_425$time = -log2mf_df_425$log2mf*0.425
log2mf_df_35$time = -log2mf_df_35$log2mf*0.35
log2mf_df_norm_35 = log2mf_df_35 %>% select(-c(log2mf, time)) %>%
  mutate(truth = signal_func_hgRNA(log2mf_df_35$time)) %>%
  mutate(across(everything(), ~ ./sum(.))) %>%
  mutate(log2mf = log2mf_df_35$log2mf) %>%
  mutate(time = log2mf_df_35$time)

log2mf_df_comb = cbind(log2mf_df_55 %>% select(time_55 = time, ICM, Epi, PS),
                       log2mf_df_425 %>% select(time_425 = time, APS),
                       log2mf_df_35 %>% select(time_35 = time, DE, Node))

log2mf_df_comb_norm = log2mf_df_comb %>% select(-c(time_55, time_425, time_35)) %>%
  mutate(average = rowMeans(log2mf_df_comb %>% select(-c(time_55, time_425, time_35)))) %>%
  mutate(truth = signal_func_hgRNA(log2mf_df_comb$time_55)) %>%
  mutate(across(everything(), ~ ./sum(.))) %>%
  mutate(time_55 = log2mf_df_comb$time_55,
         time_425 = log2mf_df_comb$time_425,
         time_35 = log2mf_df_comb$time_35)
plot_tb_55 = gather(data = log2mf_df_comb_norm %>% select(time_55, ICM, Epi, PS),
                    key = 'condition', value = 'count', -c(time_55))
plot_tb_425 = gather(data = log2mf_df_comb_norm %>% select(time_425, APS),
                    key = 'condition', value = 'count', -c(time_425))
plot_tb_35 = gather(data = log2mf_df_comb_norm %>% select(time_35, DE, Node),
                    key = 'condition', value = 'count', -c(time_35))
ggplot(data = plot_tb_55) +
  geom_line(aes(x = time_55, y = count, group = condition, color = condition)) +
  geom_line(data = plot_tb_425, aes(x = time_425, y = count, group = condition, color = condition)) +
  geom_line(data = plot_tb_35, aes(x = time_35, y = count, group = condition, color = condition)) +
  ylab('Mosaic Fraction Density') +
  xlab('Time (days)')

sqrt(sum((time_df_norm$truth - time_df_norm$average_nonzero)^2))

plot_tb = nested_density_raw$nested_density %>% select(condition, avg_density)
plot_tb$table = map(plot_tb$avg_density, function(density) {
  tibble(time_bin = density$table$bin, count = density$table$count)
})
plot_tb = select(plot_tb, condition, table) %>% unnest(cols = table)
ggplot(data = plot_tb) +
  geom_line(aes(x = time_bin, y = count, group = condition, color = condition))
