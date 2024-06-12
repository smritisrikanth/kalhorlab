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

#load('phy.rda')
phy = readRDS('gast_phylo.rds')

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

raw_data = load_simulated_data('./one_cell_hgRNA_mut_signal_flat/', 100)
#raw_data$data = filter(raw_data$data, mosaic_fraction > 0)
raw_data = annotate_mut_node(raw_data, phy)
nested_data_original = nest_raw(raw_data_original, c('ID', 'sample'))
#nested_data_tissue = nest_raw(raw_data, c('ID', 'sample', 'tissue'))
i = -9
raw_data$data = raw_data_original$data[raw_data_original$data$probability < 2^i,]
nested_data_filtered = nest_raw(raw_data, c('ID', 'sample'))
#nested_data_filtered_tissue = nest_raw(raw_data, c('ID', 'sample', 'tissue'))
breaks <- c(-Inf, seq(-13.5, -0.5, by = 1), 1)
nested_data_original = append_density(nested_data_original, breaks)
nested_data_filtered = append_density(nested_data_filtered, breaks)
nested_data_filtered = adjust_density(nested_data_filtered, nested_data_original, c('ID', 'sample'))
nested_data = nested_data_filtered
# nested_data_tissue = append_density(nested_data_tissue, breaks)
# nested_data_filtered_tissue = append_density(nested_data_filtered_tissue, breaks)
# nested_data_filtered_tissue = adjust_density(nested_data_filtered_tissue, nested_data_tissue, c('ID', 'sample', 'tissue'))
# nested_data_tissue = nested_data_filtered_tissue

#correct cum_mf
node_order_list = list('TypeS'=1, 'TypeT'=2, 'TypeC'=3, 'TypeA'=4, 'TypeB'=5)
node_order_list = c(1:25)
names(node_order_list) = c(phy$node.label, phy$tip.label)
nested_data = increment_cum_mf(nested_data, node_order_list)
#snested_data_tissue = increment_cum_mf(nested_data_tissue, node_order_list, nested_by_node = FALSE, by_var = c('ID', 'sample', 'tissue'))

#cumulative adjustment, 1/MF normalization, and averaging across node type
nested_data_unadjusted = normalize_density_list(nested_data, apply_cum_adjustment = FALSE)
nested_data_adjusted = normalize_density_list(nested_data, apply_cum_adjustment = TRUE)
nested_data_adjusted_tissue = normalize_density_list(nested_data_tissue, apply_cum_adjustment = TRUE)
#nested_data_unadjusted$nested$condition = nested_data_unadjusted$nested$gr_node
nested_data_adjusted$nested$condition = nested_data_adjusted$nested$gr_node
nested_data_adjusted$nested$condition = 'all'
nested_data_adjusted_tissue$nested$condition = nested_data_adjusted_tissue$nested$tissue
#nested_density_raw = average_density_by(nested_data_unadjusted, c('condition'), raw = TRUE)
#nested_density_unadjusted = average_density_by(nested_data_unadjusted, c('condition'), normalized = TRUE)
nested_density_adjusted = average_density_by(nested_data_adjusted, c('condition'), normalized = TRUE)
nested_density_adjusted_tissue = average_density_by(nested_data_adjusted_tissue, c('condition'), normalized = TRUE)

#time conversion
#nested_density_raw = convert_mf_to_time(nested_density_raw)
#nested_density_unadjusted = convert_mf_to_time(nested_density_unadjusted)
nested_density_adjusted = convert_mf_to_time(nested_density_adjusted)
#nested_density_adjusted_tissue = convert_mf_to_time(nested_density_adjusted_tissue)

plot_tb = nested_density_adjusted$nested_density$avg_density[[1]]$table
ggplot(data = plot_tb, aes(x = -bin, y = count)) + geom_line()

#bind matrix with counts for all conditions along time_vec
doubling_time = 0.65
time_vec = seq(0,(length(breaks)-2)*doubling_time, doubling_time)

time_df = align_time_scale(nested_density_adjusted, time_vec)
time_df_norm = time_df %>% select(-time_vec) %>%
  mutate(Input = signal_func_hgRNA(time_vec)) %>%
  mutate(across(everything(), ~ ./sum(.)/doubling_time))
time_df_norm$time_vec = time_vec
time_df_norm$log2mf = nested_density_adjusted$nested_density$avg_density[[1]]$table$bin

time_df = align_time_scale(nested_density_adjusted_tissue, time_vec)
time_df_norm_tissue = time_df %>% select(-time_vec) %>%
  mutate(truth = signal_func_hgRNA(time_vec)) %>%
  mutate(across(everything(), ~ ./sum(.)))
time_df_norm_tissue$time_vec = time_vec
time_df_norm_tissue$log2mf = nested_density_adjusted$nested_density$avg_density[[1]]$table$bin


time_df_norm_comb = cbind(select(time_df_norm, time_vec, log2mf, average, average_nonzero, truth),
                          select(time_df_norm_tissue, TypeA, TypeB, TypeC))

log2mf_df = align_counts(nested_density_adjusted)
# log2mf_df_tissue = align_counts(nested_density_adjusted_tissue)
# log2mf_df_comb = cbind(select(log2mf_df, log2mf, average, average_nonzero),
#                        select(log2mf_df_tissue, TypeA, TypeB, TypeC))
# log2mf_df_comb$time = -log2mf_df_comb$log2mf*0.65
log2mf_df$time = -log2mf_df$log2mf*0.5
log2mf_df_norm = log2mf_df %>% select(-c(log2mf, time)) %>%
  mutate(truth = signal_func_hgRNA(log2mf_df$time)) %>%
  mutate(across(everything(), ~ ./sum(.)))
log2mf_df_norm$log2mf = log2mf_df$log2mf
log2mf_df_norm$time = log2mf_df$time

comp_col = c("TypeA"  = "red",
             "TypeB" = "yellow",
             "TypeC" = "blue",
             "TypeT" = "green",
             "TypeS" = "pink",
             "average" = "black",
             "average_nonzero" = "grey80",
             "truth" = "orange")

comp_col_tissue = c("A"  = "red",
             "B" = "green",
             "C" = "blue",
             "average" = "black",
             "truth" = "orange")

#plotting
plot_tb = gather(data = time_df_norm %>% select(average, average_nonzero, Input, time_vec),
                 key = 'condition', value = 'count', -c(time_vec))
colors = c('Input' = 'black', 'average' = 'green', 'average_nonzero' = 'red')
ggplot(data = plot_tb) +
  geom_line(aes(x = time_vec, y = count, group = condition, color = condition)) +
  scale_color_manual(values = colors) +
  xlab('Time (days)') +
  ylab('Mosaic Fraction Density') +
  theme_pubr() +
  theme(legend.position = "right") + labs(colour = 'Signal Type') +
  ggtitle('Signal Recovered From hgRNA Mutations Within Samples With Multiple Cell Type')


plot_tb = nested_density$nested_density %>% select(condition, avg_density)
plot_tb$table = map(plot_tb$avg_density, function(density) {
  tibble(time_bin = density$table$bin, count = density$table$count)
})
plot_tb = select(plot_tb, condition, table) %>% unnest(cols = table)
ggplot(data = plot_tb) +
  geom_line(aes(x = time_bin, y = count, group = condition, color = condition))


