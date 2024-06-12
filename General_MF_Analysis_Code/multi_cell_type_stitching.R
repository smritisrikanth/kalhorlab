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

signal_func = function(x_vec) {
  val = 1e-9 * ((1/6)* pmin(pmax(0, x_vec-1), 4) + (-1/6) * pmax(0, x_vec-5))
  out_vec = pmax(val, 0) #+ 2*10^(-9)/3
  #out_vec[x_vec < 1] = 0
  out_vec
}

raw_data = load_simulated_data('./three_cell_somatic_mut_ss5000_015_2_signal4/', 1, mat_name = 'mf_table')

raw_data = annotate_mut_node(raw_data, phy)
nested_data = nest_raw(raw_data, c('ID', 'sample', 'gr_node'))
breaks <- c(-Inf, seq(-13.5, -0.5, by = 1), 1)
bin_sizes = rev(diff(breaks))
bin_sizes[is.infinite(bin_sizes)] = min(bin_sizes)
nested_data = append_density(nested_data, breaks)

#correct cum_mf
node_order_list = list('TypeS'=1, 'TypeT'=2, 'TypeC'=2, 'TypeA'=3, 'TypeB'=3)
nested_data = increment_cum_mf(nested_data, node_order_list)

#cumulative adjustment, 1/MF normalization, and averaging across node type
nested_data = normalize_density_list(nested_data, apply_cum_adjustment = FALSE)
nested_data$nested$condition = nested_data$nested$gr_node
nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)

log2mf_df = align_counts(nested_density)
plot_tb = gather(data = log2mf_df,
                 key = 'condition', value = 'count', -c(log2mf))
ggplot(data = plot_tb) +
  geom_line(aes(x = -log2mf, y = count, group = condition, color = condition)) +
  ylab('Mosaic Fraction Density') +
  xlab('-log2mf')

single_lineage = c('TypeS', 'TypeT', 'TypeA')
change_points = map_dbl(1:(length(single_lineage)-1), function(i) {
  node = single_lineage[i]
  next_node = single_lineage[i+1]
  list_proportion = map_dbl(3:length(-log2mf_df$log2mf)-2, function(p) {
    area1 = AUC(x = -log2mf_df$log2mf[1:p], y = -log2mf_df[[node]][1:p])
    area2 = AUC(x = -log2mf_df$log2mf, y = -log2mf_df[[node]])
    area3 = AUC(x = -log2mf_df$log2mf[p+1:length(-log2mf_df$log2mf)],
                y = -log2mf_df[[next_node]][p+1:length(-log2mf_df$log2mf)])
    area4 = AUC(x = -log2mf_df$log2mf, y = -log2mf_df[[next_node]])
    area1/area2 + area3/area4
  })
  which.max(list_proportion)
})
names(change_points) = single_lineage[1:length(single_lineage)-1]


log2mf_df_norm = log2mf_df %>% select(-c(log2mf)) %>%
  mutate(across(everything(), ~ ./sum(.))) %>%
  mutate(log2mf = log2mf_df$log2mf)



