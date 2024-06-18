library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
library(randomcoloR)
library(readxl)
setwd('~/Documents/qfm2')
devtools::load_all()

setwd("~/Documents/R/MosaicFractionAnalysis/R")
devtools::load_all()


setwd('~/Documents/R/MF_Signal_Simulation/lineage specific simulation/')

# load('single_cell_tb_full_tree_mr1.rda')
# load('bulk_tb_full_tree_mr1.rda')
load('3_13_edges.rda')
load('3_13_edge_tb_with_pseudotime_and_program.rda')


raw_data <- structure(list(data = mf_table, MFcomputed = T), class = 'raw')
raw_data$data = raw_data$data %>% select(tissue = CellType, sequence, mosaic_fraction,
                                         node_to_time, node_from_time,
                                         to_type, from_type)
  
raw_data$data  = raw_data$data %>% filter(mosaic_fraction > 0) %>% nest(data = -c(sequence)) %>%
  mutate(avg_log2mf = map_dbl(data, function(data) {log2(mean(data$mosaic_fraction))})) %>% unnest(cols = data)
raw_data$data$ID = 'somatic'
phy = readRDS('~/Documents/R/MF_Signal_Simulation/lineage specific simulation/phy_com_scaled.rds')
#raw_data = annotate_mut_node(raw_data, phy)
raw_data$data$gr_node = raw_data$data$to_type
raw_data$data$sample = raw_data$data$tissue
nested_data = nest_raw(raw_data, c('ID', 'sample', 'gr_node'))
breaks <- c(-Inf, seq(-10, -0.5, by = 0.5),1)
bin_sizes = rev(diff(breaks))
bin_sizes[is.infinite(bin_sizes)] = min(bin_sizes)
nested_data = append_density(nested_data, breaks)

nested_data = normalize_density_list(nested_data, apply_cum_adjustment = FALSE)
nested_data$nested$condition = nested_data$nested$gr_node
nested_density_raw = average_density_by(nested_data, c('condition'), raw = TRUE)
nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)
nested_density = convert_mf_to_time(nested_density)



for (n in edges$in_node) {
  if (n == 'Root') {
    edges$start_time[edges$in_node == n] = 0
  } else {
    edges$start_time[edges$in_node == n] = edges$start_time[edges$out_node == n] +
      edges$edge_length[edges$out_node == n]
  }
}


#make each edge program combination a table with time values added,
#unnest one at a time, and then plot with ggplot
edges$time = map2(edges$start_time,
                  edges$edge_length,
                  function(s, l) {
                    seq(s + 1/200, s+l, 1/200)
                  })

cell_state_col = read_excel('./gas_states_no_conv_mod.xlsx') %>%
  select(name, color)
color_values = cell_state_col$color
names(color_values) = cell_state_col$name

# color_values = c(randomColor(count = length(edges$out_node)))
# names(color_values) = edges$out_node

plot_tb = unnest(edges, cols = time) %>% select(in_node, out_node, time)
for (i in 3:3) {
  name = paste0('program', i)
  plot_tb = plot_tb %>% mutate(edges %>% unnest(cols = name) %>% select(name))
}
plot_tb = plot_tb %>% gather(key = 'program', value = 'value', -c(in_node, out_node, time))

median_program = median(unlist(edges %>% select(program3) %>% unnest(cols = program3)))
base_mean = log(10^(1)) - median_program
base = rnorm(1,base_mean,1)
base = 2.11

plot_tb$value = exp(base + plot_tb$value)

program3 = ggplot(data = plot_tb,
                  aes(x = time, y = value, color = out_node)) +
  scale_colour_manual(values = color_values) +
  geom_line() + facet_wrap(~ program)

plot_tb = nested_density$nested_density %>% select(condition, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    loess = loess(formula = log(count+0.0001) ~ time,
                  data = density$table,
                  span = 1/4)
    table = tibble(time = seq(0, 10, by = 0.01),
                   count = exp(predict(loess, newdata = tibble(time = seq(0, 10, by = 0.01))))-0.0001)
    table
  })) %>%
  unnest(cols = table)
bulk_tb = plot_tb

MF_recovered = ggplot(data = plot_tb,
       aes(x = time, y = count, color = condition)) +
  scale_colour_manual(values = color_values) +
  geom_line() + ylab('Mutation Rate (Mutations/Cell Cycle/Cell') + xlab('Time (Days)') +
  theme_pubr() + theme(axis.text = element_text(size = 10),
                       axis.title = element_text(size = 12),
                       plot.title = element_text(size = 9),
                       legend.position = 'right') +
  labs('colour' = 'Cell Type')

MF_recovered

png(filename = '~/Documents/R/MF_Signal_Simulation/pics/mf_recovered.png', width = 4, height = 4, units = "in", res = 1000)
print(MF_recovered)
dev.off()

programs+MF_recovered

png(filename = '~/Documents/R/MF_Signal_Simulation/pics/side_by_side.png', width = 10, height = 4, units = "in", res = 1000)
print(programs+MF_recovered)
dev.off()

program3 / MF_recovered



plot_tb = nested_density$nested_density %>% select(condition, avg_density) %>%
  mutate(table = map(avg_density, function(density) {
    density$table
  })) %>%
  unnest(cols = table)

ggplot(data = plot_tb %>% filter(condition == 'A/C')) +
  scale_colour_manual(values = color_values) +
  geom_line(aes(x = -bin, y = time, color = condition)) +
  geom_point(data = raw_data$data %>% select(node_to_time, node_from_time, mosaic_fraction, gr_node) %>%
               filter(gr_node == 'A/C'),
             aes(x = -log2(mosaic_fraction), y = node_to_time/2 + node_from_time/2))

