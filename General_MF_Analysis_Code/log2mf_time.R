raw_data$data
plot_tb = raw_data$data %>% select(avg_log2mf, node_time, gr_node) %>% unique()
node_list = c('ICM', 'Epi', 'PS', 'APS', 'Node')
commitment_bias_list = c(0.75, 0.5, 0.5, 0.4, 0.2)
commitment_bias_list_cum = cumprod(commitment_bias_list)
names(commitment_bias_list_cum) = node_list
plot_tb = plot_tb %>% filter(gr_node %in% node_list)
plot_tb_nest = nest(plot_tb, data = -gr_node)
plot_tb_nest$com_bias = commitment_bias_list_cum[plot_tb_nest$gr_node]
plot_tb_nest$data = map2(plot_tb_nest$data, plot_tb_nest$com_bias, function(data, b) {
  data$avg_log2mf = data$avg_log2mf + log2(b)
  data
})
plot_tb_unnest = unnest(plot_tb_nest, cols = data)
before = ggplot(data = plot_tb, aes(y = node_time, x = -avg_log2mf, color = gr_node)) +
  geom_smooth(method = 'lm')
after = ggplot(data = plot_tb_unnest, aes(y = node_time, x = -avg_log2mf, color = gr_node)) +
  geom_smooth(method = 'lm')
before + after
ggplot(data = plot_tb, aes(x = node_time, y = avg_log2mf, color = gr_node)) +
  geom_point() + geom_smooth() +
  facet_wrap(~ gr_node)


raw_data = load_simulated_data('./three_cell_somatic_mut_ss5000_signal4/', 100, mat_name = 'mf_table')
raw_data = annotate_mut_node(raw_data, phy)
plot_tb = raw_data$data %>% select(avg_log2mf, node_time, gr_node) %>% unique()
ggplot(data = plot_tb, aes(y = node_time, x = -avg_log2mf, color = gr_node)) +
  geom_smooth(method = 'lm') + xlab('Negative Average Log2 Mosaic Fraction') +
  ylab('True Mutation Time') + ggtitle('Mosaic Fraction vs. True Mutation Time For Different Cell Types') +
  theme_pubr() + theme(legend.position = 'right') + labs(colour = 'Cell Type')
