raw_data = load_simulated_data('./three_cell_somatic_T_095_5_signal4/', 100, mat_name = 'mf_table')
#phy = load('./phy.rda')
raw_data = annotate_mut_node(raw_data, phy)
nested_data = nest_raw(raw_data, c('ID', 'sample', 'gr_node'))
breaks <- c(-Inf, seq(-13.5, -0.5, by = 1), 1)
nested_data = append_density(nested_data, breaks)
node_order_list = list('TypeS'=1, 'TypeT'=2, 'TypeC'=3, 'TypeA'=4, 'TypeB'=5)
nested_data = increment_cum_mf(nested_data, node_order_list)
nested_data = normalize_density_list(nested_data, apply_cum_adjustment = FALSE)
nested_data$nested$condition = nested_data$nested$gr_node
nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)

#bias correction
commitment_bias_list = c('TypeS' = 1, 'TypeT' = 0.95, 'TypeA' = 0.7*0.95, 'TypeB' = 0.3*0.95, 'TypeC' = 0.05)
nested_density$nested_density$commitment_bias = map_dbl(nested_density$nested_density$condition,
                                                        function(node) {
  commitment_bias_list[[node]]
})
nested_density_not_corrected = convert_mf_to_time(nested_density)
nested_density_corrected = convert_mf_to_time_with_bias_correction(nested_density)
nested_density = nested_density_not_corrected

doubling_time = 0.65
time_vec = seq(0,length(breaks)*doubling_time, doubling_time)

time_df = align_time_scale(nested_density, time_vec)
time_df_norm = time_df %>% select(-time_vec) %>%
  mutate(Input = signal_func(time_vec)) %>%
  mutate(across(everything(), ~ ./sum(.)/(doubling_time))) %>%
  mutate(time_vec = time_vec)
error = sqrt(AUC(x = time_vec, y = (time_df_norm$average - time_df_norm$Input)^2))


#bias correction without true time conversion
commitment_bias_list = c('TypeS' = 1, 'TypeT' = 0.95, 'TypeA' = 0.7*0.95, 'TypeB' = 0.3*0.95, 'TypeC' = 0.05)
nested_density$nested_density$commitment_bias = map_dbl(nested_density$nested_density$condition,
                                                        function(node) {
                                                          commitment_bias_list[[node]]
                                                        })
nested_density_corrected = nested_density
nested_density_corrected$nested_density$avg_density = map2(nested_density_corrected$nested_density$avg_density,
                                                 nested_density_corrected$nested_density$commitment_bias,
                                                 function(density, b) {
                                                   density$table$bin = density$table$bin + log2(b)
                                                   density$table$time = -density$table$bin
                                                   density
                                                 })
doubling_time = 0.65
log2mf_df = align_counts(nested_density) %>% mutate(time = -log2mf*doubling_time)
log2mf_df_norm = log2mf_df %>% select(-c(time, log2mf)) %>%
  mutate(Input = signal_func(log2mf_df$time)) %>%
  mutate(across(everything(), ~ ./sum(.)/(doubling_time))) %>%
  mutate(time_vec = log2mf_df$time)

time_df_norm = log2mf_df_norm
error = sqrt(AUC(x = time_df_norm$time_vec, y = (time_df_norm$average - time_df_norm$Input)^2))


plot_tb = gather(data = log2mf_df %>% select(-c(average, average_nonzero)),
                 key = 'condition', value = 'count', -c(time))
plot_tb_avg = gather(data = time_df_norm %>% select(c(time_vec, average, average_nonzero, Input)),
                     key = 'condition', value = 'count', -c(time_vec))
plot_nodes = ggplot(data = plot_tb) +
  geom_line(aes(x = time_vec, y = count, group = condition, color = condition)) +
  ylab('Mosaic Fraction Density') +
  xlab('Time (days)') + theme_pubr() +
  theme(legend.position = 'right') + labs(colour = 'Cell Type') +
  ggtitle('Uncorrected')
plot_avg = ggplot(data = plot_tb_avg) +
  geom_line(aes(x = time_vec, y = count, group = condition, color = condition)) +
  ylab('Mosaic Fraction Density') +
  xlab('Time (days)') + theme_pubr() +
  theme(legend.position = 'right') + labs(colour = 'Signal') +
  scale_colour_manual(values = c('Input' = 'black',
                                 'average' = 'green',
                                 'average_nonzero' = 'red')) +
  ggtitle(paste0('Error: ', error))


time_df = align_time_scale(nested_density_corrected, time_vec)
time_df_norm = time_df %>% select(-time_vec) %>%
  mutate(Input = signal_func(time_vec)) %>%
  mutate(across(everything(), ~ ./sum(.)/(doubling_time))) %>%
  mutate(time_vec = time_vec)
error = sqrt(AUC(x = time_vec, y = (time_df_norm$average - time_df_norm$Input)^2))

doubling_time = 0.65
time_df = align_time_scale(nested_density_corrected,
                           nested_density_corrected$nested_density$avg_density[[1]]$table$time) %>%
  mutate(time_vec = time_vec*doubling_time)
time_df_norm = time_df %>% select(-time_vec) %>%
  mutate(Input = signal_func(time_df$time_vec)) %>%
  mutate(across(everything(), ~ ./sum(.)/(doubling_time))) %>%
  mutate(time_vec = time_df$time_vec)
error = sqrt(AUC(x = time_df_norm$time_vec, y = (time_df_norm$average - time_df_norm$Input)^2))

plot_tb = gather(data = time_df %>% select(-c(average, average_nonzero)),
                 key = 'condition', value = 'count', -c(time_vec))
plot_tb_avg = gather(data = time_df_norm %>% select(c(time_vec, average, average_nonzero, Input)),
                     key = 'condition', value = 'count', -c(time_vec))
plot_nodes_corrected = ggplot(data = plot_tb) +
  geom_line(aes(x = time_vec, y = count, group = condition, color = condition)) +
  ylab('Mosaic Fraction Density') +
  xlab('Time (days)') + theme_pubr() +
  theme(legend.position = 'right') + labs(colour = 'Cell Type') +
  ggtitle('Corrected')
plot_avg_corrected = ggplot(data = plot_tb_avg) +
  geom_line(aes(x = time_vec, y = count, group = condition, color = condition)) +
  ylab('Mosaic Fraction Density') +
  xlab('Time (days)') + theme_pubr() +
  theme(legend.position = 'right') + labs(colour = 'Signal') +
  scale_colour_manual(values = c('Input' = 'black',
                                 'average' = 'green',
                                 'average_nonzero' = 'red')) +
  ggtitle(paste0('Error: ', error))

(plot_nodes + plot_avg) / (plot_nodes_corrected + plot_avg_corrected) + plot_layout(guides = 'collect')
