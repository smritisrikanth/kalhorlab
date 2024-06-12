dist_vec = as.vector(c())
mult_vec = seq(0.1,1,0.01)
cutoff_vec = seq(-15, -5, 1)
breaks <- c(-Inf, seq(-13.5, -0.5, by = 1), 1)

tb = tibble(cutoff = cutoff_vec, error = map(1:length(cutoff_vec), function(x) NULL))

for (i in cutoff_vec) {
  message(i)
  raw_data$data = raw_data_original$data[raw_data_original$data$probability < 2^i,]
  nested_data_filtered = nest_raw(raw_data, c('ID', 'sample', 'gr_node'))
  nested_data_filtered = append_density(nested_data_filtered, breaks)
  nested_data_filtered = adjust_density(nested_data_filtered, nested_data_original, c('ID', 'sample', 'gr_node'))
  nested_data = nested_data_filtered
  nested_data = increment_cum_mf(nested_data, node_order_list)
  nested_data_adjusted = normalize_density_list(nested_data, apply_cum_adjustment = TRUE)
  nested_data_adjusted$nested$condition = nested_data_adjusted$nested$gr_node
  nested_density_adjusted = average_density_by(nested_data_adjusted, c('condition'), normalized = TRUE)

  log2mf_df_comb = align_counts(nested_density_adjusted)

  for (m in mult_vec) {
    log2mf_df_comb$time = -log2mf_df_comb$log2mf*m
    log2mf_df_norm = log2mf_df_comb %>% select(-c(log2mf, time)) %>%
      mutate(truth = signal_func_hgRNA(log2mf_df_comb$time)) %>%
      mutate(across(everything(), ~ ./sum(.)))
    log2mf_df_norm$log2mf = log2mf_df_comb$log2mf
    log2mf_df_norm$time = log2mf_df_comb$time

    # plot_tb = gather(data = log2mf_df_norm %>% select(time, log2mf, truth, average)
    #                  , key = 'condition', value = 'count', -c(log2mf, time))
    # ggplot(data = plot_tb) +
    #   geom_line(aes(x = time, y = count, group = condition, color = condition)) +
    #   scale_color_manual(values = comp_col) +
    #   #xlab('Time (days)') +
    #   ylab('Mosaic Fraction Density') +
    #   theme_pubr() +
    #   theme(legend.position = "right")

    dist = sqrt(sum((log2mf_df_norm$truth - log2mf_df_norm$average_nonzero)^2))
    dist_vec = c(dist_vec, dist)
  }
  tb$error[[i-cutoff_vec[1]+1]] = tibble(multipler = mult_vec, error = dist_vec)
  dist_vec = as.vector(c())
}

for (m in mult_vec) {
  log2mf_df_comb$time = -log2mf_df_comb$log2mf*m
  log2mf_df_norm = log2mf_df_comb %>% select(-c(log2mf, time)) %>%
    mutate(truth = signal_func_hgRNA(log2mf_df_comb$time)) %>%
    mutate(across(everything(), ~ ./sum(.)))
  log2mf_df_norm$log2mf = log2mf_df_comb$log2mf
  log2mf_df_norm$time = log2mf_df_comb$time

  plot_tb = gather(data = log2mf_df_norm %>% select(time, log2mf, truth, average)
                     , key = 'condition', value = 'count', -c(log2mf, time))
  ggplot(data = plot_tb) +
    geom_line(aes(x = time, y = count, group = condition, color = condition)) +
    scale_color_manual(values = comp_col) +
    #xlab('Time (days)') +
    ylab('Mosaic Fraction Density') +
    theme_pubr() +
    theme(legend.position = "right")

  dist = sqrt(sum((log2mf_df_norm$truth - log2mf_df_norm$average_nonzero)^2))
  dist_vec = c(dist_vec, dist)
}

tb$min_error = map_dbl(tb$error, function(data) {
  data$multipler[which.min(data$error)]
})
plot_tb = unnest(tb, cols = error)

argmin_multiplier_error_tb = plot_tb %>% group_by(cutoff) %>%
  summarise(min_error_multiplier = multipler[which.min(error)], min_error = min(error))

ggplot(data = plot_tb, aes(x = multipler, y = error, group = cutoff, color = cutoff)) +
  geom_line() +
  geom_vline(data = argmin_multiplier_error_tb,
             aes(xintercept = min_error_multiplier, color = cutoff)) +
  ggtitle(paste0('min_cutoff: ', argmin_multiplier_error_tb$cutoff[which.min(argmin_multiplier_error_tb$min_error)],
                 '; min_multiplier: ', argmin_multiplier_error_tb$min_error_multiplier[which.min(argmin_multiplier_error_tb$min_error)]))

ggplot(data = argmin_multiplier_error_tb, aes(x = cutoff, y = min_error_multiplier)) +
  geom_line() +
  geom_point()
