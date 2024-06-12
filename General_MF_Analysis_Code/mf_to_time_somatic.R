mult_vec = seq(0.1,1,0.01)

doubling_time = 0.65
time_vec = seq(0,length(breaks)*doubling_time, doubling_time)
time_df = align_time_scale(nested_density, time_vec)

time_df_norm = time_df %>% select(-time_vec) %>%
  mutate(truth = signal_func(time_vec)) %>%
  mutate(across(everything(), ~ ./sum(.)/doubling_time)) %>% mutate(time_vec = time_vec)
lm_error = sqrt(AUC(x = time_df_norm$time_vec,y = (time_df_norm$average - time_df_norm$truth)^2))
#lm_intercept = nested_density$nested_density$lm[[1]]$coefficients[1]

log2mf_df = align_counts(nested_density)

tb = tibble(multiplier = mult_vec, error = rep(0,length(mult_vec)))
i = 1
for (m in mult_vec) {
  log2mf_df$time = -log2mf_df$log2mf*m
  log2mf_df_norm = log2mf_df %>% select(-c(log2mf, time)) %>%
    mutate(truth = signal_func(log2mf_df$time)) %>%
    mutate(across(everything(), ~ ./sum(.)/m)) %>%
    mutate(time = log2mf_df$time)
  tb$error[[i]] = sqrt(AUC(x = log2mf_df_norm$time,y = (log2mf_df_norm$average - log2mf_df_norm$truth)^2))
  i = i+1
}

# lm_slope = -mean(map_dbl(nested_density$nested_density$lm, function(lm) {
#   lm$coefficients[2]
# }))

lm_slope = -nested_density$nested_density$lm[[1]]$coefficients[1]

vline_tb = tibble(label = c('Linear Model Estimate:', 'Optimal Linear Multiplier:', 'True Doubling Time:'),
                  value = c(lm_slope,
                            tb$multiplier[which.min(tb$error)], 0.65))
hline_tb = tibble(label = c('Linear Model Error:', 'Optimal Multiplier Error:'), value = c(lm_error, min(tb$error)))

ggplot(data = tb, aes(x = multiplier, y = error)) +
  geom_line() +
  geom_vline(aes(xintercept = value, color = paste0(label, ' ', value)), data = vline_tb) +
  geom_hline(aes(yintercept = value, color = paste0(label, ' ', value)), data = hline_tb) +
  xlab('Linear Multiplier For Log2MF to Time Conversion') +
  ylab('Error') +
  ggtitle('Error Between True and Recovered Signals Using Linear Time Conversion') +
  theme_pubr() + theme(legend.position = 'right') + labs(colour = 'Line Labels')
  
