save(lin_mod_tb, file = 'lin_mod_tb.rda')
load('lin_mod_tb.rda')

lin_mod_tb = raw_data$data %>% select(gr_node, mf_to_time_tb) %>%
  nest(data = -gr_node)

lin_mod_tb$fit = map(lin_mod_tb$data, function(data) {
  log2mf = c()
  node_time = c()

  for (tb in data) {
    log2mf = c(log2mf, tb$log2mf)
    node_time = c(node_time, tb$node_time)
  }
  lm(formula = node_time ~ log2mf)
})
