library(ComplexHeatmap)

method_tb = as_tibble(expand.grid(distance_method = c('manhattan', 'jaccard'),
                                  reconstruction_method = c('upgma', 'nj', 'ward2')))

method_tb$comp_matrix = map2(method_tb$distance_method,
                             method_tb$reconstruction_method,
                             function(d,r) {
                               data = param_tb[param_tb$distance_method == d &
                                                 param_tb$reconstruction_method == r,] %>%
                                 select(filtering, dropout, percent_complete) %>%
                                 group_by(filtering, dropout) %>%
                                 summarise(avg_percent_complete = mean(percent_complete))
                               
                               comp_matrix = pivot_wider(data,
                                                         names_from = filtering,
                                                         values_from = avg_percent_complete)
                             })
pdf(file = 'filtering_dropout_matrices.pdf', onefile = T)
for (i in 1:nrow(method_tb)) {
  comp_matrix = method_tb$comp_matrix[[i]]
  print(Heatmap(comp_matrix[,-1], cluster_rows = F, cluster_columns = F))
}
Heatmap()
dev.off()

print(Heatmap(method_tb$comp_matrix[[6]][,-1], cluster_rows = F, cluster_columns = F))

njjaccard = param_tb[param_tb$distance_method == 'jaccard' &
                       param_tb$reconstruction_method == 'nj',] %>%
  select(filtering, dropout, percent_complete) %>%
  group_by(filtering, dropout) %>%
  summarise(avg_percent_complete = mean(percent_complete))

comp_matrix = pivot_wider(njjaccard,
                          names_from = filtering,
                          values_from = avg_percent_complete)


Heatmap(comp_matrix[,-1], cluster_rows = F, cluster_columns = F, )
