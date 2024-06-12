colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
raw_data = filter_recurring_spacers(raw_data, indelphi_table, 0.001, filter_missing = FALSE)


raw_data$data$tissue_group = substr(raw_data$data$tissue, 1, 1)
raw_data$data$tissue_group[raw_data$data$tissue_group == 'C'] = 
  substr(raw_data$data$tissue[raw_data$data$tissue_group == 'C'], 1, 2)
raw_data$data$tissue_group[raw_data$data$tissue_group == 'CL' |
                             raw_data$data$tissue_group == 'CR'] = 'C'

tissue_list = unique(raw_data$data$tissue)
tissue_group_list = unique(raw_data$data$tissue_group)

summarise = raw_data$data %>% mutate(sample = tissue_group) %>%
  group_by(sample, mouse, ID, sequence) %>%
  summarise(avg_MF = mean(mosaic_fraction), reads = mean(reads))

mouse_list = unique(summarise$mouse)

tree_list = map(mouse_list, function(m) {
  build_tree(data = summarise %>% ungroup()
             %>% filter(mouse == m), dist = 'jaccard')
})
pdf(file = '~/Documents/R/David/mouse_trees.pdf', onefile = T)
for (i in (1:length(tree_list))) {
  plot(tree_list[[i]])
  title(mouse_list[[i]])
}
dev.off()
consensus_tree = ape::consensus(tree_list, p = 0.5, rooted = TRUE)

dist_mat_list = map(mouse_list, function(m) {
  dist_mat(data = summarise %>% ungroup() %>% mutate(sample = tissue_group)
           %>% filter(mouse == m), dist = 'jaccard')
})
cum_mat = matrix(0,nrow = nrow(dist_mat_list[[1]]), ncol = ncol(dist_mat_list[[1]]))
for (mat in dist_mat_list) {
  cum_mat = cum_mat + mat
}
dev.off()

avg_mat = cum_mat/length(dist_mat_list)
average_tree = phytools::midpoint.root(nj(as.dist(avg_mat)))



#manhattan
tree_list = map(mouse_list, function(m) {
  build_tree(data = summarise %>% ungroup() %>% filter(mouse == m) %>% select(-avg_MF),
             dist = 'manhattan')
})
pdf(file = '~/Documents/R/David/mouse_trees.pdf', onefile = T)
for (tree in tree_list) {
  plot(tree)
}
dev.off()
consensus_tree = ape::consensus(tree_list, p = 0.5, rooted = TRUE)
plot(consensus_tree)

dist_mat_list = map(mouse_list, function(m) {
  as.matrix(dist_mat(data = summarise %>% ungroup() %>% filter(mouse == m), dist = 'manhattan'))
})
cum_mat = matrix(0,nrow = nrow(dist_mat_list[[1]]), ncol = ncol(dist_mat_list[[1]]))
for (mat in dist_mat_list) {
  cum_mat = cum_mat + as.matrix(mat)
}

avg_mat = cum_mat/length(dist_mat_list)
average_tree = phytools::midpoint.root(nj(as.dist(avg_mat)))
average_tree = name_nodes(average_tree)
