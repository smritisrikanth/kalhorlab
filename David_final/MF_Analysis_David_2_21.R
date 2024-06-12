setwd('~/Documents/R/David/')
raw_data = load_samples('~/Documents/R/David/5-pair_filtering/')

parent = load_parent('./PB21-founder_filteredpairs.txt')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)

raw_data = filter_unmutated(raw_data)

sample_metadata = as_tibble(stringr::str_match(unique(raw_data$data$sample), "(\\d+)-(.*?)_filteredpairs\\.txt"))
colnames(sample_metadata) = c('sample', 'mouse', 'tissue')
raw_data = append_metadata(raw_data, sample_metadata)
colnames(mut_rate_data)[3] = 'ID'
raw_data = append_metadata(raw_data, mut_rate_data, cols = c('Class'), join_by = 'ID')
#raw_data$data = unique(raw_data$data)
animal_age_list = c('5728' = 1.7, '5729' = 1.7, '6554' = 0.58, '6555' = 0.58, '6351' = 1.36,
                    '6352' = 1.36, '6559' = 1.1, '65576558' = 1.1, '5662' = 2.4, '65606561' = 1.1)
raw_data$data$age = map_dbl(raw_data$data$mouse, function(mouse) {
  if (mouse %in% names(animal_age_list)) {
    animal_age_list[[mouse]]
  } else {
    return(0.5)
  }
})

mice = unique(raw_data$data$mouse)
mice = mice[mice != 6553]

raw_data$data$tissue_group = substr(raw_data$data$tissue, 1, 1)
dist_mat_list = map(mice, function(m) {
  message(m)
  mat = dist_mat(data = raw_data$data %>% filter(mouse == m) %>%
                   mutate(sample = tissue_group), dist = 'jaccard')
  mat
})
names_6 = rownames(as.matrix(dist_mat_list[[14]]))
names_7 = rownames(as.matrix(dist_mat_list[[2]]))
names_5 = rownames(as.matrix(dist_mat_list[[10]]))
dist_mat_list_to_avg = map(dist_mat_list, function(mat) {
  orig = mat
  mat = as.matrix(mat)
  if (length(intersect(names(orig), names_7)) == 7) {
    ordered_mat = mat[names_7, names_7]
    return(ordered_mat)
  } else {
    return(NULL)
  }
})
null_dist_mats = map_lgl(dist_mat_list_to_avg, function(mat) is.null(mat))
.dist_mat_list = dist_mat_list_to_avg[!null_dist_mats]
cum_mat = Reduce('+', .dist_mat_list)


avg_mat = cum_mat/length(.dist_mat_list)
avg_phy = nj(as.dist(avg_mat))
avg_phy = name_nodes(avg_phy)
avg_phy = phytools::midpoint.root(avg_phy)

raw_data$data$tissue = substr(raw_data$data$tissue, 1, 1)
raw_data = annotate_mut_node(raw_data, avg_phy)

nested_data = nest_raw(raw_data, c('ID', 'sample', 'mouse', 'tissue', 'age', 'gr_node'))
breaks <- c(-Inf, seq(-10, 0, by = 0.5))
bin_sizes = rev(diff(breaks))
bin_sizes[is.infinite(bin_sizes)] = min(bin_sizes)
nested_data = append_density(nested_data, breaks)
colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
raw_data_filtered = filter_recurring_spacers(raw_data, indelphi_table, 0.0001, filter_missing = FALSE)
nested_data_filtered = nest_raw(raw_data_filtered, c('ID', 'sample', 'mouse', 'tissue', 'age', 'gr_node'))
nested_data_filtered = append_density(nested_data_filtered, breaks)
nested_data = adjust_density(nested_data_filtered, nested_data, c('ID', 'sample', 'mouse', 'tissue', 'age', 'gr_node'))

nested_data = normalize_density_list(nested_data, apply_cum_adjustment = T)
nested_density = average_density_by(nested_data, c('ID', 'mouse', 'gr_node'), normalized = T)

mice = unique(nested_density$nested_density$mouse)
IDs = unique(nested_density$nested_density$ID)

pdf('~/Documents/R/David/MF_plots.pdf')
for (m in mice) {
  message(m)
  for (i in IDs) {
    plot_tb = nested_density$nested_density %>% select(mouse, ID, gr_node, avg_density) %>%
      filter(mouse == m) %>% filter(ID == i) %>%
      mutate(table = map(avg_density, function(density) {
        density$table
      })) %>%
      unnest(cols = table)
    plot_tb$group = paste0(plot_tb$mouse, plot_tb$ID, plot_tb$gr_node)
    if(nrow(plot_tb) == 0) {
      next
    }
    plot = ggplot(data = plot_tb) +
      geom_line(aes(x = -bin, y = count, group = group, color = gr_node)) +
      facet_wrap(~ mouse + ID)
    print(plot)
  }
}
dev.off()
