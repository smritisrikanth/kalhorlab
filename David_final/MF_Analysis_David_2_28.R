#functions
dist_mat <- function(data,dist=c('jaccard','manhattan')){
  if(dist == 'jaccard'){ # calculate the jaccard index distance 
    ji_df = data.frame(expand.grid(file1=unique(data$sample),file2=unique(data$sample),d=NA))
    if (nrow(ji_df) < 3) {
      return(NULL)
    }
    for(i in 1:nrow(ji_df)){
      if(ji_df$file1[i] != ji_df$file2[i]){
        file1 = data %>% filter(sample == ji_df$file1[i]) %>% 
          transmute(ID = paste0(ID,'_',sequence)) %>% pull(ID) %>% unique(.)
        file2 = data %>% filter(sample == ji_df$file2[i]) %>% 
          transmute(ID = paste0(ID,'_',sequence)) %>% pull(ID) %>% unique(.)
        s = length(intersect(file1,file2)) / length(union(file1,file2))
        ji_df$d[i]=1-s
      }else{
        ji_df$d[i]=0
      }
    }
    D = ji_df %>% pivot_wider(names_from = file1, values_from = d) %>%
      column_to_rownames(var='file2')
    return(D)
  } else{ # add manhattan distance here 
    log_binary_matrix = normalized_matrix(data, value_choice = 'logMF')
    return(dist(log_binary_matrix,method = "manhattan"))
    
  }
  phy = name_nodes(phy)
  phy = phytools::midpoint.root(phy)
  return(phy)
}

normalized_matrix <- function(data, value_choice = c("logMF", "binary")){
  data_2 = data %>% transmute(ID = paste0(ID,'_',sequence),reads=mosaic_fraction,sample=sample) %>% 
    pivot_wider(names_from = ID,values_from = reads) %>%
    column_to_rownames(.,var='sample') %>%
    mutate_all(~ifelse(is.na(.), 0, .))
  normalized_matrix <- t(apply(data_2, 1, function(x) x/sum(x)))
  if(value_choice == 'logMF'){
    log_matrix <- log2(normalized_matrix + 1)
  }else{
    log_matrix <- 1*(log(normalized_matrix + 1) > 0)
  }
  return(log_matrix)
}

setwd('~/Documents/R/David/')

setwd('~/Documents/R/Jasmine/')
indelphi_table <- readr::read_tsv('MARC1_indelphi_result.txt', col_names = T)
mut_rate_data <- read_csv('SupplementaryTable3 - TableS3-hgRNA Activities.csv')

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
raw_data$data$tissue_group[raw_data$data$tissue_group == 'C'] = 
  substr(raw_data$data$tissue[raw_data$data$tissue_group == 'C'], 1, 2)
raw_data$data$tissue_group[raw_data$data$tissue_group == 'CL' |
                             raw_data$data$tissue_group == 'CR'] = 'C'

colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
raw_data = filter_recurring_spacers(raw_data, indelphi_table, 0.0001, filter_missing = FALSE)

dist_mat_list = map(mice, function(m) {
  message(m)
  mat = dist_mat(data = raw_data$data %>% filter(mouse == m) %>%
                   mutate(sample = tissue_group), dist = 'jaccard')
  mat
})

map(dist_mat_list, function(mat) {
  names(mat)
})
names_7 = rownames(as.matrix(dist_mat_list[[6]]))
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


#get diff average tree with mice with max num tissues

raw_data$data$tissue = raw_data$data$tissue_group
raw_data = annotate_mut_node(raw_data, avg_phy)

#print trees for every mouse
pdf('~/Documents/R/David/2_28_mouse_trees_old.pdf')
for (m in mice) {
  print(m)
  mat = dist_mat(data = raw_data$data %>% filter(mouse == m) %>%
                   mutate(sample = tissue_group), dist = 'jaccard')
  phy = nj(as.dist(mat))
  phy = name_nodes(phy)
  phy = phytools::midpoint.root(phy)
  
  plot(phy)
  title(m)
}
dev.off()

#old data
raw_data = load_samples('~/Documents/R/David/')

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
raw_data$data$tissue_group[raw_data$data$tissue_group == 'C'] = 
  substr(raw_data$data$tissue[raw_data$data$tissue_group == 'C'], 1, 2)
raw_data$data$tissue_group[raw_data$data$tissue_group == 'CL' |
                             raw_data$data$tissue_group == 'CR'] = 'C'

dist_mat_list = map(mice, function(m) {
  message(m)
  mat = dist_mat(data = raw_data$data %>% filter(mouse == m) %>%
                   mutate(sample = tissue_group), dist = 'jaccard')
  mat
})

map_dbl(dist_mat_list, function(mat) {
  length(names(mat))
})
names = rownames(as.matrix(dist_mat_list[[6]]))
dist_mat_list_to_avg = map(dist_mat_list, function(mat) {
  orig = mat
  mat = as.matrix(mat)
  if (length(intersect(names(orig), names)) == 7) {
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


pdf('~/Documents/R/David/mouse_trees_old.pdf')
for (m in mice) {
  print(m)
  mat = dist_mat(data = raw_data$data %>% filter(mouse == m) %>%
                   mutate(sample = tissue_group), dist = 'jaccard')
  phy = nj(as.dist(mat))
  phy = name_nodes(phy)
  phy = phytools::midpoint.root(phy)
  
  plot(phy)
  title(m)
}
dev.off()








