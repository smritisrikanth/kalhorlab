library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
library(randomcoloR)
setwd('~/Documents/qfm2')
devtools::load_all()

#functions
# Manhhatan

dist_mat <- function(raw_data, dist_metric) {
  

  raw_data$mf_wide = raw_data$data %>%
    ungroup() %>%
    transmute(sample, ID_sequence = paste0(ID, "_", sequence), mosaic_fraction) %>% 
    group_by(sample, ID_sequence) %>% summarise(mosaic_fraction = mean(mosaic_fraction)) %>%
    spread(key = ID_sequence, value = mosaic_fraction, fill = 0)
  raw_data$mf_mat = as.matrix(raw_data$mf_wide[-1])
  rownames(raw_data$mf_mat) = raw_data$mf_wide[[1]]
  
  if (dist_metric == 'manhattan') {
    dmat = dist(raw_data$mf_mat, method = "manhattan")
    return(dmat)
  } else {
    mf_ind = 1 * (raw_data$mf_mat > 0)
    computePairwiseJaccard <- function(binaryMatrix) {
      numRows <- nrow(binaryMatrix)
      jaccardMatrix <- matrix(1, nrow = numRows, ncol = numRows)
      for (i in 1:(numRows - 1)) {
        for (j in (i + 1):numRows) {
          intersection <- sum(binaryMatrix[i, ] & binaryMatrix[j, ])
          union <- sum(binaryMatrix[i, ] | binaryMatrix[j, ])
          if (union != 0) {
            jaccardMatrix[i, j] <- intersection / union
            jaccardMatrix[j, i] <- jaccardMatrix[i, j]  # Matrix is symmetric
          }
        }
      }
      rownames(jaccardMatrix) = colnames(jaccardMatrix) = rownames(binaryMatrix)
      # Return the Jaccard matrix
      return(jaccardMatrix)
    }
    jaccard_dist = 1 - computePairwiseJaccard(mf_ind)
    #diag(jaccard_dist) = NA
    return(jaccard_dist)
  }
  
}


setwd("~/Documents/R/MosaicFractionAnalysis/R")
devtools::load_all()

setwd('~/Documents/R/Jasmine/')
indelphi_table <- readr::read_tsv('MARC1_indelphi_result.txt', col_names = T)
mut_rate_data <- read_csv('SupplementaryTable3 - TableS3-hgRNA Activities.csv')

raw_data = load_samples('~/Documents/R/David/5-pair_filtering/')
raw_data$data = raw_data$data[raw_data$data$sample != 'PB21-founder_filteredpairs.txt',]

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


#individual mice heatmaps
# m = mice[[2]]
# 
# ID_sample_tb = raw_data$data %>% filter(mouse == m) %>%
#   count(ID, sample) %>% pivot_wider(names_from = ID,values_from = n)
# ID_sample_matrix = as.matrix(ID_sample_tb[,-1])
# colnames(ID_sample_matrix) = names(ID_sample_tb)[-1]
# rownames(ID_sample_matrix) = ID_sample_tb[[1]]
# ID_sample_matrix[is.na(ID_sample_matrix)] = 0
# ID_sample_matrix = log2(ID_sample_matrix+1)
# Heatmap(ID_sample_matrix)
# 
# parent = load_parent('./PB21-founder_filteredpairs.txt')
# raw_data_id_filtered = id_filter(raw_data, parent)
# 
# ID_sample_tb = raw_data_id_filtered$data %>% filter(mouse == m) %>%
#   count(ID, sample) %>% pivot_wider(names_from = ID,values_from = n)
# ID_sample_matrix = as.matrix(ID_sample_tb[,-1])
# colnames(ID_sample_matrix) = names(ID_sample_tb)[-1]
# rownames(ID_sample_matrix) = ID_sample_tb[[1]]
# ID_sample_matrix[is.na(ID_sample_matrix)] = 0
# ID_sample_matrix = log2(ID_sample_matrix+1)
# Heatmap(ID_sample_matrix)



#filter IDs that are most likely contamination
for (m in mice) {
  m_data = raw_data$data %>% filter(mouse == m)
  
  sample_id_ind_wide = select(m_data, ID, sample) %>%
    distinct() %>%
    mutate(value = 0.) %>%
    spread(key = ID, value = value, fill = 1.)
  
  IDs_to_filter = colnames(sample_id_ind_wide[,-1])[colSums(sample_id_ind_wide[,-1]) > 1]
  
  raw_data$data = raw_data$data %>% filter(!(mouse == m & ID %in% IDs_to_filter))
}

#rest of raw_data processing
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)

raw_data = filter_unmutated(raw_data)

raw_data$data$tissue_group = substr(raw_data$data$tissue, 1, 1)
raw_data$data$tissue_group[raw_data$data$tissue_group == 'C'] = 
  substr(raw_data$data$tissue[raw_data$data$tissue_group == 'C'], 1, 2)
raw_data$data$tissue_group[raw_data$data$tissue_group == 'CL' |
                             raw_data$data$tissue_group == 'CR'] = 'C'

colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
raw_data_filtered = filter_recurring_spacers(raw_data, indelphi_table, 0.0001, filter_missing = FALSE)

mice = unique(raw_data$data$mouse)
mice = mice[mice != 6553]

#average trees
dist_mat_list = map(mice, function(m) {
  message(m)
  x = raw_data
  x$data = raw_data_filtered$data %>% filter(mouse == m) %>%
    mutate(sample = tissue_group)
  mat = dist_mat(raw_data = x, dist_metric = 'manhattan')
  mat = as.matrix(mat)
  mat
})
map(dist_mat_list, function(mat) {
  rownames(mat)
})
names = rownames(as.matrix(dist_mat_list[[14]]))
dist_mat_list_to_avg = map(dist_mat_list, function(mat) {
  orig = mat
  mat = as.matrix(mat)
  if (length(intersect(rownames(orig), names)) == 7) {
    ordered_mat = mat[names, names]
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
plot(avg_phy)

avg_phy = upgma(as.dist(avg_mat))
avg_phy = name_nodes(avg_phy)
plot(avg_phy)

cosine_distance_matrix <- function(raw_data) {
  raw_data$mf_wide = raw_data$data %>%
    ungroup() %>%
    transmute(sample, ID_sequence = paste0(ID, "_", sequence), mosaic_fraction) %>% 
    group_by(sample, ID_sequence) %>% summarise(mosaic_fraction = mean(mosaic_fraction)) %>%
    spread(key = ID_sequence, value = mosaic_fraction, fill = 0)
  raw_data$mf_mat = as.matrix(raw_data$mf_wide[-1])
  rownames(raw_data$mf_mat) = raw_data$mf_wide[[1]]
  
  mf_ind = 1 * (raw_data$mf_mat > 0)
  computePairwiseCosine <- function(binaryMatrix) {
    numRows <- nrow(binaryMatrix)
    cosineMatrix <- matrix(1, nrow = numRows, ncol = numRows)
    for (i in 1:(numRows - 1)) {
      for (j in (i + 1):numRows) {
        intersection <- sum(binaryMatrix[i, ] & binaryMatrix[j, ])
        product_cardinalities <- sum(binaryMatrix[i, ]) * sum(binaryMatrix[j, ])
        if (product_cardinalities != 0) {
          cosineMatrix[i, j] <- intersection / product_cardinalities
          cosineMatrix[j, i] <- cosineMatrix[i, j]  # Matrix is symmetric
        }
      }
    }
    rownames(cosineMatrix) = colnames(cosineMatrix) = rownames(binaryMatrix)
    # Return the Jaccard matrix
    return(cosineMatrix)
  }
  cosine_dist = 1 - computePairwiseCosine(mf_ind)
  return(cosine_dist)
}




