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

raw_data = load_samples_true_pairs('~/Documents/R/Jessica/240111_Bulk_Brain_Files/')
raw_data$data = raw_data$data[raw_data$data$sample != 'PB21-founder_filteredpairs.txt',]

sample_metadata = as_tibble(stringr::str_match(unique(raw_data$data$sample), "(.*?)_(.*?)_(\\d+)_(\\d+_\\d+)-(.*?)_truepairs.txt"))[1:5]
colnames(sample_metadata) = c('sample', 'side', 'region', 'rep', 'tech_rep')
raw_data = append_metadata(raw_data, sample_metadata)
raw_data$data = raw_data$data %>% group_by(ID, sequence, region) %>%
  summarise(reads = sum(reads))
raw_data$data = raw_data$data[!is.na(raw_data$data$region),]
raw_data$data = raw_data$data %>% mutate(sample = region)

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

sample_id_ind_wide = ungroup(raw_data$data) %>% select(ID, sample) %>%
  distinct() %>%
  mutate(value = 0.) %>%
  spread(key = ID, value = value, fill = 1.)

IDs_to_filter = colnames(sample_id_ind_wide[,-1])[colSums(sample_id_ind_wide[,-1]) > 1]

raw_data$data = raw_data$data %>% filter(!(ID %in% IDs_to_filter))

#rest of raw_data processing
parent = load_parent('./../../David/PB21-founder_filteredpairs.txt')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)

raw_data = filter_unmutated(raw_data)

colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
raw_data_filtered = filter_recurring_spacers(raw_data, indelphi_table, 0.0001, filter_missing = FALSE)

mat = dist_mat(raw_data = raw_data, dist = 'jaccard')
phy = nj(as.dist(mat))
phy = name_nodes(phy)
phy = phytools::midpoint.root(phy)
plot(phy)

phy = upgma(as.dist(mat))
phy = name_nodes(phy)
phy = phytools::midpoint.root(phy)
plot(phy)

