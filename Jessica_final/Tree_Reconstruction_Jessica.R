library(ggpubr)
library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('~/Documents/qfm2')
devtools::load_all()

setwd("~/Documents/R/MosaicFractionAnalysis/R")
devtools::load_all()

setwd('~/Documents/R/Jasmine/')
indelphi_table <- readr::read_tsv('MARC1_indelphi_result.txt', col_names = T)

raw_data = load_samples_true_pairs('~/Documents/R/Jessica/240111_Bulk_Brain_Files/')
sample_metadata = as_tibble(stringr::str_match(unique(raw_data$data$sample), "(.*?)_(.*?)_(\\d+)_(\\d+_\\d+)-(.*?)_truepairs.txt"))[1:5]
colnames(sample_metadata) = c('sample', 'side', 'region', 'rep', 'tech_rep')
raw_data = append_metadata(raw_data, sample_metadata)
raw_data$data = raw_data$data %>% group_by(ID, sequence, side, region) %>%
  summarise(reads = sum(reads))
raw_data$data = raw_data$data[!is.na(raw_data$data$region),]
raw_data$data = raw_data$data %>% mutate(sample = paste0(region, '_', side))
parent = load_parent('./../../David/PB21-founder_filteredpairs.txt')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)
raw_data = filter_unmutated(raw_data)
raw_data$data = unique(raw_data$data)
raw_data$data = raw_data$data[!is.na(raw_data$data$sample),]

mat = dist_mat(data = raw_data$data, dist = 'jaccard')
phy = nj(as.dist(mat))
#phy = ape::as.phylo(hclust(as.dist(mat), method = 'average'))
#phy = ape::as.phylo(hclust(as.dist(mat), method = 'ward.D2'))
phy = name_nodes(phy)
phy = phytools::midpoint.root(phy)