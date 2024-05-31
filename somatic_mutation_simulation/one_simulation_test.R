sink(stdout(), type="message")

library(tidyverse)
library(r2r)
library(qfm)
library(furrr)
#library(rgl)
#library(treespace)
library(TreeDist)
setwd('/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo')
print('success')

# get job id from system 
job_id = 49

# construct parameter table 

param_tb <- read.table('/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/output/param_tb.txt',header = T)
print('success')

#load file
filename <- paste('./input/yi_output/', param_tb$input_file[job_id], sep = "")
print('success')
print(filename)

if (file.exists(filename)) {
    load(filename)
} else {
    print('file not found')
    stop()
}
print('success')
#initialize parameters
total_bp <- 10^9
mut_rate <- param_tb$mutation_rate[job_id]
num_mut <- 0
mut_list <- c()
node_mut_list <- hashmap()

in_node <- count_graph$phylo_edges$in_node
print('success')
edge_length <- count_graph$phylo_edges$length
node_mut_list[[names(in_node)[1]]] <- character()
i <- 1

#generate mutations
for (d in names(in_node)) {
  node <- in_node[[d]]
  if (is.null(node_mut_list[[node]])) {
    node_mut_list[[node]] <- character()
  }
  node_mut_list[[d]] <- node_mut_list[[node]]
  for (i in 1:rbinom(1,total_bp,1-exp(-mut_rate * edge_length[i]))) {
    mut <- paste0(sample(c(letters, 1:9), size = 20, replace = T), collapse = "")
    mut_list <- c(mut_list, mut)
    node_mut_list[[d]] <- c(node_mut_list[[d]], mut)
    num_mut <- num_mut + 1
  }
  i <- i+1
}

#make cell mutation table of tips
cell_mut_tb = bind_rows(map(tr$tip.label, function(tip_id) {
  tibble(cell = tip_id,
         mut = node_mut_list[[tip_id]],
         value = 1)
}))

#subset to mutations that occur in more than one cell
occurrences <- table(cell_mut_tb$mut)
subset_mut <- occurrences[occurrences > 1]
cell_mut_tb <- subset(cell_mut_tb, mut %in% names(subset_mut))

#formatting and reconstruction
m <- spread(cell_mut_tb, mut, value, fill = 0)
m$type <- m$cell
m <- m[,c(1,ncol(m),3:ncol(m))]

chr_mat = as.matrix(m[-c(1:2)])
rownames(chr_mat) = m$cell

plan(multisession, workers = 16)
print('success')

mut_p = estimate_mut_p(chr_mat, t_total = res2$total_time - tr$root.edge)
mut_p$mut_rate = rep(mut_rate, ncol(m)-2)

tr2 = phylotime(chr_mat, t_total = res2$total_time - tr$root.edge)
print('phylotime done')
tr2$root.edge <- res2$total_time - tr$root.edge

#calculate distances
kc0 <- KendallColijn(tr,tr2)
print('success')

#save results
result <- paste(param_tb$sample_size[job_id], param_tb$sampling[job_id], mut_rate, filename, job_id, num_mut, kc0, kc0, sep = '  ')
print('success')
system(paste("echo ",result,' >> /home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/output/results.txt', sep = ""))
output_filename <- paste('./output/', param_tb$input_file[job_id], '_', mut_rate, '_', param_tb$num_sim_2[job_id], '.rda', sep = "")
print(output_filename)
save(tr2, list = c("tr2"), file = output_filename)