sink(stdout(), type="message")

library(tidyverse)
library(r2r)
library(qfm)
library(furrr)
#library(rgl)
#library(treespace)
library(TreeDist)
library(igraph)
library(ggraph)
setwd('/home/ssrikan2/data-kreza1/smriti/qfm2')
devtools::load_all()

setwd('/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo')

job_id = 49
#job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# construct parameter table 

param_tb <- read.table('/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/correctness_output/param_tb_2.txt',header = T)
filename <- paste('./input/yi_output/', param_tb$input_file[job_id], sep = "")
if (file.exists(filename)) {
    print(filename)
    load(filename)
} else {
    print('file not found')
    stop()
}

#functions

l1 <- p1$PS
l2 <- p2$`Node-7`
compare_lists(l1,l2)

compare_lists <- function (l1, l2) {
  if ((length(intersect(l1,l2)) == length(l1)) && (length(intersect(l1,l2)) == length(l2))) {
    return(TRUE)
  }
  return(FALSE)
}

get_partitions <- function(gr) {
        obj = list_dd_and_tips_mod2(gr)
        dd_list = obj$dd
        tips_list = obj$tips
        gr_tips = gr$tip.label
        names(gr_tips) = gr_tips
        tips_list = c(tips_list, gr_tips)
        part_list = map(dd_list, function(dd_vec) {
                tips_list[dd_vec]
        })
        part_list
}

setwd('~/Documents/R')
load('sample_size_200_fixed_83.rda')
setwd('~/Documents/qfm2')
phy = readRDS("intermediate_data/gast_phylo.rds")
devtools::load_all()

t1 <- phy
t2 <- res1$gr

plot(phy)
plot(res1$gr)

compare_partitions <- function (t1, t2) {
  p1 <- get_partitions(t1)
  p2 <- get_partitions(t2)
  overall_count = 0
  for (a in 1:length(p1)) {
    for (b in 1:length(p2)) {
      node1 = p1[[a]]
      node2 = p2[[b]]
      if (length(node1) != length(node2)) {
        next
      }
      count = 0
      for (i in 1:length(node1)) {
        for (j in 1:length(node2)) {
          l1 <- node1[[i]]
          l2 <- node2[[j]]
          
          if (compare_lists(l1,l2)) {
            count <- count + 1
          }
        }
      }
      if (count == length(node1)) {
        overall_count <- overall_count + 1
      }
    }
  }
  return(c(overall_count, length(p1), length(p2)))
}

compare_partitions(t1,t2)

#calculations
p <- get_partitions(phy)
p0 <- get_partitions(res0$gr)
p1 <- get_partitions(res1$gr)

r0 <- compare_partitions(p, p0)
r1 <- compare_partitions(p, p1)

print(r0[1])
print(r1[1])
print(r0[3])
print(r1[3])

r0_correct <- r0[1]/r0[3]
r1_correct <- r1[1]/r1[3]
r0_complete <- r0[1]/r0[2]
r1_complete <- r1[1]/r1[2]

correct_tree = ape::read.tree(text = "(((he,ta),yo),pl);")
drop.tip(correct_tree, tip = "ta")



#save results
result0 <- paste(param_tb$sample_size[job_id], param_tb$sampling[job_id], filename, 'res0', job_id, r0_correct, r0_complete, sep = '  ')
result1 <- paste(param_tb$sample_size[job_id], param_tb$sampling[job_id], filename, 'res1', job_id, r1_correct, r1_complete, sep = '  ')

system(paste("echo ",result0,' >> /home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/correctness_output/results_2.txt', sep = ""))
system(paste("echo ",result1,' >> /home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/correctness_output/results_2.txt', sep = ""))
