# install.packages("uwot")
library(uwot)
library(ComplexHeatmap)
library(qfm)
setwd('~/Documents/qfm2/')
devtools::load_all()

set_branch_lengths <- function(phy, edge_lengths, root_length) {
  phy$edge.length = edge_lengths
  phy$root.edge = root_length
  phy
}

bft <- function(phy, edge_length_tb, root_edge_len) {
  start_node = names(list_dd_and_tips_mod2(phy)$dd)
  start_node = c('Root', start_node)
  list = map(start_node, function(n) {
    if (n == 'Root') {
      dd = phy$node.label[1]
    } else {
      dd = list_dd_and_tips_mod2(phy)$dd[[n]]
    }
    as_tibble(expand_grid(in_node = n, out_node = dd))
  })
  df = as.data.frame(do.call(rbind, list))
  df = left_join(df, edge_length_tb, by = c('in_node', 'out_node'))
  n = ncol(df)
  df[df$in_node == 'Root', (n-2):n] = c(0, root_edge_len, root_edge_len)
  df = nest(df, data = -in_node)
  df$start = 0
  df
}

# simulate_program <- function(edges, name, steps) {
#   for (node in edges$in_node) {
#     start = edges$start[edges$in_node == node]
#     data = edges$data[edges$in_node == node]
#     data[[name]] = map(data$edge_length, function(l) {
#       diffusion(steps, l, start)
#     })
#     while(!program_correlation(data[[name]])) {
#       data[[name]] = map(data$edge_length, function(l) {
#         diffusion(steps, l, start)
#       })
#     }
#   }
# }

traverse_edges <- function(edges, num_program) {
  for (node in edges$in_node) {
    print(node)
    if (node == 'Root') {
      starts = rep(0,num_program)
    } else {
      starts = edges$start[edges$in_node == node][[1]]
    }
    edges = simulate_edge(edges, num_program, node, steps, starts)
    node_data = edges$data[edges$in_node == node][[1]]
    for (j in 1:length(node_data$out_node)) {
      ind = which(colnames(node_data) == 'program1')
      start_vec = map(edges$data[edges$in_node == node][[1]][j,ind:(num_program+ind-1)], function(x) x[[1]][length(x[[1]])])
      next_node = node_data$out_node[[j]]
      edges$start[edges$in_node == next_node] = list(start_vec)
    }
  }
  edges
}

simulate_edge <- function(edges, num_program, node, steps, starts) {
  edges_original = edges
  for (i in 1:num_program) {
    print(i)
    name = paste0('program', i)
    node_data = edges$data[edges$in_node == node][[1]]
    df = simulate_program(node_data, name, steps, starts[i][[1]])
    if (is.null(df)) {
      edges = edges_original
      return(simulate_edge(edges, num_program, node, steps, starts))
    } else {
      edges$data[edges$in_node == node][[1]] = df
    }
  }
  return(edges)
}

simulate_program <- function(node_data, name, steps, start) {
  node_data[[name]] = map(node_data$edge_length, function(l) {
    diffusion(steps, l, start)
  })

  cor = program_correlation(node_data, name)

  while(any(cor)) {
    for (a in 1:nrow(node_data)) {
      k = 0
      while (cor[a] & k < 200) {
        l = node_data$edge_length[[a]]
        node_data[[name]][[a]] = diffusion(steps, l, start)
        cor = program_correlation(node_data, name)
        k = k+1
      }
      if (k >= 200) {
        return(NULL)
      }
    }
    # node_data[[name]] = map(node_data$edge_length, function(l) {
    #   diffusion(steps, l, start)
    # })
  }
  node_data
}

# program_correlation <- function(node_data, name) {
#   for (p in names(node_data)[-c(1,2)]) {
#     if (p != name) {
#       cors = map_dbl(1:nrow(node_data), function(a) {
#         cor(node_data[[p]][[a]], node_data[[name]][[a]])
#       })
#       if (any(cors > 0.5)) {
#         return(FALSE)
#       }
#     } else {
#       return(TRUE)
#     }
#   }
#   return(TRUE)
# } 

program_correlation <- function(node_data, name) {
  vec = as.vector(c())
  for (a in 1:nrow(node_data)) {
    cors = map_dbl(names(node_data)[-c(1:4)], function(p) {
      if (p == name) {
        return(0)
      }
      cor(node_data[[p]][[a]], node_data[[name]][[a]])
    })
    val = any(cors > 0.4)
    vec = c(vec, val)
  }
  vec
}

#set to 200 steps per day
diffusion <- function(steps, branch_length, start = 0) {
  steps = branch_length*200
  walk = rep(0,steps)
  velocity = rep(0,steps)
  velocity[1] = rnorm(1,0,0.2)
  eta = runif(1,0,0.5)
  if (start == 0) {
    walk[1] = log(runif(1,0,1.5))
  } else {
    walk[1] = start
  }
  
  for (t in 1:(steps-1)) {
    walk[t+1] = walk[t] + velocity[t]
    e = rnorm(1,0,branch_length/steps)
    velocity[t+1] = eta*velocity[t] + e
  }
  return(walk)
  # return(tibble(time = seq(branch_length/steps,
  #                          branch_length,
  #                          branch_length/steps), program = walk))
}

check_parallel_branch_correlation <- function(edges) {
  for (a in 1:(nrow(edges)-1)) {
    print(a)
    for (b in (a+1):nrow(edges)) {
      print(b)
      if (edges$in_node[a] == edges$in_node[b]) {
        cors = map_dbl(1:num_genes, function(g) {
          name = paste0('gene', g)
          cor(edges[[name]][[a]], edges[[name]][[b]])
        })
        while (mean(cors < 0) < 0.45) {
          gene_coefs[[name]] = stats::rgamma(num_program, 0.05)
          edges[[name]] = map(1:nrow(edges), function(n) {
            list = map(1:num_program, function(j) {
              w = gene_coefs[[name]][j]
              p = paste0('program', j)
              edges[[p]][[n]] * w
            })
            colSums(do.call(rbind, list))
          })
          cors = map_dbl(1:num_genes, function(g) {
            name = paste0('gene', g)
            cor(edges[[name]][[a]], edges[[name]][[b]])
          })
        }
      }
    }
  }
  edges
}


#run
setwd('~/Documents/R/MF_Signal_Simulation/lineage specific simulation/')
# load('../MF Code/phy.rda')
phy = readRDS('./phy_com_scaled.rds')
#phy = readRDS('./../gast_phylo.rds')
steps = 200
num_program = 10

#edge_length_tb = readRDS('./phy_com_edges.rds')[,-3]
edge_length_tb = readRDS('./phy_edges_len.rds')
colnames(edge_length_tb) = c('in_node', 'out_node', 'start_time', 'end_time')
edge_length_tb$edge_length = edge_length_tb$end_time - edge_length_tb$start_time
root_edge = list('Root', edge_length_tb$in_node[[1]], 0, edge_length_tb$start_time[[1]],
                 edge_length_tb$start_time[[1]])
edge_length_tb = rbind(root_edge, edge_length_tb)

#edges2 = bft(phy, edge_length_tb, min(edge_length_tb$start_time))
edges = nest(edge_length_tb, data = -in_node)
edges$start = 0
st = proc.time()
edges = traverse_edges(edges, num_program)
#edges2 = traverse_edges(edges2, num_program)
proc.time() - st
# edges$start_time = 0

edges = unnest(edges, cols = data)

for (n in edges$in_node) {
  if (n == 'Root') {
    edges$start_time[edges$in_node == n] = 0
  } else {
    edges$start_time[edges$in_node == n] = edges$start_time[edges$out_node == n] +
      edges$edge_length[edges$out_node == n]
  }
}


#make each edge program combination a table with time values added,
#unnest one at a time, and then plot with ggplot
edges$time = map2(edges$start_time,
                  edges$edge_length,
                  function(s, l) {
  seq(s + 1/200, s+l, 1/200)
})

plot_tb = unnest(edges, cols = time) %>% select(in_node, out_node, time)
for (i in 3:3) {
  name = paste0('program', i)
  plot_tb = plot_tb %>% mutate(edges %>% unnest(cols = name) %>% select(name))
}
plot_tb = plot_tb %>% gather(key = 'program', value = 'value', -c(in_node, out_node, time))

median_program = plot_tb %>% group_by(program) %>% summarise(median = median(value))
base_mean = log(10^(1)) - median_program$median
median_program$base = rnorm(1,base_mean,1)
median_program$base = 2.11

plot_tb = left_join(plot_tb, median_program)
plot_tb$new_value = exp(plot_tb$base + plot_tb$value)

programs = ggplot(data = plot_tb,
       aes(x = time, y = new_value, color = out_node)) +
  scale_colour_manual(values = color_values) +
  geom_line() + xlab('Time (Days)') + ylab('Relative Gene Expression (AU)') +
  theme_pubr() + theme(axis.text = element_text(size = 10),
                       axis.title = element_text(size = 12),
                       plot.title = element_text(size = 9),
                       legend.position = 'none')
  
programs

png(filename = '~/Documents/R/MF_Signal_Simulation/pics/program3.png', width = 4, height = 4, units = "in", res = 1000)
print(programs)
dev.off()

max_programs = map_dbl(1:100, function(i) {
  name = paste0('program', i)
  max(unlist(edges %>% select(name) %>% unnest(cols = name)))
})

num_genes = 250
gene_coefs = as.data.frame(tibble(gene1 = stats::rgamma(num_program, 0.05)))
for (i in 2:num_genes) {
  name = paste0('gene', i)
  gene_coefs[[name]] = stats::rgamma(num_program, 0.05)
}

#exponentiate relative gene expression value to make it positive
for (i in 1:num_genes) {
  name = paste0('gene', i)
  edges[[name]] = map(1:nrow(edges), function(n) {
    list = map(1:num_program, function(j) {
      w = gene_coefs[[name]][j]
      p = paste0('program', j)
      edges[[p]][[n]] * w
    })
    colSums(do.call(rbind, list))
  })
}

edges = check_parallel_branch_correlation(edges)

#actually just go straight from gene program; give each program a base value, plot the exponentiated value over time
#integrate it over a day and see if it's approximately 10^-3
#base + relative should be ~ ln(10^-5) (base value should be pretty negative, maybe around -12)
#max abs expression should prob not be more than 10^-4
for (i in 1:num_genes) {
  name = paste0('gene', i)
  base = exp(rnorm(1,0.8,1))
  edges[[name]] = map(1:nrow(edges), function(n) {
    edges[[name]][[n]] * base
  })
  too_high = map_lgl(1:nrow(edges), function(n) {
    any(edges[[name]][[n]] > 5000)
  })
  while (any(too_high)) {
    base = exp(rnorm(1,0.8,1))
    edges[[name]] = map(1:nrow(edges), function(n) {
      edges[[name]][[n]] * base
    })
    too_high = map(1:nrow(edges), function(n) {
      any(edges[[name]][[n]] > 5000)
    })
  }
}



plot_tb = unnest(edges, cols = time) %>% select(in_node, out_node, time)
for (i in 1:6) {
  name = paste0('gene', i)
  plot_tb = bind_cols(plot_tb,   edges %>% unnest(cols = name) %>% select(name) %>% rename(x = name) %>%
                        mutate(x = x/max(abs(x))))
}

plot_tb = plot_tb %>% gather(key = 'gene', value = 'expression_level', -c(in_node, out_node, time))

ggplot(data = plot_tb,
       aes(x = time, y = expression_level, color = out_node)) +
  geom_line() + facet_wrap(~ gene)

#cell sampling
set.seed(0)
list_cell_gene_matrices = map(1:nrow(edges), function(i) {
  steps = edges[i,]$edge_length * 200
  sampled_cells = sample(1:steps, 30, replace = F)
  list_cell_gene_rows = map(names(edges)[(4+num_program+4):ncol(edges)],
                            function(gene) {
                              vec = edges[[gene]][[i]][sampled_cells]
                            })
  matrix = do.call(cbind, list_cell_gene_rows)
  rownames(matrix) = sampled_cells + edges$start_time[i]*steps
  matrix
})

cell_gene_matrix = do.call(rbind, list_cell_gene_matrices)
cell_gene_matrix = apply(cell_gene_matrix, 2, function(x) (x - mean(x)) / sd(x))
pca = stats::prcomp(cell_gene_matrix)
p_comp = pca$x[,1:20]


umap_data = as.data.frame(umap(p_comp))
colnames(umap_data) = c('UMAP1', 'UMAP2')
# umap_data$branch = rep(as.character(edges$out_node), each = 50)
umap_data$time = as.numeric(rownames(cell_gene_matrix))
umap_data$cell_type = rep(as.character(edges$out_node), each = 30)

umap_time = ggplot(data = umap_data,
       aes(x = UMAP1,
           y = UMAP2,
           color = time)) + geom_point() +
  scale_colour_distiller()

umap_cell = ggplot(data = umap_data,
                   aes(x = UMAP1,
                       y = UMAP2,
                       color = cell_type)) + geom_point()

png(filename = '~/Documents/R/MF_Signal_Simulation/pics/umap_time.png', width = 4, height = 4, units = "in", res = 1000)
print(umap_time)
dev.off()

png(filename = '~/Documents/R/MF_Signal_Simulation/pics/umap_cell.png', width = 4, height = 4, units = "in", res = 1000)
print(umap_cell)
dev.off()

pca_data = as.data.frame(pca$x[,1:2])
pca_data$time = as.numeric(rownames(cell_gene_matrix))
pca_data$cell_type = rep(as.character(edges$out_node), each = 30)

cell_state_col = read_excel('./gas_states_no_conv_mod.xlsx') %>%
  select(name, color)
color_values = cell_state_col$color
names(color_values) = cell_state_col$name

pca_time = ggplot(data = pca_data,
       aes(x = PC1,
           y = PC2,
           color = time)) +
  geom_point() + scale_colour_distiller()
pca_cell = ggplot(data = pca_data,
                  aes(x = PC1,
                      y = PC2,
                      color = cell_type)) +
  geom_point() +
  scale_colour_manual(values = color_values) +
  labs('colour' = 'Cell Type') +
  theme_pubr() + theme(axis.text = element_text(size = 10),
                       axis.title = element_text(size = 12),
                       plot.title = element_text(size = 9),
                       legend.position = 'none')

pca_time + pca_cell

png(filename = '~/Documents/R/MF_Signal_Simulation/pics/pca_time.png', width = 4, height = 4, units = "in", res = 1000)
print(pca_time)
dev.off()

png(filename = '~/Documents/R/MF_Signal_Simulation/pics/pca_cell.png', width = 4, height = 4, units = "in", res = 1000)
print(pca_cell)
dev.off()



pca_subset = pca_data[pca_data$time < 3*200 & pca_data$time > 0*200,]
pca_subset_plot_03 = ggplot(data = pca_subset,
                  aes(x = PC1,
                      y = PC2,
                      color = cell_type,
                      alpha = time)) +
  geom_point() + scale_colour_manual(values = color_values) +
  scale_alpha_continuous() + theme_pubr() +
  theme(legend.position = 'right')
pca_subset_plot_03

pca_subset = pca_data[pca_data$time < 6*200 & pca_data$time > 3*200,]
pca_subset_plot_36 = ggplot(data = pca_subset,
                            aes(x = PC1,
                                y = PC2,
                                color = cell_type,
                                alpha = time)) +
  geom_point() + scale_colour_manual(values = color_values) +
  scale_alpha_continuous() + theme_pubr() +
  theme(legend.position = 'right')
pca_subset_plot_36

pca_subset = pca_data[pca_data$time < 9*200 & pca_data$time > 6*200,]
pca_subset_plot_69 = ggplot(data = pca_subset,
                            aes(x = PC1,
                                y = PC2,
                                color = cell_type,
                                alpha = time)) +
  geom_point() + scale_colour_manual(values = color_values) +
  scale_alpha_continuous() + theme_pubr() +
  theme(legend.position = 'right')
pca_subset_plot_69

pca_subset_plot_03 + pca_subset_plot_36 + pca_subset_plot_69

png(filename = '~/Documents/R/MF_Signal_Simulation/lineage specific simulation/pca_over_time.png', width = 15, height = 6, units = "in", res = 1000)
print(pca_subset_plot_03 + pca_subset_plot_36 + pca_subset_plot_69)
dev.off()


library(Rtsne)
cell_gene_matrix_mat = cell_gene_matrix
cell_gene_matrix = as.data.frame(cell_gene_matrix)
tsne_data <- Rtsne(cell_gene_matrix, check_duplicates = F)
tsne_data = as.data.frame(tsne_data$Y, col.names = c('t-SNE 1', 't-SNE 2'))
tsne_data$time = as.numeric(rownames(cell_gene_matrix_mat))
tsne_data$cell_type = rep(as.character(edges$out_node), each = 200)

tsne_time = ggplot(data = tsne_data,
       aes(x = V1,
           y = V2,
           color = time)) +
  geom_point() + scale_colour_distiller() +
  xlab('t-SNE 1') + ylab('t-SNE 2')

tsne_cell = ggplot(data = tsne_data,
              aes(x = V1,
                  y = V2,
                  color = cell_type)) +
  geom_point() + xlab('t-SNE 1') + ylab('t-SNE 2')

png(filename = '~/Documents/R/MF_Signal_Simulation/pics/tsne_time.png', width = 4, height = 4, units = "in", res = 1000)
print(tsne_time)
dev.off()

png(filename = '~/Documents/R/MF_Signal_Simulation/pics/tsne_cell.png', width = 4, height = 4, units = "in", res = 1000)
print(tsne_cell)
dev.off()



load('./../mutation simulation output/count_graph_1.rda')





#pca at different time points
data = edge_tb %>% select(to_type, program3) %>% unnest(cols = program3)
base = 2.11
data$base = base
data$new_value = exp(data$base + data$value)
data_1 = data %>% filter(abs(time - 5) <= 0.001)
data_2 = data %>% filter(abs(time - 8) <= 0.001)

data_subset = data %>% filter(row_number() %% 1000 == 0)

cell_program = ggplot(data = data_subset) +
  geom_point(aes(x = time, y = new_value, color = to_type)) +
  scale_color_manual(values = color_values) + xlab('Time (Days)') +
  ylab('Relative Gene Expression (AU)') +
  theme_pubr() + theme(axis.text = element_text(size = 10),
                       axis.title = element_text(size = 12),
                       plot.title = element_text(size = 9),
                       legend.position = 'none')

cell_program

png(filename = '~/Documents/R/MF_Signal_Simulation/pics/cell_program.png', width = 4, height = 4, units = "in", res = 1000)
print(cell_program)
dev.off()





