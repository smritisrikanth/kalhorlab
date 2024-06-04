#' calculates average density across several densities
#'
#' @param density_list list of density objects with the same set of bins
#' @param colname name of column to average
#' @return a single density representing the average across the list
average_density <- function(density_list, colname = "count") {
  matrix <- as.vector(c())
  matrix_cum <- as.vector(c())
  for (i in density_list) {
    matrix <- cbind(matrix, i$table[[colname]])
    matrix_cum <- cbind(matrix_cum, i$table$cum_mf)
  }
  matrix <- as.numeric(matrix)
  matrix <- matrix(matrix, ncol = length(matrix), nrow = length(density_list[[1]]$table[[colname]]))
  matrix_cum <- as.numeric(matrix_cum)
  matrix_cum <- matrix(matrix_cum, ncol = length(matrix_cum), nrow = length(density_list[[1]]$table$cum_mf))
  avg_counts = tibble(bin = density_list[[1]]$table$bin,
                      count = rowMeans(matrix), cum_mf = rowMeans(matrix_cum))
  stdv = apply(matrix, 1, sd)
  se = stdv/sqrt(ncol(matrix))
  avg_density = structure(list(table = avg_counts, se = se), class = 'mf_density')
  avg_density
}

#' assigns each unique mutation in list of mutation (ID sequence combination) tip combinations with node in gr
#'
#' @param mut_tip_list list of mutation tip combinations
#' @param gr tree from which nodes are to be assigned to mutation tip combinations
#' @return table of mutation tip combinations and assigned nodes
assign_mut_node <- function(mut_tip_list, gr) {
  gr_tip_list = list_dd_and_tips_mod2(gr)$tips
  gr_tips = as.list(gr$tip.label)
  names(gr_tips) = gr_tips
  gr_tip_list = c(gr_tip_list, gr_tips)
  gr_tip_size = map_dbl(gr_tip_list, length)
  out = bind_rows(map(1:length(mut_tip_list), function(j) {
    state_vec = mut_tip_list[[j]]
    candid_states = names(gr_tip_list)[map_lgl(gr_tip_list,
                                               function(y) all(state_vec %in% y))]
    candid_sizes = gr_tip_size[candid_states]
    assertthat::assert_that(sum(candid_sizes == min(candid_sizes)) ==
                              1)
    temp_state = candid_states[which.min(candid_sizes)]
    assertthat::assert_that(length(temp_state) == 1)
    tibble(mut = names(mut_tip_list)[j],
           gr_node = temp_state)
  }))
  out
}

#' creates ks_test_result object with results of two sample KS test between two mosaic fraction distributions
#'
#' @param mf_distribution1 mf_distribution
#' @param mf_distribution2 mf_distribution to be compared to mf_distribution1
#' @return ks_test_result object with results of two sample KS test
ks_test <- function(mf_distribtion1, mf_distribution2) {
  ks_test_result
}
