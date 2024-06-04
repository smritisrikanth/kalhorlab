#' loads sample data into a raw object
#'
#' @param file_paths string path to location of files
#' @return raw holding raw data for each allele
#' @export
load_samples <- function(file_paths) {
  setwd(file_paths)
  all_files = list.files(file_paths, pattern = "filteredpairs\\.txt")
  all_files = all_files[-1]
  all_samples = as_tibble(all_files)
  all_samples$data = map(all_samples$value, function(fn) {
    readr::read_tsv(paste0("./", fn), col_names = c('ID', 'sequence', 'reads', '%ID', '%total'))
  })
  all_samples$id_summary = map(all_samples$data, function(tb) {
    tb %>% group_by(ID) %>%
      summarise(total = sum(reads)) %>%
      mutate(frac = total / sum(total))
  })
  all_samples_sequences = transmute(all_samples, sample = value, data) %>%
    unnest(cols = data)

  raw <- structure(list(data = all_samples_sequences, MFcomputed = FALSE), class = 'raw')
  raw
}

#' loads sample data into a raw object, true pairs version
#'
#' @param file_paths string path to location of files
#' @return raw holding raw data for each allele
#' @export
load_samples_true_pairs <- function(file_paths) {
  setwd(file_paths)
  all_files = list.files(file_paths, pattern = "truepairs\\.txt")
  all_files = all_files[-1]
  all_samples = as_tibble(all_files)
  all_samples$data = map(all_samples$value, function(fn) {
    readr::read_tsv(paste0("./", fn), col_names = c('ID', 'sequence', 'reads', '%ID', '%total'))
  })
  all_samples$id_summary = map(all_samples$data, function(tb) {
    tb %>% group_by(ID) %>%
      summarise(total = sum(reads)) %>%
      mutate(frac = total / sum(total))
  })
  all_samples_sequences = transmute(all_samples, sample = value, data) %>%
    unnest(cols = data)

  raw <- structure(list(data = all_samples_sequences, MFcomputed = FALSE), class = 'raw')
  raw
}

#' loads parent data into a parent_table object
#'
#' @param filepath string path to location of file
#' @return parent_table holding parent data
#' @export
load_parent <- function(filepath) {
  table = read.csv(filepath, sep = '\t', col.names = c('ID', 'sequence', 'reads', '%ID', '%total'))
  parent_table = structure(list(table = table), class = 'parent_table')
}

#' filters IDs based on existence within founder set of ID sequence combinations
#'
#' @param raw raw holding raw data for each allele, including the following columns: ID, sequence
#' @param parent_table parent_table holding founder set of IDs and spacer sequences, including the following columns: ID, sequence
#' @return raw holding raw data for each allele
#' @export
id_filter <- function(raw, parent_table) {
  missing_IDs <- unique(raw$data$ID[!(raw$data$ID %in% parent_table$table$ID)])
  raw$data <- raw$data[!(raw$data$ID %in% missing_IDs),]
  raw$filtered_IDs <- missing_IDs
  message('filtered IDs:')
  message(list(missing_IDs))
  raw
}

#' identifies mutated alleles based on existence within parental set of spacers
#'
#' @param raw parent_table holding raw data for each allele, including the following columns: ID, sequence
#' @param parent_table parent_table holding founder set of IDs and spacer sequences, including the following columns: ID, sequence
#' @return raw holding raw data for each allele, including the following added column: mutated
#' @export
annotate_parental_spacers <- function(raw, parent_table) {
  muts = paste0(raw$data$ID, '_', raw$data$sequence)
  muts_parent = paste0(parent_table$table$ID, '_', parent_table$table$sequence)
  raw$data$mutated = !(muts %in% muts_parent)
  raw
}

#' filters ID spacer combinations based on indelphi calculated emergence probability
#'
#' @param raw raw holding raw data for each allele, including the following columns: ID, sequence
#' @param emergence_table emergence_table holding founder set of IDs and spacer sequences, including the following columns: ID, sequence, probability
#' @param prob_cutoff double probability cutoff to subset mutations by
#' @param include_missing boolean indicating whether alleles not found in emergence_table should be filtered out
#' @return raw holding raw data for each allele, including the following added column: probability
#' @export
filter_recurring_spacers <- function(raw, emergence_table, prob_cutoff, filter_missing = TRUE) {
  emergence_table_subsetted = emergence_table %>% transmute(ID = ID,
                                                 sequence = sequence,
                                                 probability = probability)
  emergence_table_subsetted = emergence_table_subsetted %>%
    group_by(ID, sequence) %>%
    summarise(probability = sum(probability))
  raw$data <- left_join(raw$data, emergence_table_subsetted, by = c('ID', 'sequence'))
  if (filter_missing) {
    filtered_recurring_alleles = cbind(raw$data$ID[!is.na(raw$data$probability) & raw$data$probability <= prob_cutoff],
                                       raw$data$sequence[!is.na(raw$data$probability) & raw$data$probability <= prob_cutoff])
    raw$data = raw$data[!is.na(raw$data$probability) & raw$data$probability <= prob_cutoff,]
  } else {
    filtered_recurring_alleles = cbind(raw$data$ID[is.na(raw$data$probability) | raw$data$probability <= prob_cutoff],
                                       raw$data$sequence[is.na(raw$data$probability) | raw$data$probability <= prob_cutoff])
    raw$data = raw$data[is.na(raw$data$probability) | raw$data$probability <= prob_cutoff,]
    raw$data$probability[is.na(raw$data$probability)] = min(emergence_table_subsetted$probability)
  }
  raw$filtered_recurring_alleles = filtered_recurring_alleles
  message('filtered recurring alleles')
  #message(list(filtered_recurring_alleles))
  raw
}

#' calculates mosaic fraction for each mutated allele in each sample
#'
#' @param raw raw holding raw data for each allele, including the following columns: sample, ID, sequence, reads
#' @return raw holding raw data for each allele, including the following added columns: mosaic_fraction and log2MF
#' @export
compute_mosaic_fraction <- function(raw) {
  raw$data = raw$data %>% group_by(ID, sample) %>%
    mutate(mosaic_fraction = reads/sum(reads))
  raw$data$mosaic_fraction[!raw$data$mutated] = NA
  raw$data$log2MF = log2(raw$data$mosaic_fraction)
  raw$MFcomputed = TRUE
  raw
}

#' groups raw data by ID and sample
#'
#' @param raw raw holding raw data for each allele, including the following columns: sample, ID
#' @return  nested_data holding raw data for each allele in both raw and nested forms, nested by ID and sample with all additional data in dataframe data
#' @export
nest_raw <- function(raw, by_var) {
  nested <- nest(raw$data, data = -by_var)
  nested_data <- structure(list(raw = raw, nested = nested), class = 'nested_data')
  nested_data
}

#' calculates number of unique mutations and mutated fraction for each ID-sample combination
#'
#' @param nested_data nested_data object
#' @return nested_data object, including the following added columns: num_mut, mutated_fraction
#' @export
compute_mutation_information <- function(nested_data) {
  nested_data$nested$num_mut = map_dbl(nested_data$nested$data, function(data) {
    nrow(data)
  })
  nested_data$nested$mutated_fraction = map_dbl(nested_data$nested$data, function(data) {
    sum(data$reads[data$mutated])/sum(data$reads)
  })
  message('overall mutated fraction')
  message(mean(nested_data$nested$mutated_fraction))
  nested_data
}

#' filters unmutated alleles out of raw data for mosaic fraction analysis
#'
#' @param raw raw holding raw data for each allele, including the following columns: mutated
#' @return raw holding both filtered and unfiltered versions of raw data for each allele
#' @export
filter_unmutated <- function(raw) {
  raw$unfiltered <- raw$data
  raw$data <- raw$data[raw$data$mutated,]
  raw$filtered_mutated = TRUE
  raw
}

#' calculates mosaic fraction density curve given mosaic fraction data and desired bins
#'
#' @param mf_vec vector of mosaic fraction values to create density curve to represent the distribution of mosaic fractions
#' @param breaks vector of values to define desired bins to separate mosaic fraction values into
#' @return mf_density representing density curve
#' @export
calculate_density <- function(mf_vec, breaks) {
  mf_vec = as.numeric(mf_vec)
  bin_sizes = diff(breaks)
  bin_sizes[is.infinite(bin_sizes)] = min(bin_sizes)
  midpoints = rev(breaks[2:length(breaks)] - (bin_sizes/2))

  mf_bin = cut(log2(mf_vec), breaks = breaks)
  mf_count = as.numeric(rev(table(mf_bin)))
  mf_in_bin <- rep(c(0), length(mf_count))
  for (mf in mf_vec) {
    for (i in 1:length(midpoints)) {
      if (log2(mf) <= (midpoints[i] + bin_sizes[i]/2) && log2(mf) > (midpoints[i] - bin_sizes[i]/2)) {
        mf_in_bin[i] <- mf_in_bin[i] + mf
        break
      }
    }
  }
  counts = tibble(bin = midpoints, raw_count = mf_count,
         count = mf_count/sum(mf_count)/bin_sizes,
         mf_in_bin = mf_in_bin, cum_mf = cumsum(mf_in_bin))
  mf_density <- structure(list(table = counts, breaks = breaks), class = 'mf_density')
  mf_density
}

#' #' calculates mosaic fraction density curve given mosaic fraction data and existing calculated density
#' #'
#' #' @param mf_vec vector of mosaic fraction values to create density curve to represent the distribution of mosaic fractions
#' #' @param mf_density mf_density object representing previously calculated density curve
#' #' @return mf_density representing new density curve with old cumulative mosaic fraction information
#' recalculate_density <- function(mf_vec, mf_density) {
#'   counts = calculate_density(mf_vec, mf_density$breaks)$table
#'   mf_density$table$count = counts$count
#'   mf_density
#' }

#' calculates mosaic fraction distribution given mosaic fraction data and desired bins
#'
#' @param mf_vec vector of mosaic fraction values to create density curve to represent the distribution of mosaic fractions across bins
#' @param breaks vector of values to define desired bins to separate mosaic fraction values into
#' @return mf_distribution representing distribution of mosaic fractions across bins
#' @export
calculate_distribution <- function(mf_vec, breaks) {
  return(mf_distribution)
}

#' appends a density object representing the mosaic fraction curve for each ID-sample combination
#'
#' @param nested_data nested_data object
#' @param breaks vector of values to define desired bins to separate mosaic fraction values into
#' @return nested_data object with density column added
#' @export
append_density <- function(nested_data, breaks) {
  nested_data$nested$density = map(nested_data$nested$data, function(data) {
    calculate_density(data$mosaic_fraction, breaks)
  })
  nested_data
}

#' adjusts existing density object representing the mosaic fraction curve for each ID-sample combination
#'
#' @param nested_data_filtered nested_data object with an existing density column
#' @param nested_data nested_data object with an existing density column
#' @return nested_data object with density column with adjusted counts
#' @export
adjust_density <- function(nested_data_filtered, nested_data, by_var) {
  colnames(nested_data_filtered$nested)[colnames(nested_data_filtered$nested) == 'density'] = 'filtered_density'
  nested_data_rename = nested_data$nested %>% select(by_var, 'density')
  nested_data_filtered$nested = left_join(nested_data_filtered$nested, nested_data_rename, by = by_var)
  nested_data_filtered$nested$density = map2(nested_data_filtered$nested$filtered_density, nested_data_filtered$nested$density, function(filtered_density, density) {
    density$table$unfiltered_count = density$table$count
    density$table$count = filtered_density$table$count
    density$table$raw_count = filtered_density$table$raw_count
    density
  })
  nested_data_filtered
}

#' calculates normalized mosaic fraction density curve given previously calculated density
#'
#' @param mf_density mf_density object with latest_count column
#' @return mf_density representing same density with norm_count column
#' @export
normalize_density <- function(mf_density, apply_cum_adjustment = TRUE) {
  count_vec = mf_density$table$raw_count
  bin = mf_density$table$bin
  cum_mf = mf_density$table$cum_mf
  norm_counts_vec = rep(0, length(count_vec))
  for (i in 1:length(count_vec)) {
    if (count_vec[i] == 0) {
      norm_counts_vec[i] <- 0
    } else if (i == 1) {
      norm_counts_vec[i] <- count_vec[i]/(2^(-1*bin[i]))
    } else {
      if (apply_cum_adjustment) {
        norm_counts_vec[i] <- count_vec[i]/(1-cum_mf[i-1])/(2^(-1*bin[i]))
      } else {
        norm_counts_vec[i] <- count_vec[i]/(2^(-1*bin[i]))
      }
    }
  }
  mf_density$table$norm_count = norm_counts_vec
  mf_density
}

#' calculates normalized mosaic fraction density curves given previously calculated densities
#'
#' @param nested_data nested_data object with adjusted density column
#' @return nested_data with normalized counts added to each mf_density in density column
#' @export
normalize_density_list <- function(nested_data, apply_cum_adjustment = TRUE) {
  nested_data$nested$density = map(nested_data$nested$density, function(density) {
    normalize_density(density, apply_cum_adjustment)
  })
  nested_data
}

#' calculates average density for each unique combination of a given set of variables
#'
#' @param nested_data nested_data object including the following columns: density, all columns in by_var
#' @param by_var columns by which to nest
#' @param normalized logical value indicating whether to average normalized counts instead of counts
#' @param unfiltered logical value indicating whether to average unfiltered counts instead of counts
#' @return nested_density object with average_density column
#' @export
average_density_by <- function(nested_data, by_var, normalized = FALSE, unfiltered = FALSE, cum_mf = FALSE, raw = FALSE) {
  if (normalized) {
    col_to_avg = 'norm_count'
  } else if (unfiltered) {
    col_to_avg = 'unfiltered_count'
  } else if (cum_mf) {
    col_to_avg = 'cum_mf'
  } else if (raw) {
    col_to_avg = 'raw_count'
  }
    else {
    col_to_avg = 'count'
  }
  nested_density = nest(nested_data$nested, density_list = -by_var)
  nested_density$avg_density = map(nested_density$density_list, function(density_list) {
    average_density(density_list$density, col_to_avg)
  })
  nested_density_data = structure(list(nested = nested_data, nested_density = nested_density), class = 'nested_density')
  nested_density_data$averaged_col = col_to_avg
  nested_density_data
}

#' plots averaged curves in a nested density object grouped, colored, and faceted by specified variables
#'
#' @param nested_density nested_density object NOT WORKING
#' @param group variable to group by
#' @param color variable to color by
#' @param facet variable to facet by
#' @param color_values manually specified color scale to plot curves by if desired
#' @export
plot_averaged_curves <- function(nested_density,
                                 group_str = NULL,
                                 color_str = NULL,
                                 facet_str = NULL,
                                 color_values = NULL) {
  plot_tb = nested_density$nested_density %>% select(c(group_str, color_str, facet_str, 'avg_density'))
  plot_tb$table = map(plot_tb$avg_density, function(density) {
    density$table
  })
  plot_tb = unnest(plot_tb, cols = table)
  ylabel = paste0('Mosaic Fraction Density: ', nested_density$averaged_col)
  if (!is.null(facet_str)) {
    facet = paste0('~', facet_str)
    print(ggplot(data = plot_tb,
                 aes_string(x = "bin",
                            y = "count",
                            group = group_str,
                            color = color_str)) +
            scale_colour_manual(values = color_values) +
            facet_wrap(as.formula(facet)) +
            geom_line() +
            theme(legend.position = "right")) +
            ylab(ylabel)
  } else {
    print(ggplot(data = plot_tb,
                 aes_string(x = "bin",
                            y = "count",
                            group = group_str,
                            color = color_str)) +
            geom_line() +
            theme(legend.position = "right")) +
            ylab(ylabel) +
            ylim(0,1)
  }
}

#' plots heatmap depicting binary detection of each ID within each sample
#'
#' @param raw raw object with ID and sample columns
#' @param sample_annotation vector of variables by which to annotate heatmap rows (samples)
#' @export
plot_id_sample_heatmap <- function(raw, sample_annotation, ...) {
  sample_id_ind_wide = select(raw$data, ID, sample) %>%
    distinct() %>%
    mutate(value = 1.) %>%
    spread(key = ID, value = value, fill = 0.)

  temp = select(raw$data, ID, sample, !!(sample_annotation)) %>%
    distinct()
  row_annt_args = map(sample_annotation, function(annt_chr) {
    out_vec = temp[[annt_chr]]
    names(out_vec) = temp$sample
    out_vec
  })
  names(row_annt_args) = sample_annotation

  sample_id_ind_mat = as.matrix(sample_id_ind_wide[-1])
  rownames(sample_id_ind_mat) = sample_id_ind_wide[[1]]
  Heatmap(sample_id_ind_mat,
          show_row_names = F,
          show_column_names = F,
          left_annotation = do.call(rowAnnotation, map(row_annt_args, function(x) {
            x[rownames(sample_id_ind_mat)]
          })),
          ...)

}

#' appends any desired metadata to raw data based on specified columns
#'
#' @param raw raw holding raw data for each allele, including the following columns: cols
#' @param metadata metadata with columns of interest to be appended to raw data, including the following columns: cols
#' @param join_by vector of column names to join (left join) metadata by
#' @param cols columns of interest in metadata to append to raw data
#' @return raw holding raw data for each allele with all additional columns in metadata
#' @export
append_metadata <- function(raw, metadata, join_by = NULL, cols = NULL) {
  if (!is.null(cols)) {
    metadata = metadata[,c(join_by, cols)]
  }
  raw$data = left_join(raw$data, metadata, by = join_by)
  raw
}

#' assigns node from gr to each unique mutation in raw data
#'
#' @param raw raw object with ID, sequence, and tissue columns present in data
#' @param gr tree from which nodes are to be assigned
#' @return raw object with addition column gr_node in data
#' @export
annotate_mut_node <- function(raw, gr) {
  raw$data$ID_sequence = paste0(raw$data$ID, "_", raw$data$sequence)
  mut_tip_list = split(raw$data$tissue, raw$data$ID_sequence)
  raw$data = left_join(raw$data,
                     dplyr::rename(assign_mut_node(mut_tip_list, gr), ID_sequence = mut))
  raw
}


#' increments cumulative mutated fraction values for daughter nodes
#'
#' @param nested_data nested_data object
#' @param node_order_list list of nodes numbered numerically from root(1) to tips
#' @return nested_data object with updated cumulative mutated fraction values
#' @export
increment_cum_mf <- function(nested_data, node_order_list, nested_by_node = TRUE, by_var = NULL) {
  if (!nested_by_node) {
    nested_data$nested$data = map(nested_data$nested$data, function(data) {
      data$gr_node_num = map_dbl(data$gr_node, function(node) {node_order_list[[node]]})
      data
    })
    nested_data$nested = unnest(nested_data$nested, cols = data)
  } else {
    nested_data$nested$gr_node_num = map_dbl(nested_data$nested$gr_node, function(node) {node_order_list[[node]]})
  }
  nested_id_sample = nest(nested_data$nested, outer_data = -c(ID, sample))
  nested_id_sample$outer_data = map(nested_id_sample$outer_data, function(outer_data) {
    outer_data = outer_data[order(outer_data$gr_node_num),]
    if (nrow(outer_data) == 1) {
      return(outer_data)
    }
    cum_mf_total = outer_data$density[[1]]$table$cum_mf[length(outer_data$density[[1]]$table$cum_mf)]
    for (i in 2:nrow(outer_data)) {
      if (nrow(outer_data) < i) {
        print(outer_data)
      }
      outer_data$density[[i]]$table$cum_mf = outer_data$density[[i]]$table$cum_mf + cum_mf_total
      cum_mf_total = outer_data$density[[i]]$table$cum_mf[length(outer_data$density[[1]]$table$cum_mf)]
    }
    outer_data
  })
  nested_data$nested = unnest(nested_id_sample, cols = outer_data)
  if (!nested_by_node) {
    nested_data$nested = nest(nested_data$nested, data = -c(by_var, density))
  }
  nested_data
}

#' NEEDS TO BE UPDATED
#'
#' @param nested_data nested_data object
#' @param node_order_list list of nodes numbered numerically from root(1) to tips
#' @return nested_data object with updated cumulative mutated fraction values
#' @export
convert_mf_to_time <- function(nested_density) {
  nested_density$nested_density$lm = map(nested_density$nested_density$density_list, function (list) {
    log2mf = as.vector(c())
    node_time = as.vector(c())
    for (tb in list$data) {
      log2mf = c(log2mf, tb$avg_log2mf)
      node_time = c(node_time, (tb$node_to_time + tb$node_from_time)/2)
    }
    lm(formula = node_time ~ log2mf)
  })
  nested_density$nested_density$avg_density = map2(nested_density$nested_density$avg_density,
                                                   nested_density$nested_density$lm,
                                                   function(density, lm) {
                                                     density$table$time = predict(lm, newdata = list(log2mf = density$table$bin))
                                                     density
                                                   })
  nested_density
}

#' NEEDS TO BE UPDATED
#'
#' @param nested_data nested_data object
#' @param node_order_list list of nodes numbered numerically from root(1) to tips
#' @return nested_data object with updated cumulative mutated fraction values
#' @export
convert_mf_to_time_with_bias_correction <- function(nested_density) {
  nested_density$nested_density$lm = map2(nested_density$nested_density$density_list,
                                          nested_density$nested_density$commitment_bias, function (list, b) {
    log2mf = as.vector(c())
    node_time = as.vector(c())
    for (tb in list$data) {
      log2mf = c(log2mf, tb$avg_log2mf)
      node_time = c(node_time, tb$node_time)
    }
    lm(formula = node_time ~ log2mf + 0)
  })
  nested_density$nested_density$avg_density = map2(nested_density$nested_density$avg_density,
                                                   nested_density$nested_density$commitment_bias,
                                                   function(density, b) {
                                                     density$table$bin = density$table$bin + log2(b)
                                                     density
                                                   })
  nested_density$nested_density$avg_density = map2(nested_density$nested_density$avg_density,
                                                   nested_density$nested_density$lm,
                                                   function(density, lm) {
                                                     density$table$time = predict(lm, newdata = list(log2mf = density$table$bin))
                                                     density
                                                   })
  nested_density
}

#' NEEDS TO BE UPDATED
#'
#' @param nested_data nested_data object
#' @param node_order_list list of nodes numbered numerically from root(1) to tips
#' @return nested_data object with updated cumulative mutated fraction values
#' @export
convert_mf_to_time_simple <- function(nested_density) {
  nested_density$nested_density$avg_density = map(nested_density$nested_density$avg_density,
                                                  nested_density$nested_density$doubling_time,
                                                   function(density, m) {
                                                     density$table$time = -density$table$bin*m
                                                     density
                                                   })
  nested_density
}

#' NEEDS TO BE UPDATED
#'
#' @param nested_data nested_data object
#' @param node_order_list list of nodes numbered numerically from root(1) to tips
#' @return nested_data object with updated cumulative mutated fraction values
#' @export
align_time_scale <- function(nested_density, time_vec) {
  time_df = as.data.frame(time_vec)
  for (i in 1:length(nested_density$nested_density$avg_density)) {
    table = nested_density$nested_density$avg_density[[i]]$table
    time_df[[nested_density$nested_density$condition[[i]]]] = map_dbl(time_vec, function(t) {
      t1 = max(table$time[table$time <= t])
      t2 = min(table$time[table$time >= t])
      if (t1 == t2) {
        return(table$count[table$time == t])
      }
      if (is.infinite(t1)) {
        t1 = t
        c1 = 0
      } else {
        c1 = table$count[table$time == t1]
      }
      if (is.infinite(t2)) {
        t2 = t
        c2 = 0
      } else {
        c2 = table$count[table$time == t2]
      }
      c1 + (c2-c1)*(t-t1)/(t2-t1)
    })
  }

  #calculating overall average
  all_counts = time_df %>% select(-time_vec)
  time_df$average_nonzero = rowSums(all_counts)/rowSums(all_counts != 0)
  time_df$average = rowMeans(all_counts)
  time_df$average_nonzero[is.na(time_df$average_nonzero)] = 0
  time_df
}

#' NEEDS TO BE UPDATED
#'
#' @param nested_data nested_data object
#' @param node_order_list list of nodes numbered numerically from root(1) to tips
#' @return nested_data object with updated cumulative mutated fraction values
#' @export
align_counts <- function(nested_density) {
  df = as.data.frame(nested_density$nested_density$avg_density[[1]]$table$bin)
  names(df) = 'log2mf'
  for (i in 1:length(nested_density$nested_density$avg_density)) {
    table = nested_density$nested_density$avg_density[[i]]$table
    df[[nested_density$nested_density$condition[[i]]]] = table$count
  }

  #calculating overall average
  all_counts = df %>% select(-log2mf)
  df$average_nonzero = rowSums(all_counts)/rowSums(all_counts != 0)
  df$average = rowMeans(all_counts)
  df$average_nonzero[is.na(df$average_nonzero)] = 0
  df
}

#' NEEDS TO BE UPDATED
#'
#' @param nested_data nested_data object
#' @param node_order_list list of nodes numbered numerically from root(1) to tips
#' @return nested_data object with updated cumulative mutated fraction values
#' @export
load_simulated_data <- function(path, num_samples, mat_name = 'mut_frac_mat') {
  param_tb = tibble(sim = 1:num_samples)
  param_tb$filename = paste0('mf_table_', param_tb$sim, '.rda')
  param_tb[[mat_name]] = map(param_tb$filename, function(fn) {
    if (!file.exists(paste0(path, fn))) {
      return(NULL)
    }
    load(paste0(path, fn))
    mf_table %>% filter(mosaic_fraction > 0) %>% nest(data = -c(sequence)) %>%
      mutate(avg_log2mf = map_dbl(data, function(data) {log2(mean(data$mosaic_fraction))})) %>% unnest(cols = data)
  })
  param_tb = param_tb[!is.null(param_tb[[mat_name]]),]
  param_tb_unnest = unnest(param_tb, cols = mat_name) %>%
    select(sim, tissue = CellType, sequence, mosaic_fraction, avg_log2mf, node_time)
  # param_tb_unnest = unnest(param_tb, cols = mat_name) %>%
  #   select(sim, ID, tissue = CellType, sequence, mosaic_fraction, avg_log2mf, node_time)
  param_tb_unnest$sample = paste0(param_tb_unnest$sim, param_tb_unnest$tissue)
  param_tb_unnest = unique(param_tb_unnest)
  raw_data = structure(list(data = param_tb_unnest), class = 'raw')
  raw_data
}
