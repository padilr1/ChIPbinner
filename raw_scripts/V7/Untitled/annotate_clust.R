#!/usr/bin/env Rscript
#' Title
#' @title Annotate clusters.
#' @description Annotates clusters identified after using HDBSCAN.
#' @param histone_mark The broad histone mark used for the analysis.
#' @param window_size The window size in kilobytes. For example, use 10 for samples binned in 10 kilobyte (kb) bins.
#' @param cell_line The cell line of the samples.
#' @param baseline_samp_label The baseline sample name/label to which to compare the treated sample. This is usually the wildtype (WT) sample.
#' @param treated_samp_label The treated sample name/label.
#' @param number_of_clusters The total number of clusters identified by the HDBSCAN algorithm.
#' @param matrix_file The matrix file of enrichment scores for the two samples being compared, generated from the pre_clus() function.
#' @param pooled_bed_file The pooled BED file generated from 'pre_clus()' consisting of genomic coordinates found in both samples being compared.
#' @param hdbscan_output_file The HDBSCAN output which specifies whether each bin falls into a specific cluster or not.
#' @param out_dir The output directory.
#'
#' @return a R object consisting of multiple annotated clusters, each with their set of bins.
#' @export
#'
#' @example inst/examples/example_annotate_clust.R
annotate_clust <- function(histone_mark,
                           window_size,
                           cell_line,
                           baseline_samp_label,
                           treated_samp_label,
                           number_of_clusters,
                           matrix_file,
                           pooled_bed_file,
                           hdbscan_output_file,
                           out_dir) {
  # out dir
  out_dir <- paste0(out_dir)
  # samples info
  treated_samp_label <- paste0(treated_samp_label)
  baseline_samp_label <- paste0(baseline_samp_label)
  cell_line <- paste0(cell_line)
  mark <- paste0(histone_mark)
  # window size of the bins
  window_size <- paste0(".", window_size, "kb.")
  # number of cluster
  number_of_clusters <- as.integer(paste0(number_of_clusters))
  # read in matrix file of enrichment scores
  mat <- readr::read_csv(matrix_file, col_names = F)
  # pooled BED file of genomic coordinates
  pooled_BED <- rtracklayer::import.bed(pooled_bed_file)
  # read in HDBSCAN output file
  mat$clu <- readr::read_csv(file = hdbscan_output_file, col_names = F)$X1
  # separate workflows depending on the number of clusters
  if (number_of_clusters > 1 & number_of_clusters < 25) {
    ctr <- mat %>%
      dplyr::filter(clu != -1) %>%
      group_by(clu) %>%
      summarise_all(mean) %>%
      mutate(mu = 0.5 * (X1 + X2)) %>%
      arrange(mu)
    results_list <- list()
    for (i in 1:number_of_clusters) {
      results_list[[paste0(toupper(letters[i]))]] <- list(pooled_BED[mat$clu == ctr$clu[i]])
    }
    cons <- results_list %>%
      lapply(function(x) {
        Reduce(function(a, b) {
          a[overlapsAny(a, b)]
        }, x) %>%
          granges()
      })
  } else {
    print("Too few or too many clusters indicated. Minimum number of clusters is 2 and maximum number is 25.")
    break
  }
  save(cons, file = (sprintf("%s/cons.%s.%s.%s.%s%srda", out_dir, cell_line, baseline_samp_label, treated_samp_label, mark, window_size)))
  print("Clusters annotated!")
}
