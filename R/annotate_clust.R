#!/usr/bin/env Rscript
#' Title
#' @title Annotate clusters.
#' @description Annotates clusters identified after using HDBSCAN.
#' @param number_of_clusters an integer specifying the total number of clusters identified by the HDBSCAN algorithm.
#' @param matrix_file a character string specifying the matrix file of enrichment scores for the two samples being compared, generated from the pre_clus() function.
#' @param pooled_bed_file a character string specifying the pooled BED file generated from 'pre_clus()' consisting of genomic coordinates found in both samples being compared.
#' @param hdbscan_output_file a character string specifying the HDBSCAN output (.txt) file which specifies whether each bin falls into a specific cluster or not.
#' @param output_filename a character string specifying the output file name.
#' @param out_dir a character string specifying the output directory.
#'
#' @return a R object consisting of annotated clusters and their respective set of bins.
#' @export
#'
#' @include generate_clust.R
#'
#' @example inst/examples/example_annotate_clust.R
annotate_clust <- function(number_of_clusters,
                           matrix_file,
                           pooled_bed_file,
                           hdbscan_output_file,
                           output_filename,
                           out_dir) {
  suppressWarnings({
  # out dir
  out_dir <- paste0(out_dir)
  # number of cluster
  number_of_clusters <- as.integer(paste0(number_of_clusters))
  # read in matrix file of enrichment scores
  mat <- readr::read_csv(matrix_file, col_names = F)
  # read in pooled BED file of genomic coordinates
  pooled_BED <- rtracklayer::import.bed(pooled_bed_file)
  # read in HDBSCAN output file
  mat$clu <- readr::read_csv(file = hdbscan_output_file, col_names = F)$X1
  # separate workflows depending on the number of clusters
  if (number_of_clusters > 1 & number_of_clusters < 25) {
    ctr <- mat %>%
      dplyr::filter(clu != -1) %>%
      dplyr::group_by(clu) %>%
      dplyr::summarise_all(mean) %>%
      dplyr::mutate(mu = 0.5 * (X1 + X2)) %>%
      dplyr::arrange(mu)
    results_list <- list()
    for (i in 1:number_of_clusters) {
      results_list[[paste0(toupper(letters[i]))]] <- list(pooled_BED[mat$clu == ctr$clu[i]])
    }
    cons <- results_list %>%
      lapply(function(x) {
        Reduce(function(a, b) {
          a[IRanges::overlapsAny(a, b)]
        }, x) %>%
          GenomicRanges::granges()
      })
  } else {
    print("Too few or too many clusters indicated. Minimum number of clusters is 2 and maximum number is 25.")
    #break
  }
  save(cons, file = (sprintf("%s/%s.annotated_clusters.rda", out_dir, output_filename)))
  print("Clusters annotated!")
  })
}
