#!/usr/bin/env Rscript
annotate_clust <- function(number_of_clusters,
                           matrix_file,
                           pooled_bed_file,
                           hdbscan_output_file){
  # number of cluster
  number_of_clusters <- as.integer(paste0(number_of_clusters))
  # read in matrix file of enrichment scores
  mat <- readr::read_csv(matrix_file,col_names = F)
  # pooled BED file of genomic coordinates
  pooled_BED <- rtracklayer::import.bed(pooled_bed_file)
  # read in HDBSCAN output file
  mat$clu <- readr::read_csv(file = hdbscan_output_file, col_names = F)$X1
  # separate workflows depending on the number of clusters
  if (number_of_clusters == 2) {
    ctr <- mat %>%
      dplyr::filter(clu != -1) %>%
      group_by(clu) %>%
      summarise_all(mean) %>%
      mutate(mu = 0.5 * (X1 + X2)) %>%
      arrange(mu)

    a1 <- pooled_BED[mat$clu == ctr$clu[1]]
    b1 <- pooled_BED[mat$clu == ctr$clu[2]]

    cons <- list(
      A = list(a1),
      B = list(b1)
    ) %>%
      lapply(function(x) {
        Reduce(function(a, b) {
          a[overlapsAny(a, b)]
        }, x) %>%
          granges()
      })
  } else if (number_of_clusters == 3) {
    ctr <- mat %>%
      dplyr::filter(clu != -1) %>%
      group_by(clu) %>%
      summarise_all(mean) %>%
      mutate(mu = 0.5 * (X1 + X2)) %>%
      arrange(mu)

    a1 <- pooled_BED[mat$clu == ctr$clu[1]]
    b1 <- pooled_BED[mat$clu == ctr$clu[2]]
    c1 <- pooled_BED[mat$clu == ctr$clu[3]]

    cons <- list(
      A = list(a1),
      B = list(b1),
      C = list(c1)
    ) %>%
      lapply(function(x) {
        Reduce(function(a, b) {
          a[overlapsAny(a, b)]
        }, x) %>%
          granges()
      })
  } else if (number_of_clusters > 3 & number_of_clusters < 25) {
    ctr <- mat %>%
      dplyr::filter(clu != -1) %>%
      group_by(clu) %>%
      summarise_all(mean) %>%
      mutate(mu = 0.5 * (X1 + X2)) %>%
      arrange(mu)
    pre_cons <- list()
    cons <- list()
    for (i in 1:number_of_clusters){
     pre_cons[[paste0(letters[i],1)]] <- pooled_BED[mat$clu == ctr$clu[i]]
     cons[[toupper(letters[i])]] <- list(pre_cons[[paste0(letters[i])]])
    }
    cons <- cons %>%
      lapply(function(x) {
        Reduce(function(a, b) {
          a[overlapsAny(a, b)]
        }, x) %>%
          granges()
      })
  }
  save(cons, file = (sprintf('data/cons/cons.%s.%s.%s.%s%srda',cell_line,control,test,mark,window_size)))
}
