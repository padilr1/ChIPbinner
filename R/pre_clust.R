#!/usr/bin/env Rscript
#' Title
#' @title Pre-processing of normalized bigWig files.
#'
#' @description Finds overlapping regions between two bigWig files, removes the bottom and top 1% of bins across the two samples and finally generates a matrix of scores and BED file of coordinates to be used with HDBScan for generating clusters of bins with similar scores.
#'
#' @param out_dir a character string specifying the output directory for the matrix of scores, BED file of genomic coordinates and preliminary density-based scatterplot.
#' @param treated_samp_norm_bw a character string specifying the normalized bigWig file for the treated sample.
#' @param wildtype_samp_norm_bw a character string specifying the normalized bigWig file for the wildtype sample.
#' @param output_filename a character string specifying the file-name for the matrix and pooled BED file.
#' @param are_R_objects a logical indicating whether the inputted bigwig files are R objects. It'll use load() for the reps as opposed to reading them in via rtracklayer::import.bed(). Defaults to FALSE.
#'
#' @return a matrix file of scores, a pooled BED file of genomic coordinates for the bins and a preliminary density-based scatterplot for the two samples being compared.
#' @export
#'
#' @include norm_bw.R
#'
#' @example inst/examples/example_pre_clust.R
pre_clust <- function(out_dir,
                     treated_samp_norm_bw,
                     wildtype_samp_norm_bw,
                     output_filename,
                     are_R_objects=FALSE) {
  # directory parameters
  out_dir <- paste0(out_dir)
  # samples labels
  treated_samp_label <- basename(tools::file_path_sans_ext(treated_samp_norm_bw))
  wildtype_samp_label <- basename(tools::file_path_sans_ext(wildtype_samp_norm_bw))
  # loading function
  loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  # check if inputted files are R objects or not
  if (are_R_objects == "FALSE"){
    # load via import.bed
    treated_samp <- rtracklayer::import.bw(paste0(treated_samp_norm_bw))
    wildtype_samp <- rtracklayer::import.bw(paste0(wildtype_samp_norm_bw))
  } else if (are_R_objects == "TRUE") {
    # otherwise load via loading function
    treated_samp <- loadRData(paste0(treated_samp_norm_bw))
    wildtype_samp <- loadRData(paste0(wildtype_samp_norm_bw))
  }
  d <- list(treated_samp,wildtype_samp)
  names(d) <- c(treated_samp_label,wildtype_samp_label)
  # match genomic regions between the wildtype sample and the treated sample
  d[[wildtype_samp_label]] <- IRanges::subsetByOverlaps(d[[wildtype_samp_label]], d[[treated_samp_label]])
  d[[treated_samp_label]] <- IRanges::subsetByOverlaps(d[[treated_samp_label]], d[[wildtype_samp_label]])
  lapply(d, length)
  r <- lapply(d, function(y) y[y$score != 0]) %>%
    Reduce(function(a, b) a[IRanges::overlapsAny(a, b)], .) %>%
    GenomicRanges::granges()
  # find overlaps; write the matrix with each line corresponding to the score for a given bin; the pooled BED file contains the genomic coordinates
  split(d,c("","")) %>%
    lapply(function(x) {
      o <- lapply(x, function(y) {
        IRanges::findOverlaps(r, y) %>%
          to() %>%
          {
            y[.]
          } %>%
          score()
      }) %>%
        bind_cols() %>%
        `names<-`(c("x", "y"))
      ok <- o$x > quantile(o$x, .01) &
        o$x < quantile(o$x, .99) &
        o$y > quantile(o$y, .01) &
        o$y < quantile(o$y, .99)
      readr::write_csv(o[ok, ], sprintf("%s/%s_mat.csv", out_dir, output_filename), col_names = F)
      rtracklayer::export.bed(r[ok], sprintf("%s/%s_pooled.bed", out_dir, output_filename))
    })
  # read in matrix
  mat <- data.table::fread(sprintf("%s/%s_mat.csv", out_dir, output_filename))
  # open jpeg file
  jpeg(filename = sprintf("%s/%s_smoothScatter.jpeg", out_dir, output_filename))
  # create the plot
  smoothScatter(y = mat$V1, x = mat$V2, xlab = wildtype_samp_label, ylab = treated_samp_label)
  # close the file
  dev.off()
  # print output message
  print("bigWig files pre-processed. Output files ready for input into clustering algorithm!")
}
