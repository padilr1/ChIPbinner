#!/usr/bin/env Rscript
#' Title
#' @title Merging of normalized bigWig files.
#'
#' @description Performs merging of normalized bigWig files, which are normally replicates.
#'
#' @param are_R_objects a logical indicating whether the inputs are bigWig files or R objects.
#' @param rep1 a character string specifying the first bigWig file or R object to be merged. Usually the first replicate.
#' @param rep2 a character string specifying the second bigWig file or R object to be merged. Usually the second replicate.
#' @param merged_samples_label a character string specifying the filename of the merged bigWig file.
#' @param out_dir a character string specifying the output directory for the merged bigWig file.
#'
#' @return a single merged bigWig file.
#' @export
#'
#' @example inst/examples/example_merged_norm_bw.R
merge_norm_bw <- function(are_R_objects=FALSE,
                          rep1,
                          rep2,
                          merged_samples_label,
                          out_dir) {
  # example = inst/examples/example_merged_norm_bw.R
  # output dir
  out_dir <- paste0(out_dir)
  # loading function
  loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  # check if inputted files are R objects or not
  if (are_R_objects == "FALSE"){
  # load via import.bed
  rep1 <- rtracklayer::import.bw(paste0(rep1))
  rep2 <- rtracklayer::import.bw(paste0(rep2))
  } else if (are_R_objects == "TRUE") {
    # otherwise load via loading function
    rep1 <- loadRData(paste0(rep1))
    rep2 <- loadRData(paste0(rep2))
  }
  # merged samples label
  merged_samples_label <- paste0(merged_samples_label)
  # subset regions in both reps
  merged_bw <- IRanges::subsetByOverlaps(rep1,rep2)
  # make score NULL for merged bigwig
  merged_bw$score <- NULL
  # ensure only overlapping regions are included between the two reps
  rep1 <- IRanges::subsetByOverlaps(rep1,rep2)
  rep2 <-IRanges::subsetByOverlaps(rep2,rep1)
  # list scores
  l <- list(as.numeric(rep1$score), as.numeric(rep2$score))
  # Calculate the mean for each row
  row_means <- rowMeans(do.call(cbind, l))
  # Assign the mean values to the 'score' column in 'merged_bw'; the score for the merged bigWig file is now the mean of the the two replicates
  merged_bw$score <- row_means
  # export bw
  rtracklayer::export.bw(merged_bw, con = sprintf("%s/%s.bw", out_dir,merged_samples_label))
  # save bw as R object
  save(merged_bw, file = sprintf("%s/%s.rda", out_dir,merged_samples_label))
}
