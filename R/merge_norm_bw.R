#!/usr/bin/env Rscript
#' Title
#' @title Merge normalized bigWig files.
#'
#' @description Merge normalized bigWig files, which are normally replicates. This function should be used after implementing the 'norm_bw()' function.
#'
#' @param rep1 Replicate 1 to be merged. can be a bigWig file or an R object.
#' @param rep2 Replicate 2 to be merged. can be a bigWig file or an R object.
#' @param are_R_objects Boolean term (true or false) to indicate if the inputted replicates are R objects. It'll use load() for the reps as opposed to reading them in via rtracklayer::import.bed(). Default to FALSE.
#' @param merged_samples_label Sample label of the merged bigWig file.
#' @param out_dir Output directory for merged bigWig files.
#'
#' @return a single merged bigWig file.
#' @export
#'
#' @example
merge_norm_bw <- function(rep1,
                          rep2,
                          are_R_objects=FALSE,
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
  rep1 <- rtracklayer::import.bed(paste0(rep1))
  rep2 <- rtracklayer::import.bed(paste0(rep2))
  } else if (are_R_objects == "TRUE") {
    # otherwise load via loading function
    rep1 <- loadRData(rep1)
    rep2 <- loadRData(rep1)
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
  l <- list(rep1$score, rep2$score)
  # Calculate the mean for each row
  row_means <- rowMeans(do.call(cbind, l))
  # Assign the mean values to the 'score' column in 'merged_bw'; the score for the merged bigWig file is now the mean of the the two replicates
  merged_bw$score <- row_means
  # export bw
  rtracklayer::export.bw(merged_bw, con = sprintf("%s/%s.norm.bw", out_dir,merged_samples_label))
  # save bw as R object
  save(merged_bw, file = sprintf("%s/%s.rda", out_dir,merged_samples_label))
}
