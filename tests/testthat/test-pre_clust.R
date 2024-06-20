test_that("processing bigWig files prior to clustering works", {
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
  library(ChIPbinner)
  # downsampled norm bigwigs
  pre_clust(
    out_dir = testthat::test_path("testdata"),
    treated_samp_norm_bw = system.file("extdata/data","downsampled.Cal27.NSD1_KO.H3K36me2.10kb.norm_bw.rda", package = "ChIPbinner"),
    wildtype_samp_norm_bw = system.file("extdata/data","downsampled.Cal27.WT.H3K36me2.10kb.norm_bw.rda", package = "ChIPbinner"),
    output_filename = "downsampled.Cal27.WT.NSD1_KO.H3K36me2.10kb",
    are_R_objects = TRUE
  )
  # complete norm bigwigs
  # pre_clust(
  #   out_dir = testthat::test_path("testdata"),
  #   treated_samp_norm_bw = testthat::test_path("testdata","Cal27_NSD1_KO_H3K36me2_10kb_norm_bw.rda"),
  #   wildtype_samp_norm_bw = testthat::test_path("testdata","Cal27_WT_H3K36me2_10kb_norm_bw.rda"),
  #   output_filename = "complete_Cal27.WT.NSD1_KO.H3K36me2.10kb",
  #   are_R_objects = TRUE
  # )
  expect_snapshot_output(x = "bigWig files pre-processed. Output files ready for input into clustering algorithm!")
})
