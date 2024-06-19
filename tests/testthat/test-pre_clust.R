test_that("processing bigWig files prior to clustering works", {
  # load libraries
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
  # library(lattice)
  # library(gridExtra)
  # library(DiffBind)
  # library(idr2d)
  # library(patchwork)
  # pre-process bigWig files
  pre_clust(path_for_norm_bw = testthat::test_path("testdata"),
                       out_dir = testthat::test_path("testdata"),
                       window_size = 10,
                       treated_samp_label = "NSD1_KO",
                       baseline_samp_label = "WT",
                       cell_line = "Cal27",
                       histone_mark = "H3K36me2")
  expect_snapshot_output(x = "bigWig files pre-processed. Output files ready for input into clustering algorithm!")
})
