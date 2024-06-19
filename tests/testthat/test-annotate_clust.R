test_that("annotating clusters works", {
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
  # library(lattice)
  # library(gridExtra)
  # library(DiffBind)
  # library(idr2d)
  # library(patchwork)
  annotate_clust(histone_mark = "H3K36me2",
                 window_size=10,
                 cell_line = "Cal27",
                 baseline_samp_label = "WT",
                 treated_samp_label = "NSD1_KO",
                 number_of_clusters = 3,
                 matrix_file = testthat::test_path("testdata","Cal27.WT.NSD1_KO.H3K36me2.10kb.mat.csv"),
                 pooled_bed_file = testthat::test_path("testdata","Cal27.WT.NSD1_KO.H3K36me2.10kb.pooled.bed"),
                 hdbscan_output_file= testthat::test_path("testdata","clus.Cal27.WT.NSD1_KO.H3K36me2.10kb.5000.5000.gz"),
                 out_dir=testthat::test_path("testdata"))
  expect_snapshot_output(x = "Clusters annotated!")
})
