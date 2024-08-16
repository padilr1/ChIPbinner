test_that("generating genic/intergenic scatterplots works", {
  # load libraries
  # library(data.table)
  # library(tidyverse)
  # library(GenomicRanges)
  # library(rtracklayer)
  # library(ggrepel)
  # library(hexbin)
  # library(viridis)
  # library(lwgeom)
  # library(patchwork)
  # library(matrixStats)
  # library(isoband)
  # library(sf)
  # library(MASS)
  # library(pals)
  # generate genic/intergenic scatterplot
  ## sampled norm bw
  # testthat::test_path("testdata","downsampled.Cal27.WT.H3K36me2.10kb.bed.norm_bw.rda")
  # testthat::test_path("testdata","downsampled.Cal27.NSD1_KO.H3K36me2.10kb.bed.norm_bw.rda")
  # system.file("data/downsampled.Cal27.WT.H3K36me2.10kb.norm_bw.rda", package = "ChIPbinner")
  genic_intergenic_scatterplot(
    out_dir = testthat::test_path("testdata"),
    genome_assembly = "hg38",
    cell_line = "Cal27",
    wildtype_samp_norm_bw = system.file("extdata/data","downsampled.Cal27.WT.H3K36me2.10kb.norm_bw.rda", package = "ChIPbinner"),
    treated_samp_norm_bw = system.file("extdata/data","downsampled.Cal27.NSD1_KO.H3K36me2.10kb.norm_bw.rda", package = "ChIPbinner"),
    are_R_objects = TRUE,
    histone_mark = "H3K36me2",
    output_filename = "downsampled_Cal27.WT_NSD1KO.H3K36me2.10kb",
    title_of_plot = "Cal27 ChIP-seq H3K36me2",
    xaxis_label = "WT",
    yaxis_label = "NSD1_KO",
    max_x = 1,
    max_y = 1,
    min_x = -5,
    min_y = -5,
    pow = 1.25,
    show_scales = FALSE,
    show_legend = TRUE,
    legend_pos = "left"
  )
  ## complete norm bw
  # genic_intergenic_scatterplot(out_dir = testthat::test_path("testdata"),
  #                              genome_assembly = "hg38",
  #                              cell_line = "Cal27",
  #                              wildtype_samp_norm_bw = testthat::test_path("testdata","Cal27_WT_H3K36me2_10kb_norm_bw.rda"),
  #                              treated_samp_norm_bw = testthat::test_path("testdata","Cal27_NSD1_KO_H3K36me2_10kb_norm_bw.rda"),
  #                              are_R_objects = TRUE,
  #                              histone_mark = "H3K36me2",
  #                              output_filename = "Cal27.WT_NSD1KO.H3K36me2.10kb",
  #                              title_of_plot = "Cal27 ChIP-seq H3K36me2",
  #                              xaxis_label = "WT",
  #                              yaxis_label = "NSD1_KO",
  #                              max_x = 1,
  #                              max_y = 1,
  #                              min_x = -5,
  #                              min_y = -5,
  #                              pow=1.25,
  #                              show_scales = F,
  #                              show_legend = T,
  #                              legend_pos="left")
  expect_snapshot_output(x = "Generated genic/intergenic scatterplot!")
})
