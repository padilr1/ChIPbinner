test_that("Generating density-based scatterplots work", {
  # load libraries
  ## required libs
  library(data.table)
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggrepel)
  library(hexbin)
  library(viridis)
  library(lwgeom)
  library(patchwork)
  library(pals)
  library(MASS)
  library(lattice)
  library(gridExtra)
  library(DiffBind)
  library(sf)
  library(lwgeom)
  library(ggrastr)
  library(isoband)
  # generate density-based scatterplot
  ## downsampled
  # density_based_scatterplot(
  #   out_dir = testthat::test_path("testdata"),
  #   genome_assembly = "hg38",
  #   are_R_objects = TRUE,
  #   output_filename = "downsampled",
  #   wildtype_samp_norm_bw = testthat::test_path("testdata","downsampled.Cal27.WT.H3K36me2.10kb.bed.norm_bw.rda"),
  #   treated_samp_norm_bw =  testthat::test_path("testdata","downsampled.Cal27.NSD1_KO.H3K36me2.10kb.bed.norm_bw.rda"),
  #   cell_line = "Cal27",
  #   histone_mark = "H3K36me2",
  #   annotated_clusters = testthat::test_path("testdata", "downsampled.annotated_clusters.rda"),
  #   number_of_clusters = 3,
  #   title_of_plot = "H3K36me2",
  #   pow = 1.1,
  #   show_legend = TRUE,
  #   min = -5,
  #   max = 2,
  #   bin_size = 50,
  #   show_scales = FALSE,
  #   xaxis_label = "WT",
  #   yaxis_label = "NSD1_KO",
  #   height_of_figure = 6,
  #   width_of_figure = 15
  # )
  ## complete norm bw
  density_based_scatterplot(
    out_dir = testthat::test_path("testdata"),
    genome_assembly = "hg38",
    are_R_objects = TRUE,
    output_filename = "complete",
    wildtype_samp_norm_bw = system.file("extdata/data","Cal27_WT_H3K36me2_10kb_norm_bw.rda", package = "ChIPbinner"),
    treated_samp_norm_bw = system.file("extdata/data","Cal27_NSD1_KO_H3K36me2_10kb_norm_bw.rda", package = "ChIPbinner"),
    cell_line = "Cal27",
    histone_mark = "H3K36me2",
    annotated_clusters = system.file("extdata/data","complete_annotated_clusters.rda", package = "ChIPbinner"),
    number_of_clusters = 3,
    title_of_plot = "H3K36me2",
    pow = 1.1,
    show_legend = TRUE,
    min = -5,
    max = 2,
    bin_size = 50,
    show_scales = FALSE,
    xaxis_label = "WT",
    yaxis_label = "NSD1_KO",
    height_of_figure = 6,
    width_of_figure = 15
  )
  expect_snapshot_output(x = "Density-based scatterplots generated!")
})
