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
  # load genic and intergenic regions
  gene <- system.file("extdata", "hg38_gene.bed.gz", package = "ChIPbinner")
  igr <- system.file("extdata", "hg38_intergenic.bed.gz", package = "ChIPbinner")
  # generate density-based scatterplot
  density_based_scatterplot(path_for_norm_bw = testthat::test_path("testdata"),
    out_dir = testthat::test_path("testdata"),
    treated_samp_label = "NSD1_KO",
    baseline_samp_label = "WT",
    cell_line = "Cal27",
    histone_mark = "H3K36me2",
    annotated_clusters = testthat::test_path("testdata","cons.Cal27.WT.NSD1_KO.H3K36me2.10kb.rda"),
    number_of_clusters = 3,
    gene = gene,
    intergenic = igr,
    title_of_plot = "H3K36me2",
    window_size = 10,
    pow = 1.1,
    show_legend = TRUE,
    min=-5,
    max=2,
    bin_size=50,
    show_scales = FALSE,
    xaxis_label="WT",
    yaxis_label="NSD1_KO",
    height_of_figure=6,
    width_of_figure=15)
  expect_snapshot_output(x = "Density-based scatterplots generated!")
})
