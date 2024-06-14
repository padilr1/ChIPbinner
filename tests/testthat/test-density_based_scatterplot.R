test_that("Generating density-based scatterplots work", {
  # load libraries
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
  library(lattice)
  library(gridExtra)
  library(DiffBind)
  library(patchwork)
  library(sf)
  library(MASS)
  library(lwgeom)
  library(ggrepel)
  library(hexbin)
  library(ggrastr)
  library(viridis)
  library(pals)
  library(isoband)
  # load genic and intergenic regions
  gene <- system.file("extdata", "hg38_gene.bed", package = "ChIPbinner")
  igr <- system.file("extdata", "hg38_intergenic.bed", package = "ChIPbinner")
  # generate density-based scatterplot
  density_based_scatterplot(path_for_norm_bw = testthat::test_path("testdata"),
    out_dir = testthat::test_path("testdata"),
    treated_samp_label = "NSD1_KO",
    baseline_samp_label = "WT",
    cell_line = "Cal27",
    histone_mark = "H3K36me2",
    # matrix_file = "Cal27.WT.NSD1_KO.H3K36me2.10kb.mat.csv"
    matrix_file = testthat::test_path("testdata","Cal27.WT.NSD1_KO.H3K36me2.10kb.mat.csv"),
    # pooled_bed_file = "Cal27.WT.NSD1_KO.H3K36me2.10kb.pooled.bed"
    pooled_bed_file = testthat::test_path("testdata","Cal27.WT.NSD1_KO.H3K36me2.10kb.pooled.bed"),
    # hdbscan_output_file ="clus.Cal27.WT.NSD1_KO.H3K36me2.5000.5000.txt"
    hdbscan_output_file = testthat::test_path("testdata","clus.Cal27.WT.NSD1_KO.H3K36me2.10kb.5000.5000.txt"),
    number_of_clusters = 3,
    gene = gene,
    intergenic = igr,
    title_of_plot = "H3K36me2",
    window_size = 10,
    pow = 1.25,
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
