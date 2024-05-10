test_that("generating genic/intergenic scatterplots works", {
  # load chrom sizes for hg38 assembly
  hg38_chrom_sizes <- system.file("extdata", "hg38_chrom.sizes", package = "ChIPbinner")
  # load blacklisted regions
  blacklisted_regions <- system.file("extdata", "hg38_blacklist.bed", package = "ChIPbinner")
  # load WT sample and corresponding input
  WT <- system.file("extdata", "Cal27.WT.H3K36me2.10kb.bed", package = "ChIPbinner")
  WT_input <- system.file("extdata", "Cal27.WT_input.H3K36me2.10kb.bed", package = "ChIPbinner")
  # load NSD1_KO sample and corresponding input
  NSD1_KO <- system.file("extdata", "Cal27.NSD1_KO.H3K36me2.10kb.bed", package = "ChIPbinner")
  NSD1_KO_input <- system.file("extdata", "Cal27.NSD1_KO_input.H3K36me2.10kb.bed", package = "ChIPbinner")
  # load genic and intergenic regions
  gene <- system.file("extdata", "hg38_gene.bed", package = "ChIPbinner")
  igr <- system.file("extdata", "hg38_intergenic.bed", package = "ChIPbinner")
  # load libraries
  library(data.table)
  library(tidyverse)
  library(rtracklayer)
  library(matrixStats)
  library(isoband)
  library(sf)
  library(ggrepel)
  library(viridis)
  library(MASS)
  library(lwgeom)
  library(hexbin)
  library(pals)
  library(patchwork)
  library(gdata)
  # generate genic/intergenic scatterplot
  genic_intergenic_scatterplot(path_for_norm_bw = testthat::test_path("testdata"),
                               out_dir = testthat::test_path("testdata"),
                               gene = gene,
                               intergenic = igr,
                               cell_line = "Cal27",
                               control_label = "WT",
                               treated_samp_label = "NSD1_KO",
                               histone_mark = "H3K36me2",
                               window_size = 10,
                               title_of_plot = "Cal27 ChIP-seq H3K36me2",
                               xaxis_label = "WT",
                               yaxis_label = "NSD1_KO",
                               max_x = 1,
                               max_y = 1,
                               min_x = -5,
                               min_y = -5,
                               pow=1.25,
                               show_scales = F)
  expect_snapshot_output(x = "Generated genic/intergenic scatterplot!")
})
