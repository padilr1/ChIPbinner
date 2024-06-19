test_that("generating genic/intergenic scatterplots works", {
  # load libraries
  library(data.table)
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggrepel)
  library(hexbin)
  library(viridis)
  library(lwgeom)
  library(patchwork)
  library(matrixStats)
  library(isoband)
  library(sf)
  library(MASS)
  library(pals)
  # load genic and intergenic regions
  gene <- system.file("extdata", "hg38_gene.bed.gz", package = "ChIPbinner")
  igr <- system.file("extdata", "hg38_intergenic.bed.gz", package = "ChIPbinner")
  # generate genic/intergenic scatterplot
  genic_intergenic_scatterplot(path_for_norm_bw = testthat::test_path("testdata"),
                               out_dir = testthat::test_path("testdata"),
                               gene = gene,
                               intergenic = igr,
                               cell_line = "Cal27",
                               baseline_samp_label = "WT",
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
                               show_scales = F,
                               show_legend = T,
                               legend_pos="left")
  expect_snapshot_output(x = "Generated genic/intergenic scatterplot!")
})
