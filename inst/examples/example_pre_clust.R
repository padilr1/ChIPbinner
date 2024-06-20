library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(lattice)
library(gridExtra)
library(DiffBind)
library(idr2d)
library(patchwork)
pre_clust(
  out_dir = system.file("extdata",package = "ChIPbinner"),
  treated_samp_norm_bw = system.file("extdata/data","downsampled.Cal27.WT.H3K36me2.10kb.norm_bw.rda", package = "ChIPbinner"),
  wildtype_samp_norm_bw = system.file("extdata/data","downsampled.Cal27.NSD1_KO.H3K36me2.10kb.norm_bw.rda", package = "ChIPbinner"),
  output_filename = "downsampled.Cal27.WT.NSD1_KO.H3K36me2.10kb",
  are_R_objects = TRUE
)
