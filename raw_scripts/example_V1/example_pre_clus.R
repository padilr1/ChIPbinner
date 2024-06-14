# load libraries
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(lattice)
library(gridExtra)
library(DiffBind)
library(idr2d)
library(patchwork)
# pre-process bigWig files
pre_clus(path_for_norm_bw = system.file("extdata/norm_bw", package = "ChIPbinner"),
         out_dir = system.file("extdata", package = "ChIPbinner"),
         window_size = 10,
         treated_samp_label = "NSD1_KO",
         baseline_samp_label = "WT",
         cell_line = "Cal27",
         histone_mark = "H3K36me2")
