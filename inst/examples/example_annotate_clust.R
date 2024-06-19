#!/usr/bin/env Rscript
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(lattice)
library(gridExtra)
library(DiffBind)
library(idr2d)
library(patchwork)
annotate_clust(histone_mark = "H3K36me2",
               window_size=10,
               cell_line = "Cal27",
               baseline_samp_label = "WT",
               treated_samp_label = "NSD1_KO",
               number_of_clusters = 3,
               matrix_file = system.file("extdata","Cal27.WT.NSD1_KO.H3K36me2.10kb.mat.csv", package = "ChIPbinner"),
               pooled_bed_file = system.file("extdata","Cal27.WT.NSD1_KO.H3K36me2.10kb.pooled.bed", package = "ChIPbinner"),
               hdbscan_output_file= system.file("extdata","clus.Cal27.WT.NSD1_KO.H3K36me2.10kb.5000.5000.txt.gz", package = "ChIPbinner"),
               out_dir = system.file("extdata", package = "ChIPbinner"))

