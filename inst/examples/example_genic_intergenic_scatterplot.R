library(data.table)
library(tidyverse)
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
library(GenomicRanges)
library(rtracklayer)
genic_intergenic_scatterplot(path_for_norm_bw = system.file("extdata/norm_bw", package = "ChIPbinner"),
                             out_dir = system.file("extdata", package = "ChIPbinner"),
                             gene = system.file("extdata", "hg38_gene.bed", package = "ChIPbinner"),
                             intergenic = system.file("extdata", "hg38_intergenic.bed", package = "ChIPbinner"),
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
                             show_scales = FALSE,
                             show_legend = TRUE,
                             legend_pos="left")
