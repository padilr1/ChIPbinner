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
density_based_scatterplot(path_for_norm_bw = system.file("extdata/norm_bw", package = "ChIPbinner"),
                          out_dir = system.file("extdata", package = "ChIPbinner"),
                          treated_samp_label = "NSD1_KO",
                          baseline_samp_label = "WT",
                          cell_line = "Cal27",
                          histone_mark = "H3K36me2",
                          matrix_file = system.file("extdata","Cal27.WT.NSD1_KO.H3K36me2.10kb.mat.csv", package = "ChIPbinner"),
                          pooled_bed_file = system.file("extdata","Cal27.WT.NSD1_KO.H3K36me2.10kb.pooled.bed", package = "ChIPbinner"),
                          hdbscan_output_file = system.file("extdata","clus.Cal27.WT.NSD1_KO.H3K36me2.10kb.5000.5000.txt", package = "ChIPbinner"),
                          number_of_clusters = 3,
                          gene = system.file("extdata", "hg38_gene.bed", package = "ChIPbinner"),
                          intergenic = system.file("extdata", "hg38_intergenic.bed", package = "ChIPbinner"),
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
