# library(tidyverse)
# library(rtracklayer)
# library(GenomicRanges)
# library(data.table)
# library(lattice)
# library(gridExtra)
# library(DiffBind)
# library(patchwork)
# load chrom sizes for hg38 assembly
# hg38_chrom_sizes <- system.file("extdata", "hg38_chrom.sizes", package = "ChIPbinner")
# # load blacklisted regions
# blacklisted_regions <- system.file("extdata", "hg38_blacklist.bed", package = "ChIPbinner")
# # load WT sample and corresponding input
# WT <- system.file("extdata", "Cal27.WT.H3K36me2.10kb.bed", package = "ChIPbinner")
# WT_input <- system.file("extdata", "Cal27.WT_input.H3K36me2.10kb.bed", package = "ChIPbinner")
# # load NSD1_KO sample and corresponding input
# NSD1_KO <- system.file("extdata", "Cal27.NSD1_KO.H3K36me2.10kb.bed", package = "ChIPbinner")
# NSD1_KO_input <- system.file("extdata", "Cal27.NSD1_KO_input.H3K36me2.10kb.bed", package = "ChIPbinner")
# # load genic and intergenic regions
# gene <- system.file("extdata", "hg38_gene.bed", package = "ChIPbinner")
# igr <- system.file("extdata", "hg38_intergenic.bed", package = "ChIPbinner")
# generate normalized bigwig for WT sample
norm_bw(out_dir = system.file("extdata/norm_bw",package = "ChIPbinner"),
        chromSizes = system.file("extdata", "hg38_chrom.sizes", package = "ChIPbinner"),
        blacklist = system.file("extdata", "hg38_blacklist.bed", package = "ChIPbinner"),
        cell_line = "Cal27",
        histone_mark = "H3K36me2",
        treated_samp_label = "WT",
        treated_samp_file= system.file("extdata", "Cal27.WT.H3K36me2.10kb.bed", package = "ChIPbinner"),
        treated_samp_library_size = 64093770,
        use_control = TRUE,
        control_label = "WT_input",
        control_file = system.file("extdata", "Cal27.WT_input.H3K36me2.10kb.bed", package = "ChIPbinner"),
        control_library_size = 52047022,
        raw_count_cutoff = 0,
        window_size = 10,
        pseudocount = 1e-15,
        scaling_factor = 0.450328805)
