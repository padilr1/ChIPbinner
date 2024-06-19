library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
# library(lattice)
# library(gridExtra)
# library(DiffBind)
# library(patchwork)
norm_bw(out_dir = system.file("extdata/norm_bw",package = "ChIPbinner"),
        chromSizes = system.file("extdata", "hg38_chrom.sizes.gz", package = "ChIPbinner"),
        blacklist = system.file("extdata", "hg38_blacklist.bed.gz", package = "ChIPbinner"),
        cell_line = "Cal27",
        histone_mark = "H3K36me2",
        treated_samp_label = "WT",
        treated_samp_file= system.file("extdata", "Cal27.WT.H3K36me2.10kb.bed.gz", package = "ChIPbinner"),
        treated_samp_library_size = 64093770,
        use_control = TRUE,
        control_label = "WT_input",
        control_file = system.file("extdata", "Cal27.WT_input.H3K36me2.10kb.bed.gz", package = "ChIPbinner"),
        control_library_size = 52047022,
        raw_count_cutoff = 0,
        window_size = 10,
        pseudocount = 1e-15,
        scaling_factor = 0.450328805)
