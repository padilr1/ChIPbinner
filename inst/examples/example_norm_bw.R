library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
norm_bw(out_dir = system.file("extdata",package = "ChIPbinner"),
        genome_assembly = "hg38",
        chip_samp_binned_file = system.file("extdata", "downsampled.Cal27.WT.H3K36me2.10kb.bed.gz", package = "ChIPbinner"),
        chip_samp_library_size = 64093770,
        use_control = TRUE,
        control_binned_file = system.file("extdata", "downsampled.Cal27.WT_input.H3K36me2.10kb.bed.gz", package = "ChIPbinner"),
        control_library_size = 52047022,
        raw_count_cutoff = 0,
        pseudocount = 1e-15,
        scaling_factor = 0.450328805)
