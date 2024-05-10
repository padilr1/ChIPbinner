library(ChIPbinner)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(lattice)
library(gridExtra)
library(DiffBind)
library(patchwork)
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
# # params
# chrom.sizes=paste0(hg38_chrom_sizes)
# # binned_files_path <- paste0(binned_files_path)
# addition <- 1e-15
# raw_count_cutoff <- 0
# # sample info
# cell_line="Cal27"
# mark="H3K36me2"
# treated_samp_label=paste0("WT") #watch spaces in names
# treated_samp_file=import.bed(paste0(WT))
# #### reference files ####
# # blacklist
# # blacklist=paste0(blacklisted)
# bl <- import.bed(blacklisted_regions)
# # window size
# window_size=paste0(".","10","kb.")
# # scaling factor
# if(is.null(scaling_factor)){
#   to_scale = 1} else {
#     to_scale = as.numeric(paste0(scaling_factor))
#   }
# scaling_factor=0.450328805
# control_file <- import.bed(WT_input)
# d <- list(treated_samp_file,control_file)
# control_label="WT_input"
# names(d) <- c(treated_samp_label,control_label)
# keep <- fread(chrom.sizes) %>%
#   setNames(c("chr","seqlength"))
# gn <- keep %>% {Seqinfo(.$chr, .$seqlength)}
# # loop through each file and keep only signal for each 10kb bin
# raw <- lapply(d,function(x){
#   inp <- as.data.frame(x) %>% dplyr::select("seqnames","start","end","name")
#   colnames(inp) <-c('chr', 'start', 'end', 'score')
#   out <- dplyr::semi_join(inp,keep,by="chr")
#   final <- out$score
# })
# t <- d[[treated_samp_label]] %>% as.data.frame()
# #   makeGRangesFromDataFrame()
# bs_final <- d[[treated_samp_label]] %>% as.data.frame() %>% dplyr::select(1:3) %>%
#   setNames(c('chr','start','end')) %>%
#   dplyr::semi_join(keep,by="chr") %>%
#   dplyr::filter(chr != "chrM") %>%
#   makeGRangesFromDataFrame()
# # get max values
# mxs_final <- bind_cols(raw) %>% apply(1, max)
# # kfinal
# k_final <- mxs_final > raw_count_cutoff & !overlapsAny(bs_final,bl)
# test <- d[[WT]]
# test <- t[k_final]
# raw <- lapply(d,function(x){
#   inp <- as.data.frame(x) %>% dplyr::select("seqnames","start","end","name")
#   colnames(inp) <-c('chr', 'start', 'end', 'score')
#   out <- dplyr::semi_join(inp,keep,by="chr")
#   out <- makeGRangesFromDataFrame(out,keep.extra.columns = TRUE)
#   out <- out[k_final]
# })
# pre_bl <- lapply(d,function(x){
#   inp <- as.data.frame(x) %>% dplyr::select("seqnames","start","end","name")
#   colnames(inp) <-c('chr', 'start', 'end', 'score')
#   imd <- dplyr::semi_join(inp,keep,by="chr")
#   out <- imd[,1:3] %>% mutate(start = start + 1)
#   out <- makeGRangesFromDataFrame(out,keep.extra.columns = TRUE,seqinfo = gn)
#   out <- out[k_final]
#   # out <- out %>%
#   #   makeGRangesFromDataFrame(seqinfo = gn)
# })
# bw <- lapply(d,function(x){
#   inp <- as.data.frame(x) %>% dplyr::select("seqnames","start","end","name")
#   colnames(inp) <-c('chr', 'start', 'end', 'score')
#   imd <- dplyr::semi_join(inp,keep,by="chr")
#   out <- imd[,1:3] %>% mutate(start = start + 1)
#   out <- out %>%
#     makeGRangesFromDataFrame(seqinfo = gn)
#   out <- out[k_final]
#   # out <- out %>%
#   #   makeGRangesFromDataFrame(seqinfo = gn)
#   k <- !overlapsAny(out, bl)
#   out <- out[!overlapsAny(out, bl)]
# })
# raw[[treated_samp_label]]$score <- as.numeric(raw[[treated_samp_label]]$score)
# raw[[control_label]]$score <- as.numeric(raw[[control_label]]$score)
# to_scale=as.numeric(0.450328805)
# treated_samp_library_size=64093770
# control_library_size=52047022
# bw[[treated_samp_label]]$score <- (log2(((raw[[treated_samp_label]]$score*to_scale)/treated_samp_library_size + addition)/(raw[[control_label]]$score/control_library_size + addition)))
# bw[[treated_samp_label]] <- bw[[treated_samp_label]][!is.na(bw[[treated_samp_label]]$score)]
# out_dir=paste0("~/Documents/ChIPbinner/tests/testthat")
# cell_line=paste0("Cal27")
# mark=paste0("H3K36me2")
# return(export.bw(bw[[treated_samp_label]],con=sprintf("%s/%s/%s.%s.%s%snorm.bw",out_dir,mark,cell_line,treated_samp_label,mark,window_size)))
#### running the function ####
# generate normalized bigwig for WT sample
norm_bw(out_dir = "~/Documents/ChIPbinner/data/HNSCC_ChIPseq_H3K36me2",
        chromSizes = hg38_chrom_sizes,
        blacklist = blacklisted_regions,
        cell_line = "Cal27",
        histone_mark = "H3K36me2",
        treated_samp_label = "WT",
        treated_samp_file= WT,
        treated_samp_library_size = 64093770,
        use_input = TRUE,
        control_label = "WT_input",
        control_file = WT_input,
        control_library_size = 52047022,
        raw_count_cutoff = 0,
        window_size = 10,
        addition = 1e-15,
        scaling_factor = 0.450328805)
# generate normalized bigwig for NSD1-KO sample
norm_bw(out_dir = "~/Documents/ChIPbinner/data/HNSCC_ChIPseq_H3K36me2",
        chromSizes = hg38_chrom_sizes,
        blacklist = blacklisted_regions,
        cell_line = "Cal27",
        histone_mark = "H3K36me2",
        treated_samp_label = "NSD1_KO",
        treated_samp_file= NSD1_KO,
        treated_samp_library_size = 18598272,
        use_input = TRUE,
        control_label = "NSD1_KO_input",
        control_file = NSD1_KO_input,
        control_library_size = 18674189,
        raw_count_cutoff = 0,
        window_size = 10,
        addition = 1e-15,
        scaling_factor = 0.192272095)
# generate genic/intergenic scatterplot
genic_intergenic_scatterplots(path_for_norm_bw = "~/Documents/ChIPbinner/data/HNSCC_ChIPseq_H3K36me2",
  out_dir = "~/Documents/ChIPbinner/data/HNSCC_ChIPseq_H3K36me2/results_genic_intergenic_scatterplots",
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
