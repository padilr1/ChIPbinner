## code to prepare `binned_BED_files` dataset goes here

usethis::use_data(binned_BED_files, overwrite = TRUE)

library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
# downsample bed files
set.seed(123)
WT <- fread("~/Documents/ChIPbinner/inst/extdata/Cal27.WT.H3K36me2.10kb.bed.gz") %>% slice_sample(n=1000,replace=FALSE) %>% `names<-`(c('chr', 'start','end','score')) %>% mutate(strand="*") %>% mutate(bin_id = paste0("bin_",1:nrow(.))) %>% dplyr::select(c("chr","start","end","bin_id","score","strand"))
write_tsv(WT,file="~/Documents/ChIPbinner/inst/extdata/downsampled.Cal27.WT.H3K36me2.10kb.bed",col_names = FALSE,quote = "none")
WT <- import.bed("downsampled.Cal27.WT.H3K36me2.10kb.bed.gz")
# subset by overlaps
## input
WT_input <-rtracklayer::import.bed("~/Documents/ChIPbinner/inst/extdata/Cal27.WT_input.H3K36me2.10kb.bed.gz") %>% subsetByOverlaps(.,WT)
export.bed(WT_input,con="~/Documents/ChIPbinner/inst/extdata/downsampled.Cal27.WT_input.H3K36me2.10kb.bed.gz")
## NSD1_KO
NSD1_KO <- rtracklayer::import.bed("~/Documents/ChIPbinner/inst/extdata/Cal27.NSD1_KO.H3K36me2.10kb.bed.gz") %>% subsetByOverlaps(.,WT)
export.bed(NSD1_KO,con="downsampled.Cal27.NSD1_KO.H3K36me2.10kb.bed.gz")
## NSD1_KO input
NSD1_KO_input <- rtracklayer::import.bed("~/Documents/ChIPbinner/inst/extdata/Cal27.NSD1_KO_input.H3K36me2.10kb.bed.gz") %>% subsetByOverlaps(.,WT)
export.bed(NSD1_KO_input,con="downsampled.Cal27.NSD1_KO_input.H3K36me2.10kb.bed.gz")
# re-read WT
WT_input <- rtracklayer::import.bed("~/Documents/ChIPbinner/inst/extdata/downsampled.Cal27.WT_input.H3K36me2.10kb.bed.gz")
WT <- rtracklayer::import.bed("~/Documents/ChIPbinner/inst/extdata/Cal27.WT.H3K36me2.10kb.bed.gz") %>% subsetByOverlaps(.,WT_input)
export.bed(WT,con="downsampled.Cal27.WT.H3K36me2.10kb.bed.gz")
