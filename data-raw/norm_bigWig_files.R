## code to prepare `norm_bigWig_files` dataset goes here

usethis::use_data(norm_bigWig_files, overwrite = TRUE)

library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(data.table)

Cal27_WT_H3K36me2_10kb_norm_bw <- import.bw("~/Documents/ChIPbinner/inst/extdata/norm_bw/Cal27.WT.H3K36me2.10kb.norm.bw")
Cal27_NSD1_KO_H3K36me2_10kb_norm_bw <- import.bw("~/Documents/ChIPbinner/inst/extdata/norm_bw/Cal27.NSD1_KO.H3K36me2.10kb.norm.bw")

usethis::use_data(Cal27_WT_H3K36me2_10kb_norm_bw,internal=FALSE,compress="xz",overwrite = TRUE)
usethis::use_data(Cal27_NSD1_KO_H3K36me2_10kb_norm_bw,internal=FALSE,compress="xz",overwrite = TRUE)

# put downsampled norm bigwig files in /data
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
downsampled.Cal27.WT.H3K36me2.10kb.norm_bw <- loadRData("~/Documents/ChIPbinner/tests/testthat/testdata/downsampled.Cal27.WT.H3K36me2.10kb.bed.norm_bw.rda")
downsampled.Cal27.NSD1_KO.H3K36me2.10kb.norm_bw <- loadRData("~/Documents/ChIPbinner/tests/testthat/testdata/downsampled.Cal27.NSD1_KO.H3K36me2.10kb.bed.norm_bw.rda")
usethis::use_data(downsampled.Cal27.WT.H3K36me2.10kb.norm_bw,downsampled.Cal27.NSD1_KO.H3K36me2.10kb.norm_bw,internal=FALSE,compress="xz",overwrite = TRUE)

