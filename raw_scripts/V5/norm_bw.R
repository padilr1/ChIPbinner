#!/usr/bin/env Rscript
#+ message = FALSE, warning = FALSE
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(lattice)
library(gridExtra)
library(DiffBind)
library(patchwork)
#' Title
#'
#' @param out_dir
#' @param chromSizes
#' @param blacklist
#' @param cell_line
#' @param histone_mark
#' @param treated_samp_label
#' @param treated_samp_file
#' @param treated_samp_library_size
#' @param control_label
#' @param control_file
#' @param control_library_size
#' @param use_input
#' @param window_size
#' @param addition
#' @param raw_count_cutoff
#' @param scaling_factor
#'
#' @return
#' @export
#'
#' @examples
norm_bw <- function(out_dir,
                    chromSizes,
                    blacklist,
                    cell_line,
                    histone_mark,
                    treated_samp_label,
                    treated_samp_file,
                    treated_samp_library_size,
                    control_label=NULL,
                    control_file=NULL,
                    control_library_size=NULL,
                    use_input,
                    window_size,
                    addition,
                    raw_count_cutoff,
                    scaling_factor=NULL){
# reference files
## blacklist
blacklist=paste0(blacklist)
bl <- import.bed(blacklist)
## chrom sizes
chrom.sizes=paste0(chromSizes)
# cell line info
cell_line=paste0(cell_line)
mark=paste0(histone_mark)
treated_samp_label=paste0(treated_samp_label) #watch spaces in names
treated_samp_file=import.bed(paste0(treated_samp_file))
# other parameters
## addition to avoid division by 0
addition <- as.numeric(paste0(addition))
## raw count cut-off to remove bins with low counts across samples
raw_count_cutoff <- as.numeric(paste0(raw_count_cutoff))
## window size of the bins
window_size=paste0(".",window_size,"kb.")
## scaling factor
if(is.null(scaling_factor)){
  to_scale = 1} else {
    to_scale = as.numeric(paste0(scaling_factor))
  }
# if using input
if(is.null(control_file)){
  print("Not using input")
  treated_samp_library_size <- as.numeric(paste0(treated_samp_library_size))
  print(paste0("Treated sample sequencing depth=",treated_samp_library_size))
  print(paste0("Treated sample=",treated_samp_label))
} else {
  treated_samp_library_size <- as.numeric(paste0(treated_samp_library_size))
  control_library_size <- as.numeric(paste0(control_library_size))
  print(paste0("Treated sample=",treated_samp_label))
  print(paste0("Treated sample sequencing depth=",treated_samp_library_size))
  print(paste0("Control sample=",control_label))
  print(paste0("Control sample sequencing depth=",control_library_size))
  control_file <- import.bed(control_file)
  d <- list(treated_samp_file,control_file)
  names(d) <- c(treated_samp_label,control_label)
}
### deriving cut-off for raw counts ###
# read in chrom sizes
keep <- fread(chrom.sizes) %>%
  setNames(c("chr","seqlength"))
gn <- keep %>% {Seqinfo(.$chr, .$seqlength)}
# loop through each file and keep only signal for each 10kb bin
raw <- lapply(d,function(x){
  inp <- as.data.frame(x) %>% dplyr::select("seqnames","start","end","name")
  colnames(inp) <-c('chr', 'start', 'end', 'score')
  out <- dplyr::semi_join(inp,keep,by="chr")
  final <- out$score
})
# read in bin sizes
bs_final <- d[[treated_samp_label]] %>% as.data.frame() %>% dplyr::select(1:3) %>%
  setNames(c('chr','start','end')) %>%
  dplyr::semi_join(keep,by="chr") %>%
  dplyr::filter(chr != "chrM") %>%
  makeGRangesFromDataFrame()
# get max values
mxs_final <- bind_cols(raw) %>% apply(1, max)
# load exclusion factor based on raw count cut-off and overlap with blacklisted region
k_final <- mxs_final > raw_count_cutoff & !overlapsAny(bs_final,bl)
### generate bigwig template ###
# loop through each file and remove bins based on exclusion factor (k_final)
# for the raw bins, we're just keeping the scores
raw <- lapply(d,function(x){
  inp <- as.data.frame(x) %>% dplyr::select("seqnames","start","end","name")
  colnames(inp) <-c('chr', 'start', 'end', 'score')
  out <- dplyr::semi_join(inp,keep,by="chr")
  out <- makeGRangesFromDataFrame(out,keep.extra.columns = TRUE)
  out <- out[k_final]
})
# pre_bl <- lapply(d,function(x){
#   inp <- as.data.frame(x) %>% dplyr::select("seqnames","start","end","name")
#   colnames(inp) <-c('chr', 'start', 'end', 'score')
#   imd <- dplyr::semi_join(inp,keep,by="chr")
#   out <- imd[,1:3] %>% mutate(start = start + 1)
#   out <- makeGRangesFromDataFrame(out,keep.extra.columns = TRUE,seqinfo = gn)
#   out <- out[k_final]
# })
# final bins without the score. we will merge the scores after they are normalized and scaled
bw <- lapply(d,function(x){
  inp <- as.data.frame(x) %>% dplyr::select("seqnames","start","end","name")
  colnames(inp) <-c('chr', 'start', 'end', 'score')
  imd <- dplyr::semi_join(inp,keep,by="chr")
  out <- imd[,1:3] %>% mutate(start = start + 1)
  out <- out %>%
    makeGRangesFromDataFrame(seqinfo = gn)
  out <- out[k_final]
  # out <- out %>%
  #   makeGRangesFromDataFrame(seqinfo = gn)
  k <- !overlapsAny(out, bl)
  out <- out[!overlapsAny(out, bl)]
})
# ensure the scores per bin are numeric
raw[[treated_samp_label]]$score <- as.numeric(raw[[treated_samp_label]]$score)
raw[[control_label]]$score <- as.numeric(raw[[control_label]]$score)
# normalize scores per bin by the library depth and scale if necessary
if (use_input==TRUE){
  bw[[treated_samp_label]]$score <- (log2(((raw[[treated_samp_label]]$score*to_scale)/treated_samp_library_size + addition)/(raw[[control_label]]$score/control_library_size + addition)))
} else{
  bw[[treated_samp_label]]$score <- (log2(((raw[[treated_samp_label]]$score*to_scale)/treated_samp_library_size) + addition))
}
bw[[treated_samp_label]] <- bw[[treated_samp_label]][!is.na(bw[[treated_samp_label]]$score)]
# output directory
out_dir = paste0(out_dir)
# output bigwig file
return(export.bw(bw[[treated_samp_label]],con=sprintf("%s/%s.%s.%s%snorm.bw",out_dir,cell_line,treated_samp_label,mark,window_size)))
print("Normalized bigWig file created!")
}
