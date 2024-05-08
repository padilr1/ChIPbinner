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
# testing
bw <- import.bed("~/Documents/HNSCC_K36me2/ref/blacklist.bed") %>% as.data.frame()
#
wd = "~/Documents/HNSCC_K36me2"
chromSizes = "~/Documents/HNSCC_K36me2/ref/hg38.chrom.sizes"
binned_files_path = "~/Documents/HNSCC_K36me2/data/10kb.bins/March_2024_batch/Cal27/H3K9me3"
treatment = "PA"
chipLibrarySize = 36922491
cell_line = "Cal27"
histone_mark = "H3K9me3"
blacklist = "~/Documents/HNSCC_K36me2/ref/blacklist.bed"
raw_count_cutoff = 0
use_input = FALSE
window_size = 10
addition = 1e-15
multiplier = 1
scaling_factor = 0.42
control=NULL
inputLibrarySize=NULL
# function
generate_norm_bw <- function(wd,chromSizes,binned_files_path,treatment,control=NULL,chipLibrarySize,inputLibrarySize=NULL,cell_line,histone_mark,blacklist,addition,raw_count_cutoff,use_input,window_size,scaling_factor=NULL,multiplier){
#need to set wd
setwd(wd)
getwd()
# params
chrom.sizes=paste0(chromSizes)
binned_files_path <- paste0(binned_files_path)
addition <- as.numeric(paste0(addition))
raw_count_cutoff <- as.numeric(paste0(raw_count_cutoff))
# sample info
cell_line=paste0(cell_line)
mark=paste0(histone_mark)
samp=paste0(treatment) #watch spaces in names
#### reference files ####
# blacklist
blacklist=paste0(blacklist)
bl <- import.bed(blacklist)
# window size
window_size=paste0(".",window_size,"kb.")
# scaling factor
if(is.null(scaling_factor)){
  to_scale = 1} else {
    to_scale = as.numeric(paste0(scaling_factor))
  }
# if using input
if(is.null(control)){
  print("Not using input")
  chipLibrary <- as.numeric(paste0(chipLibrarySize))
  print(paste0("Sample sequencing depth=",chipLibrary))
  print(paste0("Treatment sample=",samp))
  s <- list.files(path = binned_files_path, pattern = window_size,full.names = FALSE,recursive = FALSE) %>%
    tibble(f = .) %>%
    separate(f,c("line","samp","mark",NA,NA),'\\.',F) %>%
    mutate(f = file.path(binned_files_path, f)) %>%
    dplyr::filter(line==cell_line) %>%
    dplyr::filter(samp==samp)
  d <- deframe(s[,c('samp', 'f')])
} else {
  input=paste0(control)
  chipLibrary <- as.numeric(paste0(chipLibrarySize))
  inputLibrary <- as.numeric(paste0(inputLibrarySize))
  print(paste0("Control sample=",control))
  print(paste0("Treatment sample=",samp))
  print(paste0("Sample sequencing depth=",chipLibrary))
  print(paste0("Input/IgG sequencing depth=",inputLibrary))
  s <- list.files(path = binned_files_path, pattern = window_size,full.names = FALSE,recursive = FALSE) %>%
    tibble(f = .) %>%
    separate(f,c("line","samp","mark",NA,NA),'\\.',F) %>%
    mutate(f = file.path(binned_files_path, f)) %>%
    dplyr::filter(line==cell_line) %>%
    dplyr::filter(samp==samp | samp==input)
  d <- deframe(s[,c('samp', 'f')])
}
###################################
# generate filter for k_final
s <- list.files(path = binned_files_path, pattern = window_size,full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c("line","samp","mark",NA,NA),'\\.',F) %>%
  mutate(f = file.path(binned_files_path, f))
# deframe
d <- deframe(s[,c('samp', 'f')])
# read in chrom sizes and chr list
keep <- fread(chrom.sizes) %>%
  setNames(c("chr","seqlength"))
gn <- keep %>% {Seqinfo(.$chr, .$seqlength)}
# loop through each file and keep only signal for each 10kb bin
raw <- lapply(d,function(x){
  inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
  out <- dplyr::semi_join(inp,keep,by="chr")
  final <- out$score
})
# read in bin sizes
bs_final <- import.bed(s$f[s$samp == samp]) %>% as.data.frame() %>% dplyr::select(1:3) %>%
  setNames(c('chr','start','end')) %>%
  dplyr::semi_join(keep,by="chr") %>%
  dplyr::filter(chr != "chrM") %>%
  makeGRangesFromDataFrame()
# get max values
mxs_final <- bind_cols(raw) %>% apply(1, max)
# kfinal
k_final <- mxs_final > raw_count_cutoff & !overlapsAny(bs_final,bl)
# k_final <- !overlapsAny(bs_final,bl)

####################################

# load exclusion factor
# load("")
# list files
# s <- list.files(path = binned_files_path, pattern = window_size,full.names = FALSE,recursive = FALSE) %>%
#   tibble(f = .) %>%
#   separate(f,c("line","samp","mark",NA,NA),'\\.',F) %>%
#   mutate(f = file.path(binned_files_path, f)) %>%
#   dplyr::filter(line==cell_line) %>%
#   dplyr::filter(samp==samp | samp==input)
# d <- deframe(s[,c('samp', 'f')])
# read in chrom sizes and chr list
keep <- fread(chrom.sizes) %>%
  setNames(c("chr","seqlength"))
gn <- keep %>% {Seqinfo(.$chr, .$seqlength)}
# loop through each file
raw <- lapply(d,function(x){
  inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
  out <- dplyr::semi_join(inp,keep,by="chr")
  out <- out[k_final]
})
pre_bl <- lapply(d,function(x){
  inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
  imd <- dplyr::semi_join(inp,keep,by="chr")
  out <- imd[,1:3] %>% mutate(start = start + 1)
  out <- out[k_final]
  out <- out %>%
    makeGRangesFromDataFrame(seqinfo = gn)
})
bw <- lapply(d,function(x){
  inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
  imd <- dplyr::semi_join(inp,keep,by="chr")
  out <- imd[,1:3] %>% mutate(start = start + 1)
  out <- out[k_final]
  out <- out %>%
    makeGRangesFromDataFrame(seqinfo = gn)
  k <- !overlapsAny(out, bl)
  out <- out[!overlapsAny(out, bl)]
})
# normalize by input and scale by a factor
# bw[[samp]]$score <- (log2(((raw[[samp]]$score/chipLibrary) + addition))*1e6)[!overlapsAny(pre_bl[[samp]], bl)]
# bw[[samp]]$score <- log2(((raw[[samp]]$score/chipLibrary) + 1e-15) / (raw[[input]]$score/inputLibrary + 1e-15))[!overlapsAny(pre_bl[[samp]], bl)]
multiplier=as.numeric(paste0(multiplier))
if (use_input==TRUE){
  bw[[samp]]$score <- ((log2((raw[[samp]]$score/chipLibrary + addition)/(raw[[input]]$score/inputLibrary + addition) ))+5)*multiplier*to_scale
} else{
  bw[[samp]]$score <- raw[[samp]]$score * to_scale * 100 (chipLibrary / length(bw))
  # bw[[samp]]$score <- (log2(((raw[[samp]]$score/chipLibrary) + addition))*multiplier)*to_scale
}
# bw[[samp]]$score <- (log2(((raw[[samp]]$score/chipLibrary) + addition))*1e6)
# drop NA values
bw[[samp]] <- bw[[samp]][!is.na(bw[[samp]]$score)]
# Create the 'norm.bw' directory if it doesn't exist
out_dir = paste0("data/norm.bw/",mark)
if (!file.exists(out_dir)) {
  dir.create(out_dir)
}
# export bw
return(export.bw(bw[[samp]],con=sprintf("data/norm.bw/%s/%s.%s.%s%snorm.bw",mark,cell_line,samp,mark,window_size)))
}
