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
# function
generate_norm_bw <- function(wd,chromSizes,binned_files_path,treatment,control,,cell_line,histone_mark,blacklist,addition,raw_count_cutoff,use_input,windowed_bed,window_size){
#need to set wd
setwd(wd)
getwd()
# params
chrom.sizes=paste0(chromSizes)
path <- paste0(binned_files_path)
addition <- as.numeric(paste0(addition))
raw_count_cutoff <- as.numeric(paste0(raw_count_cutoff))
# sample info
cell_line=paste0(cell_line)
mark=paste0(histone_mark)
samp=paste0(treatment)
input=paste0(control) #watch spaces in names
windowed_bed=paste0(windowed_bed)
# RX #
########## make sure column name for chip library is "chip" and column name for input library is "inp", as well as "rx_ratio" for rx_ratio #########
rx=paste0(table)
rx_samp=paste0(table_samp_name)
rx <- fread(rx)
chipLibrary <- rx$CR_MappedReads[rx$samp == rx_samp]
inputLibrary <- rx$IgG_MappedReads[rx$samp == rx_samp]
# rxRatio <- rx$rx_ratio[rx$samp == rx_samp]
# scalingFactor <- rxRatio
# blacklist
blacklist=paste0(blacklist)
bl <- import.bed(blacklist)
print(paste0("Sample sequencing depth=",chipLibrary))
print(paste0("Input/IgG sequencing depth=",inputLibrary))
# window size
window_size=paste0(".",window_size,"kb.")
###################################

# generate filter for k_final
s <- list.files(path = path, pattern = window_size,full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c("line","samp","mark",NA,NA),'\\.',F) %>%
  mutate(f = file.path(path, f))
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
# read in bin
bs_final <- import.bed(windowed_bed) %>% as.data.frame() %>% dplyr::select(1:3) %>%
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
s <- list.files(path = path, pattern = window_size,full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c("line","samp","mark",NA,NA),'\\.',F) %>%
  mutate(f = file.path(path, f)) %>%
  dplyr::filter(line==cell_line) %>%
  dplyr::filter(samp==samp | samp==input)
d <- deframe(s[,c('samp', 'f')])
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
if (use_input==TRUE){
  bw[[samp]]$score <- (log2((raw[[samp]]$score/chipLibrary + addition)/(raw[[input]]$score/inputLibrary + addition) ))*1e6
} else{
  bw[[samp]]$score <- (log2(((raw[[samp]]$score/chipLibrary) + addition))*1e6)
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
return(export.bw(bw[[samp]],con=sprintf("data/norm.bw/%s/%s.%s.%s.%s.norm.bw",mark,cell_line,samp,mark,window_size)))
}
