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
#' @param wd
#' @param path
#' @param rep1
#' @param rep2
#' @param merged_sample_name
#' @param cell_line
#' @param histone_mark
#' @param window_size
#'
#' @return
#' @export
#'
#' @examples
merge_norm_bw <- function(wd,path,rep1,rep2,merged_sample_name,cell_line,histone_mark,window_size){
# parameters
#need to set wd
setwd(wd)
getwd()
# params
window_size=paste0(".",window_size,"kb.")
path <- paste0(path)
# sample info
cell_line=paste0(cell_line)
mark=paste0(histone_mark)
rep1=paste0(rep1)
rep2=paste0(rep2)
samp=paste0(merged_sample_name)
# list files
s <- list.files(path = path, pattern = window_size,full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c("line","samp","mark",NA,NA,NA),'\\.',F) %>%
  mutate(f = file.path(path, f))
d <- deframe(s[,c('samp', 'f')])
# loop through each file
bw <- lapply(d,function(x){
  inp <- import.bw(x)
  return(inp)
})
# take template
merged_bw <- bw[[rep1]]
merged_bw$score <- NULL
#list scores
l <- list(bw[[rep1]]$score,bw[[rep2]]$score)
# Calculate the mean for each row
row_means <- rowMeans(do.call(cbind, l))
# Assign the mean values to the 'score' column in 'merged_bw'
merged_bw$score <- row_means
# Create the 'norm.bw' directory if it doesn't exist
# if (!file.exists("data/norm.bw/merged_reps")) {
#   dir.create("data/norm.bw/merged_reps")
# }
out_dir = paste0(wd,"/data/norm.bw/merged")
if (!file.exists(out_dir)) {
  dir.create(out_dir)
}
final_out_dir = paste0(wd,"/data/norm.bw/merged/",mark)
if (!file.exists(final_out_dir)) {
  dir.create(final_out_dir)
}
# export bw
export.bw(merged_bw,con=sprintf("data/norm.bw/merged/%s/%s.%s.%s%snorm.bw",mark,cell_line,samp,mark,window_size))
}
