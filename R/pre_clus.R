#!/usr/bin/env Rscript
#+ message = FALSE, warning = FALSE
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(lattice)
library(gridExtra)
library(DiffBind)
library(idr2d)
library(patchwork)
# start of function
#' Title
#'
#' @param wd
#' @param path_for_norm_bw
#' @param window_size
#' @param treatment
#' @param control
#' @param cell_line
#' @param histone_mark
#'
#' @return
#' @export
#'
#' @examples
pre_clus <- function(wd,path_for_norm_bw,window_size,treatment,control,cell_line,histone_mark){
# parameters
setwd(wd)
getwd()
path <- paste0(path_for_norm_bw)
# pattern = paste0(typeof_file)
window_size=paste0(".",window_size,"kb.")
# gene <- import.bed('ensembl/gene.bed')
# igr <- import.bed('ensembl/intergenic.bed')
control=paste0(control)
test=paste0(treatment)
cell_line=paste0(cell_line)
mark=paste0(histone_mark)
#
s <- list.files(path = path, pattern = window_size,full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c("line","samp","mark",NA,NA,NA),'\\.',F) %>%
  mutate(f = file.path(path, f)) %>%
  dplyr::filter(line == cell_line) %>%
  dplyr::filter(mark == mark) %>%
  dplyr::filter(samp == control | samp == test)
# filter(cond != 'MT') %>%
# filter(samp != "Cal27_KO17") %>%
# filter(samp != "Cal27_NSD2OE")
#
#
# NEED TO ORDER SAMPLE FIRST THEN PARENTAL #
odr <- c(test,control)
s <- s %>%
  dplyr::slice(match(odr,samp))
#
d <- deframe(s[,c('samp', 'f')]) %>%
  lapply(import.bw) #original ranges for KO = 283852
#
lapply(d,length)
#
#
d[[control]] <- subsetByOverlaps(d[[control]],d[[test]])
d[[test]]<- subsetByOverlaps(d[[test]],d[[control]])
#
lapply(d,length)
#
#
#
r <- lapply(d, function(y) y[y$score != 0]) %>%
  Reduce(function(a, b) a[overlapsAny(a, b)], .) %>%
  granges()
#
#
cell_line <- as.character(s$line[[1]])
mark <- as.character(s$mark[[1]])
#
# Create the 'mat.csv' directory if it doesn't exist
if (!file.exists("data/mat.csv")) {
  dir.create("data/mat.csv")
}
# Create the 'pooled.bed' directory if it doesn't exist
if (!file.exists("data/pooled.bed")) {
  dir.create("data/pooled.bed")
}
# Create the figs directory if it doesn't exist
if (!file.exists("figs/")) {
  dir.create("figs/")
}
# split
# ok <- o$x > quantile(o$x, .01) &
# o$x < quantile(o$x, .99) &
#   o$y > quantile(o$y, .01) &
#   o$y < quantile(o$y, .99)

split(d, s$line) %>%
  lapply(function(x) {
    o <- lapply(x, function(y) {
      findOverlaps(r, y) %>%
        to() %>%
        {y[.]} %>%
        score()
    }) %>%
      bind_cols() %>%
      `names<-`(c('x', 'y'))
    ok <- o$x > quantile(o$x, .01) &
      o$x < quantile(o$x, .99) &
      o$y > quantile(o$y, .01) &
      o$y < quantile(o$y, .99)
    write_csv(o[ok,], sprintf('data/mat.csv/%s.%s.%s.%s%smat.csv',cell_line,control,test,mark,window_size), col_names = F)
    export.bed(r[ok], sprintf('data/pooled.bed/%s.%s.%s.%s%spooled.bed',cell_line,control,test,mark,window_size))
  })

lfc <- fread(sprintf('data/mat.csv/%s.%s.%s.%s%smat.csv',cell_line,control,test,mark,window_size))
# 1. Open jpeg file
jpeg(filename = sprintf('figs/%s.%s.%s.%s%ssmoothScatter.jpeg',cell_line,control,test,mark,window_size))
# 2. Create the plot
smoothScatter(y = lfc$V1,x=lfc$V2,xlab = control,ylab=test)
# 3. Close the file
dev.off()
return(list(s=s,mark=mark))
}
