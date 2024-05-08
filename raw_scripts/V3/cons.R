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
# start of code
cons <- function(wd,treatment,control,cell_line,histone_mark,gene,intergenic,cluster_size,number_of_clusters,window_size){
# parameters
setwd(wd)
getwd()
control=paste0(control)
test=paste0(treatment)
cell_line=paste0(cell_line)
mark=paste0(histone_mark)
gene <- import.bed(gene)
igr <- import.bed(intergenic)
cluster_size <- paste0(cluster_size)
number_of_clusters <- as.integer(paste0(number_of_clusters))
window_size=paste0(".",window_size,"kb.")
# read in files
d <- read_csv(sprintf('data/mat.csv/%s.%s.%s.%s%smat.csv',cell_line,control,test,mark,window_size), col_names = F)
r <- import.bed(sprintf('data/pooled.bed/%s.%s.%s.%s%spooled.bed',cell_line,control,test,mark,window_size))
d$clu <- read_csv(file = sprintf('data/hdbscan.out/clus.%s.%s.%s.%s%s%s.%s.txt',cell_line,control,test,mark,window_size,cluster_size,cluster_size), col_names = F)$X1

# Create the directory if it doesn't exist
if (!file.exists("data/cons")) {
  dir.create("data/cons")
}

if ( number_of_clusters == 2 ) {
  ctr <- d %>%
    dplyr::filter(clu != -1) %>%
    group_by(clu) %>%
    summarise_all(mean) %>%
    mutate(mu = 0.5 * (X1 + X2)) %>%
    arrange(mu)

  a1 <- r[d$clu == ctr$clu[1]]
  b1 <- r[d$clu == ctr$clu[2]]

  cons <- list(A = list(a1),
               B = list(b1)) %>%
    lapply(function(x) {
      Reduce(function(a,b){a[overlapsAny(a,b)]}, x) %>%
        granges()
    })
  save(cons, file = (sprintf('data/cons/cons.%s.%s.%s.%s%srda',cell_line,control,test,mark,window_size)))
} else{
  ctr <- d %>%
    dplyr::filter(clu != -1) %>%
    group_by(clu) %>%
    summarise_all(mean) %>%
    mutate(mu = 0.5 * (X1 + X2)) %>%
    arrange(mu)

  a1 <- r[d$clu == ctr$clu[1]]
  b1 <- r[d$clu == ctr$clu[2]]
  c1 <- r[d$clu == ctr$clu[3]]

  cons <- list(A = list(a1),
               B = list(b1),
               C = list(c1)) %>%
    lapply(function(x) {
      Reduce(function(a,b){a[overlapsAny(a,b)]}, x) %>%
        granges()
    })
  save(cons, file = (sprintf('data/cons/cons.%s.%s.%s.%s%srda',cell_line,control,test,mark,window_size)))
}
}
