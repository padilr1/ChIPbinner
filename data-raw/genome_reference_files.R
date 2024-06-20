## code to prepare `genome_reference_files` dataset goes here

usethis::use_data(genome_reference_files, overwrite = TRUE)

library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(data.table)

hg38_chrom_sizes <- data.table::fread("~/Documents/ChIPbinner/inst/extdata/hg38_chrom.sizes.gz")
mm10_chrom_sizes <- data.table::fread("~/Documents/ChIPbinner/inst/extdata/mm10_chrom.sizes.gz")

hg38_bl <- rtracklayer::import.bed("~/Documents/ChIPbinner/inst/extdata/hg38_blacklist.bed.gz")
mm10_bl <- rtracklayer::import.bed("~/Documents/ChIPbinner/inst/extdata/mm10_blacklist.bed.gz")

hg38_gene <- rtracklayer::import.bed("~/Documents/ChIPbinner/inst/extdata/hg38_gene.bed.gz")
mm10_gene <- rtracklayer::import.bed("~/Documents/ChIPbinner/inst/extdata/mm10_gene.bed.gz")

hg38_igr <- rtracklayer::import.bed("~/Documents/ChIPbinner/inst/extdata/hg38_intergenic.bed.gz")
mm10_igr <- rtracklayer::import.bed("~/Documents/ChIPbinner/inst/extdata/mm10_intergenic.bed.gz")

# "gzip", "bzip2", or "xz"
# internal
usethis::use_data(hg38_bl,mm10_bl,mm10_chrom_sizes,hg38_chrom_sizes,hg38_gene,hg38_igr,mm10_gene,mm10_igr,internal=TRUE,compress="xz",overwrite = TRUE)

