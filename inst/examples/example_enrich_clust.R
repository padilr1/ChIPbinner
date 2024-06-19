#!/usr/bin/env Rscript
library(tidyverse)
library(readr)
library(rtracklayer)
library(GenomicRanges)
library(LOLA)
library(highcharter)
enrich_clust(gene=system.file("extdata", "hg38_gene.bed.gz", package = "ChIPbinner"),
             intergenic=system.file("extdata", "hg38_intergenic.bed.gz", package = "ChIPbinner"),
             assembly="hg38",
             annotated_clusters=system.file("extdata","cons.Cal27.WT.NSD1_KO.H3K36me2.10kb.rda", package = "ChIPbinner"),
             query_cluster="B",
             pooled_bed_file=system.file("extdata","Cal27.WT.NSD1_KO.H3K36me2.10kb.pooled.bed", package = "ChIPbinner"),
             functional_db="ccre",
             region="genic",
             cores=1,
             n_elements=7,
             cutoff_for_overlap = 100,
             file_plot_name="ccre.genic.enrichment_depletion",
             output_table_name="ccre.genic.enrichment_depletion",
             width_of_plot=7,
             height_of_plot=3,
             out_dir=system.file("extdata", package = "ChIPbinner"))
