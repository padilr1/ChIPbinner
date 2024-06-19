test_that("running enrichment and depletion analysis works for specified clusters", {
  library(tidyverse)
  library(readr)
  library(rtracklayer)
  library(GenomicRanges)
  library(LOLA)
  library(highcharter)
  enrich_clust(gene=system.file("extdata", "hg38_gene.bed.gz", package = "ChIPbinner"),
               intergenic=system.file("extdata", "hg38_intergenic.bed.gz", package = "ChIPbinner"),
               assembly="hg38",
               annotated_clusters=testthat::test_path("testdata","cons.Cal27.WT.NSD1_KO.H3K36me2.10kb.rda"),
               query_cluster="B",
               pooled_bed_file=testthat::test_path("testdata","Cal27.WT.NSD1_KO.H3K36me2.10kb.pooled.bed"),
               functional_db="ccre",
               region="genome_wide",
               cores=1,
               n_elements=7,
               cutoff_for_overlap = 100,
               file_plot_name="ccre.genome_wide.enrichment_depletion",
               output_table_name="ccre.genome_wide.enrichment_depletion",
               width_of_plot=7,
               height_of_plot=3,
               out_dir=testthat::test_path("testdata"))
  # enrich_clust(gene=system.file("extdata", "hg38_gene.bed", package = "ChIPbinner"),
  #              intergenic=system.file("extdata", "hg38_intergenic.bed", package = "ChIPbinner"),
  #              assembly="hg38",
  #              annotated_clusters=testthat::test_path("testdata","cons.Cal27.WT.NSD1_KO.H3K36me2.10kb.rda"),
  #              query_cluster="B",
  #              pooled_bed_file=testthat::test_path("testdata","Cal27.WT.NSD1_KO.H3K36me2.10kb.pooled.bed"),
  #              functional_db="ensembl",
  #              region="intergenic",
  #              cores=6,
  #              n_elements=7,
  #              cutoff_for_overlap = 100,
  #              file_plot_name="igr.enrichment_depletion",
  #              output_table_name="igr.enrichment_depletion",
  #              width_of_plot=8,
  #              height_of_plot=4,
  #              out_dir=testthat::test_path("testdata"))
  # enrich_clust(gene=system.file("extdata", "hg38_gene.bed", package = "ChIPbinner"),
  #              intergenic=system.file("extdata", "hg38_intergenic.bed", package = "ChIPbinner"),
  #              assembly="hg38",
  #              annotated_clusters=testthat::test_path("testdata","cons.Cal27.WT.NSD1_KO.H3K36me2.10kb.rda"),
  #              query_cluster="B",
  #              pooled_bed_file=testthat::test_path("testdata","Cal27.WT.NSD1_KO.H3K36me2.10kb.pooled.bed"),
  #              functional_db="ensembl",
  #              region="genic",
  #              cores=6,
  #              n_elements=7,
  #              cutoff_for_overlap = 100,
  #              file_plot_name="genic.enrichment_depletion",
  #              output_table_name="genic.enrichment_depletion",
  #              width_of_plot=8,
  #              height_of_plot=4,
  #              out_dir=testthat::test_path("testdata"))
  expect_snapshot_output(x = "Enrichment/depletion output successfully generated!")
})
