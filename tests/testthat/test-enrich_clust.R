test_that("running enrichment and depletion analysis works for specified clusters", {
  library(tidyverse)
  library(readr)
  library(rtracklayer)
  library(GenomicRanges)
  library(LOLA)
  library(highcharter)
  # downsampled
  enrich_clust(
    genome_assembly = "hg38",
    annotated_clusters = system.file("extdata/data","downsampled_annotated_clusters.rda", package = "ChIPbinner"),
    query_cluster = "B",
    pooled_bed_file = system.file("extdata", "downsampled.Cal27.WT.NSD1_KO.H3K36me2.10kb_pooled.bed.gz", package = "ChIPbinner"),
    functional_db = system.file("extdata/data","downsampled_hg38_ccreDB.rda", package = "ChIPbinner"),
    region = "genome_wide",
    cores = 1,
    n_elements = 7,
    cutoff_for_overlap = 10,
    file_plot_name = "downsampled.ccre.genome_wide.enrichment_depletion",
    output_table_name = "downsampled.ccre.genome_wide.enrichment_depletion",
    width_of_plot = 7,
    height_of_plot = 3,
    out_dir = testthat::test_path("testdata")
  )
  # complete
  # enrich_clust(genome_assembly="hg38",
  #              annotated_clusters=testthat::test_path("testdata","cons.Cal27.WT.NSD1_KO.H3K36me2.10kb.rda"),
  #              query_cluster="B",
  #              pooled_bed_file=testthat::test_path("testdata","Cal27.WT.NSD1_KO.H3K36me2.10kb.pooled.bed"),
  #              functional_db=testthat::test_path("testdata","downsampled_hg38_ccreDB.rda"),
  #              region="genome_wide",
  #              cores=1,
  #              n_elements=7,
  #              cutoff_for_overlap = 100,
  #              file_plot_name="downsampled.ccre.genome_wide.enrichment_depletion",
  #              output_table_name="downsampled.ccre.genome_wide.enrichment_depletion",
  #              width_of_plot=7,
  #              height_of_plot=3,
  #              out_dir=testthat::test_path("testdata"))
  expect_snapshot_output(x = "Enrichment/depletion output successfully generated!")
})
