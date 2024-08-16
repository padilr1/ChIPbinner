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
  out_dir = system.file("extdata", package = "ChIPbinner")
)
