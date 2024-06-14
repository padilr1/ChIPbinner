library(reticulate)
clust(output_file_name = "sampled.Cal27.WT.NSD1_KO.H3K36me2.10kb",
  out_dir = system.file("extdata", package = "ChIPbinner"),
  matrix_file = "sampled.Cal27.WT.NSD1_KO.H3K36me2.10kb.mat.csv",
  minpts = 100,
  minsamps = 100,
  cores = 6)
