test_that("generating normalized bigwig files from binned BED files works", {
  # generate normalized bigwig for WT sample
  norm_bw(out_dir = testthat::test_path("testdata"),
          genome_assembly = "hg38",
          immunoprecipitated_binned_file = system.file("extdata", "downsampled.Cal27.WT.H3K36me2.10kb.bed.gz", package = "ChIPbinner"),
          use_input = TRUE,
          input_binned_file = system.file("extdata", "downsampled.Cal27.WT_input.H3K36me2.10kb.bed.gz", package = "ChIPbinner"),
          raw_count_cutoff = 0,
          pseudocount = 1e-3,
          scaling_factor = 0.450328805)
  # generate normalized bigwig for NSD1-KO sample
  # norm_bw(out_dir = testthat::test_path("testdata"),
  #         genome_assembly = "hg38",
  #         chip_samp_binned_file = system.file("extdata", "downsampled.Cal27.NSD1_KO.H3K36me2.10kb.bed.gz", package = "ChIPbinner"),
  #         chip_samp_library_size = 18598272,
  #         use_control = TRUE,
  #         control_binned_file = system.file("extdata", "downsampled.Cal27.NSD1_KO_input.H3K36me2.10kb.bed.gz", package = "ChIPbinner"),
  #         control_library_size = 18674189,
  #         raw_count_cutoff = 0,
  #         pseudocount = 1e-15,
  #         scaling_factor = 0.192272095)
  expect_snapshot_output(x = "Normalized bigWig file created!")
})
