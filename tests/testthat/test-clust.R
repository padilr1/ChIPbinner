test_that("Density-based clustering works using HDBSCAN", {
  library(reticulate)
  # use_python("/Users/padilr1/opt/anaconda3/envs/r_env_V2/bin/python")
  # reticulate::source_python("~/Documents/ChIPbinner/inst/python/clus.py")
  # run HDBSCAN
  clust(
    output_file_name = "sampled.Cal27.WT.NSD1_KO.H3K36me2.10kb",
    out_dir = testthat::test_path("testdata"),
    matrix_file = "sampled.Cal27.WT.NSD1_KO.H3K36me2.10kb.mat.csv",
    minpts = 100,
    minsamps = 100,
    cores = 6
  )
  expect_snapshot_output(x = "Density-based clusters generated using HDBSCAN!")
})
