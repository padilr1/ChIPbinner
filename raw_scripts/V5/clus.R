#!/usr/bin/env Rscript
clus <- function(output_file_name,
                 out_dir,
                 matrix_file,
                 baseline_samp_label,
                 treated_samp_label,
                 minpts,
                 minsamps,
                 cores) {
  install_HDBSCAN <-
    function(...,
             envname = "r-HBDSCAN",
             new_env = identical(envname, "r-HBDSCAN")) {

      if(new_env && virtualenv_exists(envname))
        virtualenv_remove(envname)

      reticulate::py_install(packages = c("hdbscan","numpy","matplotlib","seaborn"), envname = envname, ...)

      .onLoad <- function(...) {
        use_virtualenv("r-HDBSCAN", required = FALSE)
      }
    }
  reticulate::source_python(system.file("python", "clus.py", package = "ChIPbinner"))
  # file parameters
  output_file_name <- paste0(output_file_name)
  # output directory
  setwd(paste0(out_dir))
  # samples info
  baseline_samp_label <- paste0(baseline_samp_label)
  treated_samp_label <- paste0(treated_samp_label)
  # input file
  mat_file <- paste0(matrix_file)
  # algorithm parameters
  minpts <- as.integer(paste0(minpts))
  minsamps <- as.integer(paste0(minsamps))
  # number of parallel cores to use
  cores <- as.integer(paste0(cores))
  # out_dir
  hdbscan_clustering(
    file_name = output_file_name,
    mat_file = mat_file,
    minpts = minpts,
    minsamps = minsamps,
    baseline_samp_label = baseline_samp_label,
    treated_samp_label = treated_samp_label,
    cores = cores
  )
  print("Density-based clusters generated using HDBSCAN!")
}
