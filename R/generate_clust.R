#!/usr/bin/env Rscript
#' Title
#' @title Generate density-based clusters using the HDBSCAN algorithm.
#'
#' @description This function generates density-based clusters needed to generate the density-based scatterplots and furhter downstream analysis.
#' @param output_file_name The name/label for the output from HDBSCAN.
#' @param out_dir The output directory.
#' @param matrix_file The matrix file of enrichment scores for the two samples being compared, generated from the pre_clus() function.
#' @param minpts The minimum number of data points for a cluster size. Set it to the smallest size grouping that you wish to consider a cluster
#' @param minsamps This provides a measure of how conservative you want your clustering to be. The larger the value of min_samps you provide, the more conservative the clustering – more points will be declared as noise, and clusters will be restricted to progressively more dense areas.
#' @param cores The number of parallel cores to be used.
#'
#' @return A text (.txt) output designating each bin from the matrix file as belonging to a specific cluster or not.
#' @export
#'
#' @include pre_clust.R
#'
#' @example inst/examples/example_generate_clust.R
generate_clust <- function(output_file_name,
                 out_dir,
                 matrix_file,
                 minpts,
                 minsamps,
                 cores) {
  suppressWarnings({
  install_HDBSCAN <-
    function(...,
             envname = "r-HDBSCAN",
             new_env = identical(envname, "r-HDBSCAN")) {

      if(new_env && virtualenv_exists(envname))
        virtualenv_remove(envname)
      #reticulate::install_python(version = '3.9.19')
      reticulate::py_install(packages = "hdbscan", envname = envname,python_version="3.9",...)

      # .onLoad <- function(...) {
      #   use_virtualenv("r-HDBSCAN", required = FALSE)
      # .onLoad()
      #}
    }
  # install isolated local environment
  install_HDBSCAN()
  # requires that the isolated local environment be used
  reticulate::use_virtualenv(virtualenv = "r-HDBSCAN", required = TRUE)
  # gives user info on which environment and python is being used
  print("Please see below for more information regarding the isolated virtual environment being used for this package.")
  print(reticulate::py_config())
  # source the python script
  reticulate::source_python(file=system.file("python", "clus.py", package = "ChIPbinner"))
  # reticulate::source_python("~/Documents/ChIPbinner/inst/python/clus.py")
  # file parameters
  output_file_name <- paste0(output_file_name)
  # output directory
  setwd(paste0(out_dir))
  # input file
  mat_file <- paste0(matrix_file)
  # algorithm parameters
  minpts <- as.integer(paste0(minpts))
  minsamps <- as.integer(paste0(minsamps))
  # number of parallel cores to use
  cores <- as.integer(paste0(cores))
  # out_dir
  hdbscan_clustering(file_name = output_file_name,
    mat_file = mat_file,
    minpts = minpts,
    minsamps = minsamps,
    cores = cores
  )
  print("Density-based clusters generated using HDBSCAN!")
  })
}
