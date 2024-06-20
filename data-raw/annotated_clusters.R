## code to prepare `annotated_clusters` dataset goes here

usethis::use_data(annotated_clusters, overwrite = TRUE)

# put downsampled norm bigwig files in /data
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

downsampled_annotated_clusters <- loadRData("~/Documents/ChIPbinner/tests/testthat/testdata/downsampled.annotated_clusters.rda")
usethis::use_data(downsampled_annotated_clusters,internal=FALSE,compress="xz",overwrite = TRUE)

complete_annotated_clusters <- loadRData("~/Documents/ChIPbinner/tests/testthat/testdata/complete.annotated_clusters.rda")
usethis::use_data(complete_annotated_clusters,internal=FALSE,compress="xz",overwrite = TRUE)
