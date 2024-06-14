library(dbscan)
clust <- function(){

}
mat <- fread("")

mat_file <- read_csv("~/Documents/ChIPbinner/tests/testthat/testdata/Cal27.WT.NSD1_KO.H3K36me2.10kb.mat.csv",col_names = F) %>% dplyr::sample_n(10000,replace=FALSE)
write_csv(x = mat_file,file = "~/Documents/ChIPbinner/tests/testthat/testdata/sampled.Cal27.WT.NSD1_KO.H3K36me2.10kb.mat.csv",col_names = FALSE)
set.seed(1)
cl <- dbscan::hdbscan(mat_file,minPts = 1000)
mat_file$clu <- cl$cluster
