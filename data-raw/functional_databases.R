## code to prepare `functional_databases` dataset goes here

usethis::use_data(functional_databases, overwrite = TRUE)

library(LOLA)

hg38_ensemblDB <- LOLA::loadRegionDB("~/Documents/ChIPbinner/inst/extdata/regionDB/hg38", collections = "ensembl")
hg38_ccreDB <- LOLA::loadRegionDB("~/Documents/ChIPbinner/inst/extdata/regionDB/hg38", collections = "ccre")
hg38_repeatsDB <- LOLA::loadRegionDB("~/Documents/ChIPbinner/inst/extdata/regionDB/hg38", collections = "repeats")

downsampled_hg38_ensembl <-LOLA::loadRegionDB("~/Documents/ChIPbinner/inst/extdata/regionDB/hg38", collections = "ensembl")
usethis::use_data(downsampled_hg38_ensembl,internal=FALSE,compress="xz",overwrite = TRUE)

downsampled_hg38_ccreDB <- LOLA::loadRegionDB("~/Documents/ChIPbinner/inst/extdata/regionDB/hg38", collections = "ccre")
usethis::use_data(downsampled_hg38_ccreDB,internal=FALSE,compress="xz",overwrite = TRUE)

usethis::use_data(hg38_ensemblDB,internal=FALSE,compress="xz",overwrite = TRUE)
usethis::use_data(hg38_ccreDB,internal=FALSE,compress="xz",overwrite = TRUE)
usethis::use_data(hg38_repeatsDB,internal=FALSE,compress="xz",overwrite = TRUE)

mm10_ensemblDB <- LOLA::loadRegionDB("~/Documents/ChIPbinner/inst/extdata/regionDB/mm10", collections = "ensembl")
mm10_ccreDB <- LOLA::loadRegionDB("~/Documents/ChIPbinner/inst/extdata/regionDB/mm10", collections = "ccre")
mm10_repeatsDB <- LOLA::loadRegionDB("~/Documents/ChIPbinner/inst/extdata/regionDB/mm10", collections = "repeats")

usethis::use_data(mm10_ensemblDB,internal=FALSE,compress="xz",overwrite = TRUE)
usethis::use_data(mm10_ccreDB,internal=FALSE,compress="xz",overwrite = TRUE)
usethis::use_data(mm10_repeatsDB,internal=FALSE,compress="xz",overwrite = TRUE)
