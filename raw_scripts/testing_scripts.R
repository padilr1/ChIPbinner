# for testing

# when tweaking the function, will need to re-load it
devtools::load_all()
devtools::test()

devtools::load_all()
testthat::test_file("tests/testthat/test-genic_intergenic_scatterplot.R")
