# test_that("Correlation plot works", {
#   # plot_correlation(
#   #   out_dir = testthat::test_path("testdata"),
#   #   output_filename = "Cal27_KO",
#   #   sample_labels = c("PA", "PA", "NSD1KO", "NSD1KO", "NSD2KO", "NSD2KO"),
#   #   correlation_method = "spearman",
#   #   colors = c("blue", "blue", "red", "red", "forestgreen", "forestgreen"),
#   #   # key_title = "Pearson's R",
#   #   # col_label_angle = 45,
#   #   # row_label_angle = 45,
#   #   # row_adj_label = c(0,0),
#   #   # col_adj_label = c(1,1),
#   #   "~/Documents/HNSCC_K36me2/data/norm.bw/unscaled_norm_bw/Cal27/Cal27.P8.H3K36me2.10kb.norm.bw",
#   #   "~/Documents/HNSCC_K36me2/data/norm.bw/unscaled_norm_bw/Cal27/Cal27.unedited_c6.H3K36me2.10kb.norm.bw",
#   #   "~/Documents/HNSCC_K36me2/data/norm.bw/unscaled_norm_bw/Cal27/Cal27.NSD1KO_1.H3K36me2.10kb.norm.bw",
#   #   "~/Documents/HNSCC_K36me2/data/norm.bw/unscaled_norm_bw/Cal27/Cal27.NSD1KO_17.H3K36me2.10kb.norm.bw",
#   #   "~/Documents/HNSCC_K36me2/data/norm.bw/unscaled_norm_bw/Cal27/Cal27.NSD2KO_c2.H3K36me2.10kb.norm.bw",
#   #   "~/Documents/HNSCC_K36me2/data/norm.bw/unscaled_norm_bw/Cal27/Cal27.NSD2KO_c3.H3K36me2.10kb.norm.bw"
#   expect_snapshot_output(x = "Correlation plot generated!")
# })
