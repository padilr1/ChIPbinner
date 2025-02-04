# test_that("Running differential bin analysis works", {
#   clusterSpecific_bin_stats(
#     out_dir = "~/Documents/ChIPbinner_manuscript/analysis/simulating_downsampled_dataset/sim_output/H3K36me2/differential_output",
#     genome_assembly = "hg38",
#     treated_sample_bedfiles = c("~/Documents/ChIPbinner_manuscript/analysis/simulating_downsampled_dataset/sim_output/H3K36me2/Cal27_PA_H3K36me2_rep1_downsampled_noisy_sorted.10kb.bed", "~/Documents/ChIPbinner_manuscript/analysis/simulating_downsampled_dataset/sim_output/H3K36me2/Cal27_PA_H3K36me2_rep2_downsampled_noisy_sorted.10kb.bed"),
#     wildtype_sample_bedfiles = c("~/Documents/ChIPbinner_manuscript/analysis/simulating_downsampled_dataset/sim_output/H3K36me2/Cal27_P8_k36me2.10kb.bed", "~/Documents/ChIPbinner_manuscript/analysis/simulating_downsampled_dataset/sim_output/H3K36me2/Cal27_unedited-c6_K36me2.10kb.bed"),
#     treated_condition_label = "Downsampled",
#     wildtype_condition_label = "WT",
#     annotated_clusters = "~/Documents/ChIPbinner_manuscript/analysis/simulating_downsampled_dataset/sim_output/H3K36me2/clus.depthNorm_Cal27_PA_H3K36me2_rep1_downsampled_noisy_sorted.10kb_Cal27_P8_k36me2.10kb.15000.1.annotated_clusters.rda",
#     output_filename = "test_Cal27_WT_WT_downsampled_H3K36me2",
#     return_results_for_all_bins = FALSE
#   )
#   expect_snapshot_output(x = "Stats per bin generated!")
# })
