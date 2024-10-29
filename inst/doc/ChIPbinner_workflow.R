## ----include = FALSE-----------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------------------------------------------------------------
library(ChIPbinner)

## ----include=FALSE-------------------------------------------------------------------------------------------------------------
library(formatR)
library(knitr)

## ----eval=FALSE,tidy=TRUE, tidy.opts=list(width.cutoff=60)---------------------------------------------------------------------
#  # generate normalized bigwig for WT sample
#  norm_bw(out_dir = "output_directory",
#          genome_assembly = "hg38",
#          immunoprecipitated_binned_file = "Cal27.WT.H3K36me2.10kb.bed",
#          use_input = TRUE,
#          input_binned_file = "Cal27.WT_input.H3K36me2.10kb.bed",
#          raw_count_cutoff = 0,
#          pseudocount = 1e-15,
#          scaling_factor = 0.450328805)
#  # generate normalized bigwig for NSD1-KO sample
#  norm_bw(out_dir = "output_directory",
#            genome_assembly = "hg38",
#            immunoprecipitated_binned_file = "Cal27.NSD1_KO.H3K36me2.10kb.bed",
#            use_input = TRUE,
#            input_binned_file = "Cal27.NSD1_KO_input.H3K36me2.10kb.bed",
#            raw_count_cutoff = 0,
#            pseudocount = 1e-15,
#            scaling_factor = 0.192272095)

## ----eval=FALSE, tidy=TRUE, tidy.opts=list(width.cutoff=60)--------------------------------------------------------------------
#  genic_intergenic_scatterplot(
#    out_dir = "output_directory",
#    genome_assembly = "hg38",
#    cell_line = "Cal27",
#    wildtype_samp_norm_bw = "Cal27.WT.H3K36me2.10kb.norm.bw",
#    treated_samp_norm_bw = "Cal27.NSD1_KO.H3K36me2.10kb.norm.bw",
#    are_R_objects = FALSE,
#    histone_mark = "H3K36me2",
#    output_filename = "Cal27.WT_NSD1KO.H3K36me2.10kb",
#    title_of_plot = "Cal27 ChIP-seq H3K36me2",
#    xaxis_label = "WT",
#    yaxis_label = "NSD1_KO",
#    max_x = 1,
#    max_y = 1,
#    min_x = -5,
#    min_y = -5,
#    pow = 1.25,
#    show_scales = FALSE,
#    show_legend = TRUE,
#    legend_pos = "left"
#  )

## ----eval=FALSE----------------------------------------------------------------------------------------------------------------
#    plot_PCA(
#      out_dir = "outdir",
#      output_filename = "PCA_plot",
#      colors = c("blue","blue3","red","red4","forestgreen","green1"),
#      sample_labels = c("samp_1", "samp_2", "samp_3", "samp_4", "samp_5", "samp_6"),
#      plot_title = "",
#      plot_height = 8,
#      plot_width= 12,
#      "samp1_10kb.norm.bw",
#      "samp2_10kbb.norm.bw",
#      "samp3_10kbb.norm.bw",
#      "samp4_10kbb.norm.bw",
#      "samp5_10kbb.norm.bw",
#      "samp6_10kbb.norm.bw",
#    )

## ----eval=FALSE----------------------------------------------------------------------------------------------------------------
#  plot_correlation(
#      out_dir = "outdir",
#      output_filename = "correlation_plot",
#      sample_labels = c("samp_1", "samp_2", "samp_3", "samp_4", "samp_5", "samp_6"),
#      colors = c("blue","blue3","red","red4","forestgreen","green1"),
#      plot_height = 4,
#      plot_width = 4,
#      "samp1_10kb.norm.bw",
#      "samp2_10kbb.norm.bw",
#      "samp3_10kbb.norm.bw",
#      "samp4_10kbb.norm.bw",
#      "samp5_10kbb.norm.bw",
#      "samp6_10kbb.norm.bw"
#    )

## ----eval=FALSE----------------------------------------------------------------------------------------------------------------
#  merge_norm_bw(rep1 = "WT_rep1_10kb_norm.bw",
#                rep2 = "WT_rep2_10kb_norm.bw",
#                merged_samples_label = "merged_WT",
#                out_dir = "outdir")

## ----eval=FALSE----------------------------------------------------------------------------------------------------------------
#  pre_clust(
#    out_dir = "output_directory",
#    treated_samp_norm_bw = "Cal27.WT.H3K36me2.10kb.norm.bw",
#    wildtype_samp_norm_bw = "Cal27.NSD1_KO.H3K36me2.10kb.norm.bw",
#    output_filename = "Cal27.WT.NSD1_KO.H3K36me2.10kb",
#    are_R_objects = FALSE
#  )

## ----eval=FALSE----------------------------------------------------------------------------------------------------------------
#  # using the complete matrix
#  generate_clust(
#    output_file_name = "Cal27.WT.NSD1_KO.H3K36me2.10kb",
#    out_dir = "output_directory",
#    matrix_file = "Cal27.WT.NSD1_KO.H3K36me2.10kb.mat.csv",
#    minpts = 5000,
#    minsamps = 5000,
#    cores = 6
#  )

## ----eval=FALSE----------------------------------------------------------------------------------------------------------------
#  annotate_clust(number_of_clusters = 3,
#                 matrix_file = "Cal27.WT.NSD1_KO.H3K36me2.10kb_mat.csv",
#                 pooled_bed_file = "Cal27.WT.NSD1_KO.H3K36me2.10kb_pooled.bed",
#                 hdbscan_output_file = "clus.Cal27.WT.NSD1_KO.H3K36me2.10kb.5000.5000.txt",
#                 output_filename = "Cal27.WT.NSD1_KO.H3K36me2.10kb",
#                 out_dir="output_directory")

## ----eval=FALSE----------------------------------------------------------------------------------------------------------------
#  density_based_scatterplot(
#    out_dir = "output_directory",
#    genome_assembly = "hg38",
#    are_R_objects = FALSE,
#    output_filename = "Cal27.WT.NSD1_KO.H3K36me2.10kb",
#    wildtype_samp_norm_bw ="Cal27_WT_H3K36me2_10kb.norm.bw",
#    treated_samp_norm_bw = "Cal27_NSD1_KO_H3K36me2_10kb.norm.bw",
#    cell_line = "Cal27",
#    histone_mark = "H3K36me2",
#    annotated_clusters ="Cal27.WT.NSD1_KO.H3K36me2.10kb_annotated_clusters.rda",
#    number_of_clusters = 3,
#    title_of_plot = "H3K36me2",
#    pow = 1.1,
#    show_legend = TRUE,
#    min = -5,
#    max = 2,
#    bin_size = 50,
#    show_scales = FALSE,
#    xaxis_label = "WT",
#    yaxis_label = "NSD1_KO",
#    height_of_figure = 6,
#    width_of_figure = 15
#  )

## ----eval=FALSE----------------------------------------------------------------------------------------------------------------
#  enrich_clust(
#    genome_assembly = "hg38",
#    annotated_clusters = "Cal27.WT.NSD1_KO.H3K36me2.10kb_annotated_clusters.rda",
#    query_cluster = "B",
#    pooled_bed_file = "Cal27.WT.NSD1_KO.H3K36me2.10kb_pooled.bed",
#    functional_db = "hg38_ensemblDB.rda",
#    region = "genome_wide",
#    cores = 1,
#    n_elements = 7,
#    cutoff_for_overlap = 1000,
#    file_plot_name = "ensembl.genome_wide.enrichment_depletion",
#    output_table_name = "ensembl.genome_wide.enrichment_depletion",
#    width_of_plot = 7,
#    height_of_plot = 3,
#    out_dir = "output_directory"
#    )

