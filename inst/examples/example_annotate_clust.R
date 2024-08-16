annotate_clust(number_of_clusters = 3,
               matrix_file = system.file("extdata", "downsampled.Cal27.WT.NSD1_KO.H3K36me2.10kb_mat.csv", package = "ChIPbinner"),
               pooled_bed_file = system.file("extdata", "downsampled.Cal27.WT.NSD1_KO.H3K36me2.10kb_pooled.bed.gz", package = "ChIPbinner"),
               hdbscan_output_file= system.file("extdata", "clus.downsampled.Cal27.WT.NSD1_KO.H3K36me2.10kb.100.100.txt", package = "ChIPbinner"),
               output_filename = "downsampled",
               out_dir=system.file("extdata",package = "ChIPbinner"))
