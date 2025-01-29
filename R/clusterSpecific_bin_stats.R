#!/usr/bin/env Rscript
#' Title
#' @title Extract results per bin within each cluster from a DESeq analysis
#'
#' @description Performs differential bin analysis on raw counts using DESeq2. Log2 fold changes are shrunk using the "apeglm" method. Results include base means, log2 fold changes, standard errors, p-values, adjusted p-values, genomic coordinates, and cluster assignments for each bin.
#'
#' @param out_dir a character string specifying the output directory for the sizeFactor-normalized BED files.
#' @param genome_assembly a character string specifying the genome assembly. Allowed values include "hg38" or "mm10".
#' @param treated_sample_bedfiles a vector specifying the BED file(s) for the treated sample(s) with raw counts.
#' @param wildtype_sample_bedfiles a vector specifying the BED file(s) for the wildtype sample(s) with raw counts.
#' @param treated_condition_label a character string specifying the condition for the treated sample(s).
#' @param wildtype_condition_label a character string specifying the condition for the wildtype sample(s).
#' @param annotated_clusters a character string specifying the R object containing the annotated clusters generated using annotate_clust(). If no clusters are provided, then the function proceeds without associating any bins to a cluster.
#' @param output_filename a character string specifying the file name for the resulting table to be saved on disk.
#' @param return_results_for_all_bins a logical indicating whether to also return a table with all results, including bins not associated with any cluster. Defaults to FALSE.
#' @return a csv file
#' @export
#'
#'
#' @examples
#' \dontrun{
#' clusterSpecific_bin_stats(out_dir, mm10, treated_sample_bedfiles = c("NSD1KO_rep1.bed", "NSD1KO_rep2.bed"), wildtype_sample_bedfiles = c("WT_rep1.bed", "WT_rep2.bed"), treated_condition_label = "NSD1KO", wildtype_condition_label = "WT", annotated_clusters = "input_directory/annotated_clusters.rda", output_filename = "results_per_bin_per_cluster")
#' }
clusterSpecific_bin_stats <- function(out_dir,
                                      genome_assembly,
                                      treated_sample_bedfiles,
                                      wildtype_sample_bedfiles,
                                      treated_condition_label,
                                      wildtype_condition_label,
                                      annotated_clusters = NULL,
                                      output_filename,
                                      return_results_for_all_bins = FALSE) {
  suppressWarnings({
    # directory parameters
    out_dir <- paste0(out_dir)
    if (genome_assembly == "hg38") {
      bl <- hg38_bl
      chrom_sizes <- hg38_chrom_sizes
    } else if (genome_assembly == "mm10") {
      bl <- mm10_bl
      chrom_sizes <- mm10_chrom_sizes
    }

    read_and_agg_samps <- function(samps) {
      # Use purrr::map to process each sample and return a list of data frames
      samp_list <- purrr::map(samps, function(samp) {
        samp_label <- basename(tools::file_path_sans_ext(samp)) # Get sample label from file name

        samp_raw_bed <- rtracklayer::import.bed(samp)

        samp_bed_gr <- samp_raw_bed[!IRanges::overlapsAny(samp_raw_bed, bl)]

        samp_bed <- samp_bed_gr %>%
          as.data.frame() # Import BED file as data frame

        # Rename the column for the sample
        colnames(samp_bed)[6] <- paste0(samp_label)

        # Select the relevant column (the one with the sample data)
        samp_mat <- samp_bed %>%
          dplyr::select(c(paste0(samp_label))) %>%
          dplyr::mutate(across(everything(), as.numeric))

        return(samp_mat)
      })

      # Combine all the sample matrices into one data frame
      agg_samp <- dplyr::bind_cols(samp_list)

      return(agg_samp)
    }
    agg_treated_samps <- read_and_agg_samps(samps = treated_sample_bedfiles)
    agg_wildtype_samps <- read_and_agg_samps(samps = wildtype_sample_bedfiles)
    metadata_treated_samps <- data.frame(kind = colnames(agg_treated_samps)) %>% dplyr::mutate(cond = treated_condition_label)
    metadata_wildtype_samps <- data.frame(kind = colnames(agg_wildtype_samps)) %>% dplyr::mutate(cond = wildtype_condition_label)
    if (ncol(agg_treated_samps) < 2 || ncol(agg_wildtype_samps) < 2) {
      stop("At least two replicates per condition are needed.")
    }
    if (nrow(agg_treated_samps) != nrow(agg_wildtype_samps)) {
      stop("The number of genomic regions differs between wildtype and treated samples. The same number of genomic regions should be found in both wildtype and treated samples.")
    } else if (nrow(agg_treated_samps) == nrow(agg_wildtype_samps)) {
      agg_metadata <- rbind(metadata_wildtype_samps, metadata_treated_samps) %>% tibble::column_to_rownames(., "kind")
      agg_counts <- cbind(agg_wildtype_samps, agg_treated_samps) %>% as.matrix()
    }

    norm_method <- "DESeq2"
    if (norm_method == "DESeq2") {
      ### get genomic coordinates ###
      get_genomic_coordinates <- function(samps) {
        # Import the first sample as a GRanges object
        first_sample_bed <- rtracklayer::import.bed(samps[1])

        # Filter out regions that overlap with 'bl' (assuming 'bl' is predefined)
        first_sample_bed_filtered <- first_sample_bed[!IRanges::overlapsAny(first_sample_bed, bl)]

        # Use purrr::reduce to iteratively apply subsetByOverlaps to find overlapping coordinates across all samples
        overlapping_coords <- purrr::reduce(samps[-1], function(accumulated_gr, samp) {
          samp_bed <- rtracklayer::import.bed(samp)

          # Filter out overlapping regions with 'bl' for the current sample
          samp_bed_filtered <- samp_bed[!IRanges::overlapsAny(samp_bed, bl)]

          # Find overlaps between the accumulated regions and the current sample
          overlap_result <- IRanges::subsetByOverlaps(accumulated_gr, samp_bed_filtered)

          # Return the overlapping regions so far
          return(overlap_result)
        }, .init = first_sample_bed_filtered)

        return(overlapping_coords)
      }

      samp_bed <- get_genomic_coordinates(samps = c(treated_sample_bedfiles, wildtype_sample_bedfiles)) %>%
        as.data.frame() %>%
        dplyr::select(c("seqnames", "start", "end")) %>%
        dplyr::mutate(genomic_coord = sprintf("%s.%d.%d", .$seqnames, .$start, .$end)) %>%
        dplyr::select("genomic_coord")

      agg_counts_with_coordinates <- cbind(samp_bed, agg_counts) %>% tibble::column_to_rownames("genomic_coord")

      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = agg_counts_with_coordinates,
        colData = agg_metadata,
        design = ~cond
      )
      dds <- DESeq2::DESeq(dds)

      dds$condition <- stats::relevel(dds$cond, paste0(wildtype_condition_label))

      dds <- DESeq2::nbinomWaldTest(dds)

      resLFC <- DESeq2::lfcShrink(dds = dds, coef = 2, type = "apeglm") %>%
        as.data.frame() %>%
        tibble::rownames_to_column("genomic_coord") %>%
        na.omit()

      resLFC[c("chr", "start", "end")] <- stringr::str_split_fixed(resLFC$genomic_coord, "[.]", 3)

      resLFC$genomic_coord <- NULL

      resLFC$start <- as.integer(resLFC$start)

      adjust_last_digit <- function(number) {
        last_digit <- number %% 10
        if (last_digit == 1) {
          return(number)
        } else {
          return(as.numeric(number - last_digit + 1))
        }
      }

      resLFC$start <- sapply(resLFC$start, adjust_last_digit)

      resLFC$end <- as.integer(resLFC$end)

      resLFC$genomic_coord <- sprintf("%s.%d.%d", resLFC$chr, resLFC$start, resLFC$end)

      if (is.null(annotated_clusters)) {
        readr::write_csv(resLFC, file = sprintf("%s/%s_stats_all_bins_no_clustering.csv", out_dir, output_filename), col_names = TRUE, quote = "none")
        stop("The results table for all bins without clustering information has been generated.")
      }

      if (return_results_for_all_bins == TRUE) {
        readr::write_csv(resLFC, file = sprintf("%s/%s_stats_all_bins_no_clustering.csv", out_dir, output_filename), col_names = TRUE, quote = "none")
      }

      # remove redundant columns
      resLFC$chr <- NULL
      resLFC$start <- NULL
      resLFC$end <- NULL

      # loading function
      loadRData <- function(fileName) {
        # loads an RData file, and returns it
        load(fileName)
        get(ls()[ls() != "fileName"])
      }

      cons <- loadRData(paste0(annotated_clusters))

      get_dataframe_from_cluster <- function(x) {
        df <- x %>% as.data.frame()
        df$start <- sapply(df$start, adjust_last_digit)
        df$genomic_coord <- sprintf("%s.%d.%d", df$seqnames, df$start, df$end)
        df <- df %>%
          dplyr::left_join(., resLFC, by = "genomic_coord") %>%
          na.omit()
        return(df)
      }
      cons_dfs <- lapply(cons, get_dataframe_from_cluster)

      # Function to modify each data frame in the list and cluster name
      add_cluster_name <- function(df, df_name) {
        df %>%
          mutate(cluster = df_name)
      }

      # Apply the function to each data frame in the list
      modified_dfs <- purrr::map2(cons_dfs, names(cons_dfs), add_cluster_name)

      # aggregate clusters
      agg_cluster <- dplyr::bind_rows(modified_dfs)

      # write csv
      readr::write_csv(agg_cluster, file = sprintf("%s/%s_stats_per_cluster.csv", out_dir, output_filename), col_names = TRUE, quote = "none")
    }
    # print output message
    print("Stats per bin generated!")
  })
}
