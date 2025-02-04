#!/usr/bin/env Rscript
#' Title
#' @title Extract results for statistical testing per bin between two groups within each cluster
#'
#' @description Performs differential bin analysis on normalized and/or scaled counts using reproducibility-optimized test statistics (ROTS) between two groups. Results include log fold changes, adjusted p-values, genomic coordinates, and cluster assignments for each bin.
#'
#' @param out_dir a character string specifying the output directory.
#' @param genome_assembly a character string specifying the genome assembly. Allowed values include "hg38" or "mm10".
#' @param treated_sample_bigWigFiles a vector specifying the normalized bigWig file(s) for the treated sample(s).
#' @param wildtype_sample_bigWigFiles a vector specifying the the normalized bigWig file(s) for the wildtype sample(s).
#' @param treated_condition_label a character string specifying the condition for the treated sample(s).
#' @param wildtype_condition_label a character string specifying the condition for the wildtype sample(s).
#' @param annotated_clusters a character string specifying the R object containing the annotated clusters generated using annotate_clust(). If no clusters are provided, then the function proceeds without associating any bins to a cluster.
#' @param bootstrap_value an integer specifying the number of bootstrap and permutation resamplings to estimate the null distribution of the test statistic (default 1000). Increasing B can improve the precision of your results, but at the expense of greater computational run-time.
#' @param K_value an integer indicating the largest top list size considered. It is recommended that the K value should be considerably higher than the number of features expected to be differentially expressed.
#' @param output_filename a character string specifying the file name for the resulting table to be saved on disk. An output filename will be automatically generated if none is specified.
#' @param return_results_for_all_bins a logical indicating whether to also return a table with all results, including bins not associated with any cluster. Defaults to FALSE.
#' @param precomputed_ROTS a character string specifying the R object containing the precomputed ROTS object. This parameter should only be used if ROTS has already been executed, and the goal is to match the resulting bins to their assigned clusters.
#' @references Suomi T, Seyednasrollah F, Jaakkola MK, Faux T, Elo LL. ROTS: An R package for reproducibility-optimized statistical testing. PLoS Comput Biol 2017; 13: e1005562.
#' @return a csv file of statistical comparison results for each bin assigned to their respective clusters.
#' @export
#'
#'
#' @examples
#' \dontrun{
#' clusterSpecific_bin_stats(out_dir, hg38, treated_sample_bigWigFiles,= c("NSD1KO_rep1.bw", "NSD1KO_rep2.bw"), wildtype_sample_bigWigFiles = c("WT_rep1.bw", "WT_rep2.bw"), treated_condition_label = "NSD1KO", wildtype_condition_label = "WT",bootstrap_value=1000,K_value = 5000, annotated_clusters = "input_directory/annotated_clusters.rda")
#' }
clusterSpecific_bin_stats <- function(out_dir,
                                      genome_assembly,
                                      treated_sample_bigWigFiles,
                                      wildtype_sample_bigWigFiles,
                                      treated_condition_label,
                                      wildtype_condition_label,
                                      bootstrap_value=1000,
                                      K_value,
                                      annotated_clusters = NULL,
                                      output_filename = NULL,
                                      return_results_for_all_bins = FALSE,
                                      precomputed_ROTS = NULL) {
  suppressWarnings({
    # loading function
    loadRData <- function(fileName) {
      # loads an RData file, and returns it
      load(fileName)
      get(ls()[ls() != "fileName"])
    }

    # directory parameters
    if (genome_assembly == "hg38") {
      bl <- hg38_bl
      chrom_sizes <- hg38_chrom_sizes
    } else if (genome_assembly == "mm10") {
      bl <- mm10_bl
      chrom_sizes <- mm10_chrom_sizes
    }

    read_and_agg_samps <- function(samps) {
      # Read all samples and store as GRanges
      samp_list_gr <- purrr::map(samps, function(samp) {
        samp_raw_bed <- rtracklayer::import.bw(samp)  # Read bigWig file
        samp_bed_gr <- samp_raw_bed[!IRanges::overlapsAny(samp_raw_bed, bl)]  # Remove blacklist regions
        return(samp_bed_gr)
      })

      # **Find common overlapping regions**
      common_regions <- base::Reduce(intersect, samp_list_gr)

      # Ensure all samples only contain these regions
      samp_list_filtered <- purrr::map(samp_list_gr, function(samp_bed_gr) {
        IRanges::subsetByOverlaps(samp_bed_gr, common_regions)  # Retain only overlapping regions
      })

      # Convert filtered GRanges objects to data frames
      samp_list <- purrr::map2(samp_list_filtered, samps, function(samp_bed_gr, samp) {
        samp_label <- basename(tools::file_path_sans_ext(samp))
        samp_bed <- as.data.frame(samp_bed_gr)

        # Rename the value column with sample label
        colnames(samp_bed)[6] <- paste0(samp_label)  # Assuming signal values are in column 6

        # Extract relevant sample column and ensure numeric format
        samp_mat <- samp_bed %>%
          dplyr::select(paste0(samp_label)) %>%
          dplyr::mutate(across(everything(), as.numeric))

        return(samp_mat)
      })

      # Combine all the sample matrices into one data frame
      agg_samp <- dplyr::bind_cols(samp_list)

      return(agg_samp)
    }


    agg_all_samps <- read_and_agg_samps(samps = c(treated_sample_bigWigFiles,wildtype_sample_bigWigFiles))

    treated_sampleNames <- lapply(treated_sample_bigWigFiles,function(x){basename(tools::file_path_sans_ext(x))}) %>% as.character()
    wildtype_sampleNames <- lapply(wildtype_sample_bigWigFiles,function(x){basename(tools::file_path_sans_ext(x))}) %>% as.character()

    agg_treated_samps  <- agg_all_samps %>% dplyr::select(dplyr::all_of(treated_sampleNames))
    agg_wildtype_samps <- agg_all_samps %>% dplyr::select(dplyr::all_of(wildtype_sampleNames))

    metadata_treated_samps <- data.frame(kind = colnames(agg_treated_samps)) %>% dplyr::mutate(cond = treated_condition_label)
    metadata_wildtype_samps <- data.frame(kind = colnames(agg_wildtype_samps)) %>% dplyr::mutate(cond = wildtype_condition_label)
    if (ncol(agg_treated_samps) < 2 || ncol(agg_wildtype_samps) < 2) {
      stop("At least two replicates per condition are needed.")
    } else if (nrow(agg_treated_samps) == nrow(agg_wildtype_samps)) {
      agg_metadata <- rbind(metadata_wildtype_samps, metadata_treated_samps) %>% tibble::column_to_rownames(., "kind")
      agg_counts <- cbind(agg_wildtype_samps, agg_treated_samps) %>% as.matrix()
    }

    norm_method <- "ROTS"
    if (norm_method == "ROTS") {
      ### get genomic coordinates ###
      get_genomic_coordinates <- function(samps) {
        # Import the first sample as a GRanges object
        # first_sample_bed <- rtracklayer::import.bed(samps[1])

        first_sample_bed <- rtracklayer::import.bw(samps[1])

        # Filter out regions that overlap with 'bl' (assuming 'bl' is predefined)
        first_sample_bed_filtered <- first_sample_bed[!IRanges::overlapsAny(first_sample_bed, bl)]

        # Use purrr::reduce to iteratively apply subsetByOverlaps to find overlapping coordinates across all samples
        overlapping_coords <- purrr::reduce(samps[-1], function(accumulated_gr, samp) {
          # samp_bed <- rtracklayer::import.bed(samp)
          samp_bed <- rtracklayer::import.bw(samp)

          # Filter out overlapping regions with 'bl' for the current sample
          samp_bed_filtered <- samp_bed[!IRanges::overlapsAny(samp_bed, bl)]

          # Find overlaps between the accumulated regions and the current sample
          overlap_result <- IRanges::subsetByOverlaps(accumulated_gr, samp_bed_filtered)

          # Return the overlapping regions so far
          return(overlap_result)
        }, .init = first_sample_bed_filtered)

        return(overlapping_coords)
      }

      samp_bw <- get_genomic_coordinates(samps = c(treated_sample_bigWigFiles, wildtype_sample_bigWigFiles)) %>%
        as.data.frame() %>%
        dplyr::select(c("seqnames", "start", "end")) %>%
        dplyr::mutate(genomic_coord = sprintf("%s.%d.%d", .$seqnames, .$start, .$end)) %>%
        dplyr::select("genomic_coord")

      # agg_metadata$cond <- factor(x = agg_metadata$cond,levels = c(wildtype_condition_label,treated_condition_label))

      agg_counts_with_coordinates <- cbind(samp_bw, agg_counts) %>% tibble::column_to_rownames("genomic_coord") %>% as.matrix()

      if (is.null(precomputed_ROTS)) {
        rots.out <- ROTS::ROTS(agg_counts_with_coordinates, groups = agg_metadata$cond,B = bootstrap_value, K = K_value,progress=TRUE,log = TRUE,seed=123)
      } else {
        rots.out <- loadRData(paste0(precomputed_ROTS))
      }

      # get df
      # ROTS_df_preProcessed <- rots.out$data %>% as.data.frame()
      #
      # ROTS_treated_samps  <- ROTS_df_preProcessed %>% dplyr::select(dplyr::all_of(treated_sampleNames))
      # ROTS_wildtype_samps <- ROTS_df_preProcessed %>% dplyr::select(dplyr::all_of(wildtype_sampleNames))
      #
      # ROTS_df <- cbind(ROTS_wildtype_samps,ROTS_treated_samps)

      # Saving ROTS output to a file
      utils::capture.output({
        # Capture the full summary output as a character vector
        full_summary <- utils::capture.output(summary(rots.out, fdr = 0.05))

        # Print the first 10 lines of the summary
        cat(paste(head(full_summary, n = 10), collapse = "\n"))
      }, file = sprintf("%s/ROTS_output_summary.txt", out_dir))

      # resLFC <- bind_cols(list(ROTS_df, rots.out$logfc %>% as.data.frame() %>% `names<-`(c("logFC")), rots.out$FDR %>% as.data.frame() %>% `names<-`(c("FDR")))) %>% rownames_to_column("genomic_coord")

      resLFC <- bind_cols(list(rots.out$logfc %>% as.data.frame() %>% `names<-`(c("logFC")), rots.out$FDR %>% as.data.frame() %>% `names<-`(c("FDR")))) %>% rownames_to_column("genomic_coord")

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

      resLFC <- resLFC[,c("genomic_coord","chr","start","end","logFC","FDR")]

      if (is.null(annotated_clusters)) {
        readr::write_csv(resLFC, file = sprintf("%s/stats_for_all_bins_no_clustering.csv", out_dir), col_names = TRUE, quote = "none")
        stop("The results table for all bins without clustering information has been generated.")
      }

      if (return_results_for_all_bins == TRUE) {
        readr::write_csv(resLFC, file = sprintf("%s/stats_for_all_bins_no_clustering.csv", out_dir), col_names = TRUE, quote = "none")
      }

      # remove redundant columns
      resLFC$chr <- NULL
      resLFC$start <- NULL
      resLFC$end <- NULL

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
      agg_cluster <- dplyr::bind_rows(modified_dfs) %>% dplyr::select(c("genomic_coord","seqnames","start","end","strand","logFC","FDR","cluster"))

      # write csv
      if (is.null(output_filename)) {
        readr::write_csv(agg_cluster, file = sprintf("%s/stats_per_cluster.csv", out_dir), col_names = TRUE, quote = "none")
      } else {
      readr::write_csv(agg_cluster, file = sprintf("%s/%s.csv", out_dir,output_filename), col_names = TRUE, quote = "none")
      }
    }
    # print output message
    print("Stats per bin generated!")
  })
}
