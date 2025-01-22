#!/usr/bin/env Rscript
#' Title
#' @title Normalize binned raw counts using normalization factors
#'
#' @description Normalizes binned raw counts across samples using normalization factors, which are calculated using DESeq2's median ratio method or edgeR's trimmed mean of M-values (TMM)
#' @param norm_method a character string specifying the normalization method to use. Allowed values include "DESeq2" or "edgeR".
#' @param out_dir a character string specifying the output directory for the sizeFactor-normalized BED files.
#' @param genome_assembly a character string specifying the genome assembly. Allowed values include "hg38" or "mm10".
#' @param treated_sample_bedfiles a vector specifying the BED file(s) for the treated sample(s) with raw counts.
#' @param wildtype_sample_bedfiles a vector specifying the BED file(s) for the wildtype sample(s) with raw counts.
#' @param treated_condition_label a character string specifying the condition for the treated sample(s).
#' @param wildtype_condition_label a character string specifying the condition for the wildtype sample(s).
#'
#' @return BED files with normalized counts
#' @export
#'
#'
#' @examples
#' \dontrun{
#' apply_normFactors(norm_method="DESeq2",out_dir, mm10, treated_sample_bedfiles = c("NSD1KO_rep1.bed", "NSD1KO_rep2.bed"), wildtype_sample_bedfiles = c("WT_rep1.bed", "WT_rep2.bed"), treated_condition_label = "NSD1KO", wildtype_condition_label = "WT")
#' }
apply_normFactors <- function(norm_method,
                              out_dir,
                              genome_assembly,
                              treated_sample_bedfiles,
                              wildtype_sample_bedfiles,
                              treated_condition_label,
                              wildtype_condition_label) {
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

    if (nrow(agg_treated_samps) != nrow(agg_wildtype_samps)) {
      stop("The number of genomic regions differs between wildtype and treated samples. The same number of genomic regions should be found in both wildtype and treated samples.")
    } else if (nrow(agg_treated_samps) == nrow(agg_wildtype_samps)) {
      agg_metadata <- rbind(metadata_wildtype_samps, metadata_treated_samps) %>% tibble::column_to_rownames(., "kind")
      agg_counts <- cbind(agg_wildtype_samps, agg_treated_samps) %>% as.matrix()
    }

    if (norm_method == "DESeq2") {
    # generate DESeq2 object
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = agg_counts,
      colData = agg_metadata,
      design = ~cond
    )
    # get size factors
    dds <- DESeq2::estimateSizeFactors(dds)

    print(DESeq2::sizeFactors(dds))

    # apply sizeFactor normalization to raw counts
    normalized_counts <- DESeq2::counts(dds, normalized = TRUE) %>% as.data.frame()

    } else if (norm_method == "edgeR") {

      group <- factor(agg_metadata$cond)

      y <- edgeR::DGEList(counts = agg_counts, group = group)

      y <- edgeR::calcNormFactors(y, method = "TMM")

      print(y$samples)

      normalized_counts <- edgeR::cpm(y,normalized.lib.sizes = TRUE) %>% as.data.frame()

    }

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

    ### export bed files for each sample ###

    # Loop through each column (variable) in normalized_counts
    for (samp in colnames(normalized_counts)) {
      # Get the genomic coordinates for the current sample
      samp_bed <- get_genomic_coordinates(samps = c(treated_sample_bedfiles, wildtype_sample_bedfiles)) %>% as.data.frame()

      # Assign the corresponding normalized counts to the 'name' column of the BED data frame
      # samp_bed$score <- normalized_counts[[samp]]

      samp_bed$name <- NULL

      samp_df <- samp_bed %>%
        as.data.frame() %>%
        dplyr::select(c(1:3))

      samp_df$score <- normalized_counts[[samp]]

      # Export the data to a BED file
      # rtracklayer::export.bed(object = samp_bed, con = sprintf("%s/%s_sizeFactor_norm.bed", out_dir, samp))
      readr::write_tsv(samp_df, file = sprintf("%s/%s_sizeFactorNorm.bed", out_dir, samp), col_names = FALSE, quote = "none")
    }
  })
}
