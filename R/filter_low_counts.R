#!/usr/bin/env Rscript
#' Title
#' @title Filter bins depleted in signal across all samples
#' @description Filters bins with raw read count consistently lower than a specified value across all samples
#' @param out_dir a character string specifying the output directory for the sizeFactor-normalized BED files.
#' @param genome_assembly a character string specifying the genome assembly. Allowed values include "hg38" or "mm10".
#' @param sample_bedfiles a vector specifying the BED file(s) for the samples with raw counts.
#' @param cutoff a numeric value specifying the raw read count cutoff for bins. Defaults to 100.
#' @return BED files with filtered counts
#' @export
#' @examples
#' \dontrun{
#' filter_low_counts(out_dir, hg38, sample_bedfiles = c("NSD1KO_rep1.bed", "NSD1KO_rep2.bed","WT_rep1.bed", "WT_rep2.bed"),cutoff=100)
#' }
filter_low_counts <- function(out_dir,
                              genome_assembly,
                              sample_bedfiles,
                              cutoff=100) {
  suppressWarnings({
    # directory parameters
    out_dir <- paste0(out_dir)
    if (genome_assembly == "hg38") {
      bl <- hg38_bl
      chrom_sizes <- hg38_chrom_sizes
    } else if (genome_assembly == "mm10") {
      bl <- mm10_bl
      chrom_sizes <- mm10_chrom_sizes
    } else {
      stop("Allowed values are hg38 or mm10.")
    }

    read_and_agg_samps <- function(samps) {
      # Read all samples and store as GRanges
      samp_list_gr <- purrr::map(samps, function(samp) {
        samp_raw_bed <- rtracklayer::import.bed(samp)  # Read bed file
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

    agg_all_samps <- read_and_agg_samps(samps = sample_bedfiles)

    # filter counts
    filtered_counts <- agg_all_samps %>%
      dplyr::mutate(max = apply(., 1, max)) %>%
      dplyr::filter(max > as.numeric(cutoff))

    filtered_counts$max <- NULL

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

    # Loop through each column (variable) in filtered_counts
    for (samp in colnames(filtered_counts)) {
      # Get the genomic coordinates for the current sample
      samp_bed <- get_genomic_coordinates(samps = sample_bedfiles) %>% as.data.frame()

      # Assign the corresponding normalized counts to the 'name' column of the BED data frame
      # samp_bed$score <- filtered_counts[[samp]]

      samp_bed$name <- NULL

      samp_df <- samp_bed %>%
        as.data.frame() %>%
        dplyr::select(c(1:3))

      samp_df$score <- filtered_counts[[samp]]

      # Export the data to a BED file
      # rtracklayer::export.bed(object = samp_bed, con = sprintf("%s/%s_sizeFactor_norm.bed", out_dir, samp))
      readr::write_tsv(samp_df, file = sprintf("%s/%s_filt.bed", out_dir, samp), col_names = FALSE, quote = "none")

    }

  })
}
