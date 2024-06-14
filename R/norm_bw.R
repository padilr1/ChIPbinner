#!/usr/bin/env Rscript
#+ message = FALSE, warning = FALSE
#' Title
#' @title Generate normalized bigWig files from binned BED files.
#'
#' @description Use this function to generate the normalized bigWig files needed for the rest of the package workflows.
#'
#' @param out_dir The output directory for normalized bigWig files.
#' @param chromSizes Chromosome sizes for a specific assembly. hg38 and mm10 are included in the package and can be accessed in 'extdata'.
#' @param blacklist Blacklisted regions to be removed. Blacklisted regions for hg38 and mm10 are included in the package and can be accessed in 'extdata'.
#' @param cell_line The cell line of the samples.
#' @param histone_mark The broad histone mark used for the analysis.
#' @param treated_samp_label The treated sample name/label.
#' @param treated_samp_file The binned BED file associated with the treated sample. This function assumes your file is named according to the following format: cell_line.sample_label.histone_mark.window_size_in_kb.bed
#' @param treated_samp_library_size The mapped library size for your treated sample.
#' @param control_label The name/label for your control sample. This could be the input sample (for ChIP-seq) or IgG (for Cut & Run).
#' @param control_file The binned BED file associated with the control sample. This function assumes your file is named according to the following format: cell_line.sample_label.histone_mark.window_size_in_kb.bed
#' @param control_library_size The library size for the control sample.
#' @param use_control This is a Boolean term, whether to use a control sample or not.
#' @param window_size The window size in kilobytes. For example, use 10 for samples binned in 10 kilobyte (kb) bins.
#' @param pseudocount Pseudocount to avoid division by 0. If not defined, defaults to 1e-3.
#' @param raw_count_cutoff Raw read count cutoff to exclude bins depleted in signal across all tracks. For example, if 10 is inputted, this removes raw read count consistently lower than 10. Defaults to 0.
#' @param scaling_factor Optional: quantatitive normalization/scaling of the raw binned signal by this factor. For example, this could be ChIP-Rx values or genome-wide modification percentage values obtained from mass spectrometry values.
#'
#' @return a single normalized bigWig file with binned scores.
#' @export
#'
#' @example inst/examples/example_norm_bw.R
norm_bw <- function(out_dir,
                    chromSizes,
                    blacklist,
                    cell_line,
                    histone_mark,
                    treated_samp_label,
                    treated_samp_file,
                    treated_samp_library_size,
                    control_label = NULL,
                    control_file = NULL,
                    control_library_size = NULL,
                    use_control,
                    window_size,
                    pseudocount = NULL,
                    raw_count_cutoff = NULL,
                    scaling_factor = NULL) {
  # reference files
  ## blacklist
  blacklist <- paste0(blacklist)
  bl <- rtracklayer::import.bed(blacklist)
  ## chrom sizes
  chrom.sizes <- paste0(chromSizes)
  # samples info
  cell_line <- paste0(cell_line)
  mark <- paste0(histone_mark)
  treated_samp_label <- paste0(treated_samp_label) # watch spaces in names
  treated_samp_file <- rtracklayer::import.bed(paste0(treated_samp_file))
  # other parameters
  ## pseudocount to avoid division by 0
  # pseudocount <- as.numeric(paste0(pseudocount))
  if (is.null(pseudocount)) {
    pseudocount <- as.numeric(1e-3)
  } else {
    pseudocount <- as.numeric(paste0(pseudocount))
  }
  ## raw count cut-off to remove bins with low counts across samples
  # raw_count_cutoff <- as.numeric(paste0(raw_count_cutoff))
  if (is.null(raw_count_cutoff)) {
    raw_count_cutoff <- as.numeric(0)
  } else {
    raw_count_cutoff <- as.numeric(paste0(raw_count_cutoff))
  }
  ## window size of the bins
  window_size <- paste0(".", window_size, "kb.")
  ## scaling factor
  if (is.null(scaling_factor)) {
    to_scale <- 1
  } else {
    to_scale <- as.numeric(paste0(scaling_factor))
  }
  # if using input
  if (is.null(control_file)) {
    print("Not using input")
    treated_samp_library_size <- as.numeric(paste0(treated_samp_library_size))
    print(paste0("Treated sample sequencing depth=", treated_samp_library_size))
    print(paste0("Treated sample=", treated_samp_label))
  } else {
    treated_samp_library_size <- as.numeric(paste0(treated_samp_library_size))
    control_library_size <- as.numeric(paste0(control_library_size))
    print(paste0("Treated sample=", treated_samp_label))
    print(paste0("Treated sample sequencing depth=", treated_samp_library_size))
    print(paste0("Control sample=", control_label))
    print(paste0("Control sample sequencing depth=", control_library_size))
    control_file <- import.bed(control_file)
    d <- list(treated_samp_file, control_file)
    names(d) <- c(treated_samp_label, control_label)
  }
  ### deriving cut-off for raw counts ###
  # read in chrom sizes
  keep <- fread(chrom.sizes) %>%
    setNames(c("chr", "seqlength"))
  gn <- keep %>%
    {
      GenomeInfoDb::Seqinfo(.$chr, .$seqlength)
    }
  # loop through each file and keep only signal for each 10kb bin
  raw <- lapply(d, function(x) {
    inp <- as.data.frame(x) %>% dplyr::select("seqnames", "start", "end", "name")
    colnames(inp) <- c("chr", "start", "end", "score")
    out <- dplyr::semi_join(inp, keep, by = "chr")
    final <- out$score
  })
  # read in bin sizes
  bs_final <- d[[treated_samp_label]] %>%
    as.data.frame() %>%
    dplyr::select(1:3) %>%
    setNames(c("chr", "start", "end")) %>%
    dplyr::semi_join(keep, by = "chr") %>%
    dplyr::filter(chr != "chrM") %>%
    makeGRangesFromDataFrame()
  # get max values
  mxs_final <- bind_cols(raw) %>% apply(1, max)
  # load exclusion factor based on raw count cut-off and overlap with blacklisted region
  k_final <- mxs_final > raw_count_cutoff & !overlapsAny(bs_final, bl)
  ### generate bigwig template ###
  # loop through each file and remove bins based on exclusion factor (k_final)
  # for the raw bins, we're just keeping the scores
  raw <- lapply(d, function(x) {
    inp <- as.data.frame(x) %>% dplyr::select("seqnames", "start", "end", "name")
    colnames(inp) <- c("chr", "start", "end", "score")
    out <- dplyr::semi_join(inp, keep, by = "chr")
    out <- GenomicRanges::makeGRangesFromDataFrame(out, keep.extra.columns = TRUE)
    out <- out[k_final]
  })
  # final bins without the score. we will merge the scores after they are normalized and scaled
  bw <- lapply(d, function(x) {
    inp <- as.data.frame(x) %>% dplyr::select("seqnames", "start", "end", "name")
    colnames(inp) <- c("chr", "start", "end", "score")
    imd <- dplyr::semi_join(inp, keep, by = "chr")
    out <- imd[, 1:3] %>% mutate(start = start + 1)
    out <- out %>%
      GenomicRanges::makeGRangesFromDataFrame(seqinfo = gn)
    out <- out[k_final]
    k <- !overlapsAny(out, bl)
    out <- out[!overlapsAny(out, bl)]
  })
  # ensure the scores per bin are numeric
  raw[[treated_samp_label]]$score <- as.numeric(raw[[treated_samp_label]]$score)
  raw[[control_label]]$score <- as.numeric(raw[[control_label]]$score)
  # normalize scores per bin by the library depth and scale if necessary
  if (use_control == TRUE) {
    bw[[treated_samp_label]]$score <- (log2(((raw[[treated_samp_label]]$score * to_scale) / treated_samp_library_size + pseudocount) / (raw[[control_label]]$score / control_library_size + pseudocount)))
  } else {
    bw[[treated_samp_label]]$score <- (log2(((raw[[treated_samp_label]]$score * to_scale) / treated_samp_library_size) + pseudocount))
  }
  bw[[treated_samp_label]] <- bw[[treated_samp_label]][!is.na(bw[[treated_samp_label]]$score)]
  # output directory
  out_dir <- paste0(out_dir)
  # output bigwig file
  return(rtracklayer::export.bw(bw[[treated_samp_label]], con = sprintf("%s/%s.%s.%s%snorm.bw", out_dir, cell_line, treated_samp_label, mark, window_size)))
  print("Normalized bigWig file created!")
}
