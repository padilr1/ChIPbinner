#!/usr/bin/env Rscript
#+ message = FALSE, warning = FALSE
#' Title
#' @title Generate normalized bigWig files.
#'
#' @description Generate bigWig files with normalized bin scores needed for subsequent functions in the package workflow.
#'
#' @param out_dir Output directory for normalized bigWig files.
#' @param genome_assembly Must of be one of hg38 or mm10.
#' @param chip_samp_binned_file Binned BED file associated with the sample with ChIP-seq or Cut&Run signal.
#' @param chip_samp_library_size Mapped library size for the sample with ChIP-seq or Cut&Run signal.
#' @param control_binned_file Binned BED file associated with the control sample.
#' @param control_library_size Mapped library size for the control sample.
#' @param use_control Boolean term, whether to use a control sample or not.
#' @param pseudocount Pseudocount to avoid division by 0. If not defined, defaults to 1e-3.
#' @param raw_count_cutoff Raw read count cutoff to exclude bins depleted in signal across all tracks. For example, if 10 is inputted, this removes raw read count consistently lower than 10. Defaults to 0.
#' @param scaling_factor Optional: quantatitive normalization/scaling of the raw binned signal by this factor. For example, this could be ChIP-Rx values or genome-wide modification percentage values obtained from mass spectrometry values.
#'
#' @return a single bigWig file with normalized binned scores.
#' @export
#'
#' @example inst/examples/example_norm_bw.R
norm_bw <- function(out_dir,
                    genome_assembly,
                    chip_samp_binned_file,
                    chip_samp_library_size,
                    control_binned_file = NULL,
                    control_library_size = NULL,
                    use_control,
                    pseudocount = NULL,
                    raw_count_cutoff = NULL,
                    scaling_factor = NULL) {
  suppressWarnings({
    # reference files
    ## blacklist (bl) and chrom_sizes are internal R objects that are pre-loaded with the package
    if (genome_assembly == "hg38") {
      bl <- hg38_bl
      chrom_sizes <- hg38_chrom_sizes
    } else if (genome_assembly == "mm10") {
      bl <- mm10_bl
      chrom_sizes <- mm10_chrom_sizes
    }
    # remove file extensions to generate chip samp label
    chip_samp_label <- basename(tools::file_path_sans_ext(chip_samp_binned_file))
    # import binned BED file
    chip_samp <- rtracklayer::import.bed(paste0(chip_samp_binned_file))
    # pseudocount to avoid division by 0
    if (is.null(pseudocount)) {
      pseudocount <- as.numeric(1e-3)
    } else {
      pseudocount <- as.numeric(paste0(pseudocount))
    }
    # raw count cut-off to remove bins with low counts across samples
    if (is.null(raw_count_cutoff)) {
      raw_count_cutoff <- as.numeric(0)
    } else {
      raw_count_cutoff <- as.numeric(paste0(raw_count_cutoff))
    }
    # scaling factor
    if (is.null(scaling_factor)) {
      to_scale <- 1
    } else {
      to_scale <- as.numeric(paste0(scaling_factor))
    }
    # chip library size
    chip_samp_library_size <- as.numeric(paste0(chip_samp_library_size))
    # if using control/input
    if (is.null(control_binned_file)) {
      print("Not using input")
      print(paste0("ChIP/Cut&Run sample sequencing depth=", chip_samp_library_size))
      print(paste0("ChIP/Cut&Run sample=", chip_samp_label))
    } else {
      # print chip samp info
      print(paste0("ChIP/Cut&Run sample=", chip_samp_label))
      print(paste0("ChIP/Cut&Run sample sequencing depth=", chip_samp_library_size))
      # control lib size
      control_library_size <- as.numeric(paste0(control_library_size))
      # control label
      control_label <- basename(tools::file_path_sans_ext(control_binned_file))
      # print control info
      print(paste0("Control sample=", control_label))
      print(paste0("Control sample sequencing depth=", control_library_size))
      # import control binned bed file
      control <- import.bed(control_binned_file)
      d <- list(chip_samp, control)
      names(d) <- c(chip_samp_label, control_label)
    }
    ### deriving cut-off for raw counts ###
    # read in chrom sizes
    keep <- chrom_sizes %>%
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
    bs_final <- d[[chip_samp_label]] %>%
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
    raw[[chip_samp_label]]$score <- as.numeric(raw[[chip_samp_label]]$score)
    raw[[control_label]]$score <- as.numeric(raw[[control_label]]$score)
    # normalize scores per bin by the library depth and scale if necessary
    if (use_control == TRUE) {
      bw[[chip_samp_label]]$score <- (log2(((raw[[chip_samp_label]]$score * to_scale) / chip_samp_library_size + pseudocount) / (raw[[control_label]]$score / control_library_size + pseudocount)))
    } else {
      bw[[chip_samp_label]]$score <- (log2(((raw[[chip_samp_label]]$score * to_scale) / chip_samp_library_size) + pseudocount))
    }
    bw[[chip_samp_label]] <- bw[[chip_samp_label]][!is.na(bw[[chip_samp_label]]$score)]
    # rename final object
    chip_samp_norm_bw <- bw[[chip_samp_label]]
    # output directory
    out_dir <- paste0(out_dir)
    # output bigwig file and R object
    return(c(rtracklayer::export.bw(object = chip_samp_norm_bw, con = sprintf("%s/%s.norm.bw", out_dir, chip_samp_label)),save(chip_samp_norm_bw, file = sprintf("%s/%s.norm_bw.rda", out_dir, chip_samp_label))))
    print("Normalized bigWig file created!")
  })
}
