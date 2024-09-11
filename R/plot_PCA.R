#!/usr/bin/env Rscript
#' Title
#' @title Plot PCA.
#' @description Plots the results from principal component analysis using samples that have been normalized.
#' @param ... file paths to normalized bigWig files. Users can list as many files as needed for the samples to be plotted.
#' @param out_dir a character string specifying the output directory for the scatterplot.
#' @param output_filename a character string specifying the file name for the PCA plot to be saved on disk.
#' @param colors a vector of colors corresponding to each file. The length should match the number of files.
#' @param sample_labels a vector of sample labels corresponding to the inputted files. If not supplied, defaults to the input filenames. The length should match the number of files.
#' @param plot_title a character string specifying the title of the plot.
#' @param font_size an integer specifying the font sizes for the text in the plot.
#' @param point_size an integer specifying the dot size in the plot.
#' @param plot_height a numeric specifying the height of the scatterplot. If not supplied, defaults to 5.
#' @param plot_width a numeric specifying the width of the scatterplot. If not supplied, defaults to 5.
#' @return a PCA plot
#' @export
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' file1 <- "file1.bw"
#' file2 <- "file2.bw"
#' colors <- c("red", "blue")
#' sample_labels <- c("samp1","samp2)
#' result <- plotPCA(file1, file2, sample_labels=sample_labels, colors = colors,out_dir=outdir,output_filename=output,plot_title=title)
#' }
plot_PCA <- function(...,
                     out_dir,
                     output_filename,
                     colors,
                     plot_title,
                     sample_labels = NULL,
                     font_size = NULL,
                     point_size = NULL,
                     plot_width = NULL,
                     plot_height = NULL) {
  # list of normalized bigWig files
  file_list <- list(...)

  # check if at least 2 normalized bigWig files are provided
  if (length(file_list) < 2){
    stop("Less than two bigWig files provided! Please provide at least two files. At least two samples are needed for the PCA plot.")
  }
  # set default values for plot
  font_size <- ifelse(is.null(font_size),12,as.integer(font_size))
  point_size <- ifelse(is.null(point_size),3,as.integer(point_size))
  plot_width <- ifelse(is.null(plot_width),20,as.numeric(plot_width))
  plot_height <- ifelse(is.null(plot_height),15,as.numeric(plot_height))

  if (length(colors) != length(file_list)) {
    stop("The number of colors must match the number of files.")
  }

  # Define the function to perform subsetByOverlaps on a list of files
  subsetByOverlapsForFiles <- function(file_list) {
    # Initialize a list to store the final subsetted GRanges objects with names
    subsetted_granges_list <- list()

    # Initialize a list to store GRanges objects
    granges_list <- list()

    # Loop over each file, read it and store it as a GRanges object
    for (i in seq_along(file_list)) {
      # Read the file (assuming it's in BED format)
      granges <- rtracklayer::import.bw(file_list[[i]])

      # Store the GRanges object in the list
      granges_list[[i]] <- granges
    }

    # Start with the first GRanges object
    common_regions <- granges_list[[1]]

    # Iteratively subset common regions using the other GRanges objects
    for (i in 2:length(granges_list)) {
      # Apply subsetByOverlaps between common regions and the next file
      common_regions <- IRanges::subsetByOverlaps(common_regions, granges_list[[i]])
    }

    # Subset each original GRanges object by the final common regions and keep names
    for (i in seq_along(file_list)) {
      # Subset by the common regions
      subsetted_granges <- IRanges::subsetByOverlaps(granges_list[[i]], common_regions)

      # Get the file name without directory or extension
      file_name <- tools::file_path_sans_ext(basename(file_list[[i]]))

      # Store the subsetted GRanges in the list with the file name
      subsetted_granges_list[[file_name]] <- subsetted_granges
    }

    # Return the named list of subsetted GRanges objects
    return(subsetted_granges_list)
  }

  # ensure only genomic regions found in all samples are included
  subsetted_list <- subsetByOverlapsForFiles(file_list)

  # sample labels
  if (is.null(sample_labels)){
    sample_labels <- as.character(unlist(names(subsetted_list)))
  } else {
    names(subsetted_list) <- sample_labels
  }

  # ensure length of sample labels = length of inputted files
  if (length(sample_labels) != length(subsetted_list)) {
    stop("The length of the sample labels must match the length of the input files.")
  }

  # keep only the signals/scores for each sample
  signal <- lapply(subsetted_list,function(x){
    inp <- as.data.frame(x) %>% dplyr::select("score")
    final <- inp$score
  }) %>% dplyr::bind_cols() %>%
    stats::na.omit()

  # perform principal components analysis on the given data matrix
  res <- stats::prcomp(t(signal), scale. = T, center = T)

  # summarize
  pd <- summary(res)$importance['Proportion of Variance', ] %>%
    Map(function(x, p) sprintf('%s (%.2g%%)', p, x*100), ., names(.)) %>%
    {tmp <- res$x; colnames(tmp) <- .; tmp} %>%
    as.data.frame() %>% dplyr::mutate(samp = rownames(.))
  pd$samp <- factor(x = pd$samp,levels = sample_labels)

  # plot
  p <- pd %>%
  ggplot2::ggplot(ggplot2::aes(x = pd[,1], y = pd[,2],color=samp)) +
  ggplot2::geom_point(size = point_size, show.legend = T) +
  ggplot2::geom_hline(yintercept = 0,alpha=0.1) +
  ggplot2::geom_vline(xintercept = 0,alpha=0.1) +
  ggplot2::scale_color_manual(values = colors) +
  ggplot2::xlab(label = paste(colnames(pd[1]))) +
  ggplot2::ylab(label = paste(colnames(pd[2]))) +
  ggplot2::labs(title = plot_title,colour="Samples") +
  ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill=NA, linewidth=1),
        panel.background = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size=font_size,family="Helvetica",colour = "black"),
        axis.text.y = ggplot2::element_text(size=font_size,family="Helvetica",colour = "black"),
        axis.text.x = ggplot2::element_text(size=font_size,family="Helvetica",colour = "black"),
        axis.line = ggplot2::element_line(linewidth = 0.5, linetype = "solid",colour = "black"),
        plot.title = ggplot2::element_text(hjust = 0.5,color = "black",size=font_size,family="Helvetica"),
        strip.background = ggplot2::element_rect(fill="white"),
        strip.text = ggplot2::element_text(
          size = font_size, color = "black",family = "Helvetica"),
        legend.text = ggplot2::element_text(size=font_size,family = "Helvetica",color = "black"),
        legend.title = ggplot2::element_text(size=font_size,family="Helvetica",color="black"),
        legend.background = ggplot2::element_rect(fill="white"),
        legend.key = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(2.4,2.4,2.4,2.4), "mm"),
        panel.spacing = ggplot2::unit(0.1,'cm'),
        panel.spacing.y = ggplot2::unit(0.1,'cm'),
        panel.spacing.x = ggplot2::unit(0.1,'cm'))

  # save
  ggplot2::ggsave((sprintf("%s/%s_PCA_plot.pdf", out_dir, output_filename)), p, height = plot_height, width = plot_width, device = cairo_pdf, units = "cm")

  # print success message
  print("PCA plot generated!")
}
