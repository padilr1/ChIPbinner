#!/usr/bin/env Rscript
#' Title
#' @title Generates correlation plot.
#' @description Plots a correlation matrix-based hierarchical clustering.
#' @param ... file paths to normalized bigWig files. Users can list as many files as needed for the samples to be plotted.
#' @param out_dir a character string specifying the output directory for the scatterplot.
#' @param output_filename a character string specifying the file name for the PCA plot to be saved on disk.
#' @param sample_labels a vector of sample labels corresponding to the inputted files. If not supplied, defaults to the input filenames. The length should match the number of files.
#' @param correlation_method a character string specifying which correlation coefficient is to be computed. Must be one of "pearson" (default) or "spearman".
#' @param plot_title a character string specifying the title of the plot.
#' @param key_title a character string specifying the legend title.
#' @param font_size an integer specifying the font sizes for the text in the plot.
#' @param row_label_angle an integer specifying the angle of the labels for the rows. If not supplied, defaults to 270.
#' @param col_label_angle an integer specifying the angle of the labels for the columns. If not supplied, defaults to 360.
#' @param row_adj_label a 2-element vector giving the (left-right, top-bottom) justification of row labels (relative to the text orientation). Defaults to c(0.5,0.5).
#' @param col_adj_label a 2-element vector giving the (left-right, top-bottom) justification of column labels (relative to the text orientation). Defaults to c(0.5,0.5).
#' @param colors a vector of colors corresponding to each file. The length should match the number of files.
#' @param plot_height a numeric specifying the height of the scatterplot. If not supplied, defaults to 10.
#' @param plot_width a numeric specifying the width of the scatterplot. If not supplied, defaults to 8.
#' @return a correlation plot
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
#' result <- plot_correlation(file1, file2, sample_labels = sample_labels,colors = colors,out_dir=outdir,output_filename=output,plot_title=title)
#' }
plot_correlation <- function(...,
                     out_dir,
                     output_filename,
                     sample_labels = NULL,
                     correlation_method="pearson",
                     colors,
                     plot_title = NULL,
                     key_title = NULL,
                     font_size = NULL,
                     col_label_angle = NULL,
                     row_label_angle = NULL,
                     row_adj_label=NULL,
                     col_adj_label=NULL,
                     plot_width = NULL,
                     plot_height = NULL) {
  # list of normalized bigWig files
  file_list <- list(...)

  # check if at least 2 normalized bigWig files are provided
  if (length(file_list) < 2){
    stop("Less than two bigWig files provided! Please provide at least two files. At least two samples are needed for the correlation plot.")
  }
  # set default values for plot
  plot_title <- ifelse(is.null(plot_title),NA,plot_title)
  font_size <- ifelse(is.null(font_size),1,as.integer(font_size))
  col_label_angle <- ifelse(is.null(col_label_angle),360,as.integer(col_label_angle))
  row_label_angle <- ifelse(is.null(row_label_angle),270,as.integer(row_label_angle))
  row_adj_label <- ifelse(is.null(row_adj_label),c(0.5,0.5),row_adj_label)
  col_adj_label <- ifelse(is.null(col_adj_label),c(0.5,0.5),col_adj_label)
  key_title <- ifelse(is.null(key_title),NA,as.character(key_title))
  plot_width <- ifelse(is.null(plot_width),10,as.numeric(plot_width))
  plot_height <- ifelse(is.null(plot_height),8,as.numeric(plot_height))

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
    sample_labels <- sample_labels
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

  # save
  pdf(file = (sprintf("%s/%s_correlation_plot.pdf", out_dir, output_filename)),width = plot_width,height = plot_height,bg = "white")

  # plot correlation matrix-based hierarchical clustering
 signal %>% stats::cor(method=correlation_method) %>% gplots::heatmap.2(density.info = "none",key.title=key_title,main=plot_title,key.xlab=NA,key.ylab = NA,ColSideColors = colors,RowSideColors = colors,trace="none",scale="none",keysize = 1,cexRow = font_size,cexCol = font_size,srtCol = col_label_angle,srtRow = row_label_angle,offsetRow=0,offsetCol = 0,labRow = sample_labels,labCol = sample_labels,adjRow = row_adj_label,adjCol = col_adj_label)

 # device off
 dev.off()

  # print success message
  print("Correlation plot generated!")
}
