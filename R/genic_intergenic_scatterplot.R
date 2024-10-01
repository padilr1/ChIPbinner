#!/usr/bin/env Rscript
#+ message = FALSE, warning = FALSE
#' Title
#' @title Plot a scatterplot of bins annotated as genic and intergenic.
#' @description Using two normalized bigWig files generated using 'norm.bw', the user can generate a scatterplot to compare bins across two conditions (treated sample versus baseline sample (such as wildtype)). The bins will be annotated according to their overlap with genic and intergenic regions.
#' @param out_dir a character string specifying the output directory for the scatterplot.
#' @param genome_assembly a character string specifying the genome assembly. Allowed values include "hg38" or "mm10".
#' @param cell_line a character string specifying the cell line of the samples.
#' @param histone_mark a character string specifying the broad histone mark being analyzed.
#' @param treated_samp_norm_bw a character string specifying the normalized bigwig file for the treated sample.
#' @param wildtype_samp_norm_bw a character string specifying the normalized bigwig file for the wildtype sample.
#' @param are_R_objects a logical indicating whether the inputted bigwig files are R objects. It'll use load() for the reps as opposed to reading them in via rtracklayer::import.bed(). Defaults to FALSE.
#' @param output_filename a character string specifying the file name for the resulting scatterplot to be saved on disk.
#' @param title_of_plot a character string specifying title of the plot.
#' @param max_x a numeric specifying the maximum value for the x-axis.
#' @param max_y a numeric specifying the maximum value for the y-axis.
#' @param min_x a numeric specifying the minimum value for the x-axis.
#' @param min_y a numeric specifying the minimum value for the y-axis.
#' @param pow a numeric specifying the power of. Returns the value of x to the power of y (x^y) and this is used for the scales (i.e. show bins if surpassing a certain intensity). Defaults to 1.25.
#' @param xaxis_label a character string specifying the label for the x-axis. This is normally the wildtype sample.
#' @param yaxis_label a character string specifying the label for the y-axis. This is normally the treated sample.
#' @param show_scales a logical indicating whether to show scales or not. Defaults to TRUE.
#' @param show_legend a logical indicating whether to show legend or not. Defaults to false.
#' @param legend_pos a character string specifying the legend position. Allowed values include "left" or "right". Defaults to "left".
#' @param height_of_plot a numeric specifying the height of the scatterplot. If not supplied, defaults to 5.
#' @param width_of_plot a numeric specifying the width of the scatterplot. If not supplied, defaults to 5.
#'
#' @return Returns a scatterplot of bins annotated as genic or intergenic comparing two samples.
#' @export
#'
#' @include norm_bw.R
#' @example inst/examples/example_genic_intergenic_scatterplot.R
genic_intergenic_scatterplot <- function(out_dir,
                                         genome_assembly,
                                         cell_line,
                                         histone_mark,
                                         treated_samp_norm_bw,
                                         wildtype_samp_norm_bw,
                                         are_R_objects=FALSE,
                                         output_filename,
                                         title_of_plot,
                                         max_x,
                                         max_y,
                                         pow = NULL,
                                         xaxis_label,
                                         yaxis_label,
                                         min_x = NULL,
                                         min_y = NULL,
                                         show_scales = TRUE,
                                         show_legend = FALSE,
                                         legend_pos = NULL,
                                         height_of_plot = NULL,
                                         width_of_plot = NULL) {
  suppressWarnings({
  # gene and intergenic reference files
  if (genome_assembly == "hg38") {
    gene <- hg38_gene
    igr <- hg38_igr
  } else if (genome_assembly == "mm10") {
    gene <- mm10_gene
    igr <- mm10_igr
  }
  # samples labels
  treated_samp_label <- basename(tools::file_path_sans_ext(treated_samp_norm_bw))
  wildtype_samp_label <- basename(tools::file_path_sans_ext(wildtype_samp_norm_bw))
  # title of plot
  title_of_plot <- paste0(title_of_plot)
  # x and y axis labels
  xaxis_label <- paste0(xaxis_label)
  yaxis_label <- paste0(yaxis_label)
  # parameters for max x and y
  max_x <- as.numeric(max_x)
  max_y <- as.numeric(max_y)
  # parameters for plot dimension
  ## height
  if (is.null(height_of_plot)) {
    height_of_plot <- as.numeric(5)
  } else {
    height_of_plot <- as.numeric(paste0(height_of_plot))
  }
  ## width
  if (is.null(width_of_plot)) {
    width_of_plot <- as.numeric(5)
  } else {
    width_of_plot <- as.numeric(paste0(width_of_plot))
  }
  # loading function
  loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  # check if inputted files are R objects or not
  if (are_R_objects == "FALSE"){
    # load via import.bed
    treated_samp <- rtracklayer::import.bw(paste0(treated_samp_norm_bw))
    wildtype_samp <- rtracklayer::import.bw(paste0(wildtype_samp_norm_bw))
  } else if (are_R_objects == "TRUE") {
    # otherwise load via loading function
    treated_samp <- loadRData(paste0(treated_samp_norm_bw))
    wildtype_samp <- loadRData(paste0(wildtype_samp_norm_bw))
  }
  d <- list(treated_samp,wildtype_samp)
  names(d) <- c(treated_samp_label,wildtype_samp_label)
  # only keep regions found in both
  d[[wildtype_samp_label]] <- IRanges::subsetByOverlaps(d[[wildtype_samp_label]], d[[treated_samp_label]])
  d[[treated_samp_label]] <- IRanges::subsetByOverlaps(d[[treated_samp_label]], d[[wildtype_samp_label]])
  lapply(d, length)
  r <- lapply(d, function(y) y[y$score != 0]) %>%
    Reduce(function(a, b) a[IRanges::overlapsAny(a, b)], .) %>%
    GenomicRanges::granges()
  # overlaps
  olap <- tibble::tibble(
    gene = IRanges::overlapsAny(r, gene),
    igr = IRanges::overlapsAny(r, igr)
  ) %>%
    dplyr::mutate(out = dplyr::case_when(
      gene & !igr ~ 1,
      !gene & igr ~ -1,
      TRUE ~ 0
    )) %>%
    dplyr::pull(out)
  # parameters for the plot
  lineclr <- "black"
  horz <- F
  gradientn1 <- pals::brewer.rdylbu(50)
  cramp <- colorRampPalette(c("#000000ff", "#ffffff00"), alpha = T)(5)
  leg.brks <- seq(-1, 1, length.out = 19)[seq(2, 18, by = 2)]
  leg.labs <- c(
    sprintf("Genic\u25bc"), rep("", 3), "50%",
    rep("", 3), sprintf("Genic\u25b2")
  )
  len <- 9
  pal <- pals::brewer.rdylbu(len)
  cmat <- seq(0, 255, length.out = len + 1) %>%
    {
      .[-1]
    } %>%
    round() %>%
    as.hexmode() %>%
    format(width = 2, upper.case = T) %>%
    lapply(function(x) {
      paste0(pal, x)
    }) %>%
    do.call(cbind, .) %>%
    data.frame() %>%
    `colnames<-`(1:len) %>%
    dplyr::mutate(clr = 1:dplyr::n()) %>%
    reshape2::melt(id.vars = "clr", variable.name = "opa") %>%
    dplyr::mutate(opa = as.integer(opa))
  # legend
  leg <- ggplot2::ggplot() +
    ggplot2::geom_tile(ggplot2::aes(x = opa, y = clr, fill = value),
      data = cmat
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(
      x = "# of bins \u25ba",
      y = "% genic \u25ba"
    ) +
    ggplot2::coord_fixed(expand = F) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "white", fill = NA, size = 0.1),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(
        family = "Helvetica",
        color = "black",
        size = 10
      )
    )
  leg
  ggplot2::ggsave(filename = sprintf("%s/legend.png", out_dir), plot = leg, height = 1, width = 1, device = "png", dpi = 600, bg = "white")
  # pow
  if (is.null(pow)) {
    pow <- as.numeric(1.25)
  } else {
    pow <- as.numeric(paste0(pow))
  }
  # make plot
  ps <- split(d, c(cell_line,cell_line)) %>%
    lapply(function(x) {
      pdat <- lapply(x[2:1], function(y) {
        IRanges::findOverlaps(r, y) %>%
          IRanges::to() %>%
          {
            y[.]
          } %>%
          GenomicRanges::score()
      }) %>%
        dplyr::bind_cols() %>%
        `names<-`(c("x", "y")) %>%
        dplyr::mutate(r = olap) %>%
        dplyr::filter(
          x > quantile(x, .01),
          x < quantile(x, .99),
          y > quantile(y, .01),
          y < quantile(y, .99)
        )
      # colour gradient for bins
      hex <- hexbin::hexbin(pdat$x, pdat$y, xbins = 100, IDs = T)
      pdat$cell <- hex@cID
      hex <- data.frame(hexbin::hcell2xy(hex),
        cell = hex@cell,
        count = hex@count
      )

      t1 <- "Intergenic vs genic ratio"
      pdat2 <- pdat %>%
        dplyr::group_by(cell) %>%
        dplyr::summarise(prop = mean(r, na.rm = T)) %>%
        dplyr::ungroup() %>%
        dplyr::right_join(hex, by = "cell") %>%
        dplyr::mutate(
          logcount = log10(count),
          ttl = t1
        )
      if (is.null(min_x)) {
        min_x <- as.numeric(min(pdat2$x[pdat2$count > 10]))
      } else {
        min_x <- as.numeric(paste0(min_x))
      }
      if (is.null(min_y)) {
        min_y <- as.numeric(min(pdat2$y[pdat2$count > 10]))
      } else {
        min_y <- as.numeric(paste0(min_y))
      }
      lim <- data.frame(
        x = c(
          min_x,
          max_x
        ) %>%
          scales::expand_range(mul = .05),
        y = c(
          min_y,
          max_y
        ) %>%
          scales::expand_range(mul = .05)
      )
      yrang <- diff(lim$y)
      tdat <- data.frame(
        x = rep(median(pdat$x), 2),
        y = rep(median(pdat$y), 2),
        c = c(0, 1),
        ttl = t1
      )
      # pow = power of; returns the value of x to the power of y (x^y); used for the scales (i.e. show bins if surpassing a certain intensity)
      pow <- pow
      # if user wants to show the scale for x and y axes
      if (show_scales == TRUE) {
        pm <- ggplot2::ggplot() +
          ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = c),
            alpha = 0, data = tdat,
            show.legend = T
          ) +
          ggplot2::geom_hex(ggplot2::aes(x = x, y = y, fill = prop, alpha = count, color = prop),
            stat = "identity", color = NA, data = pdat2, size = 5
          ) +
          ggplot2::scale_fill_gradientn(
            colors = gradientn1, name = "% Genic",
            breaks = leg.brks, labels = leg.labs,
            limits = c(-1, 1)
          ) +
          ggplot2::scale_color_gradientn(
            colors = gradientn1, name = "% Genic",
            breaks = leg.brks, labels = leg.labs,
            limits = c(-1, 1)
          ) +
          ggplot2::scale_alpha(
            range = c(0.01, 5), name = "Number of bins", guide = F,
            trans = scales::trans_new(
              "square",
              function(x) {
                x^(pow)
              },
              function(x) x^(1 / pow)
            )
          ) +
          ggplot2::coord_cartesian(
            expand = F,
            xlim = lim$x,
            ylim = lim$y
          ) +
          ggplot2::labs(
            x = sprintf("%s \u25ba", xaxis_label),
            y = sprintf("%s \u25ba", yaxis_label), title = title_of_plot
          ) +
          ggplot2::annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = Inf, color = "gray") +
          ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = "white"),
            panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1),
            legend.position = "none",
            panel.grid = ggplot2::element_blank(),
            plot.background = ggplot2::element_blank(),
            legend.background = ggplot2::element_blank(),
            legend.margin = ggplot2::margin(0.015, 0, 0, 0, unit = "npc"),
            plot.title = ggplot2::element_text(hjust = 0.5, color = "black", size = 12, family = "Helvetica"),
            strip.text = ggplot2::element_text(color = "black"),
            strip.background = ggplot2::element_rect(fill = "black"),
            axis.title = ggplot2::element_text(family = "Helvetica", color = "black")
          )
      } else {
        pm <- ggplot2::ggplot() +
          ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = c),
            alpha = 0, data = tdat,
            show.legend = T
          ) +
          ggplot2::geom_hex(ggplot2::aes(x = x, y = y, fill = prop, alpha = count, color = prop),
            stat = "identity", color = NA, data = pdat2, size = 5
          ) +
          ggplot2::scale_fill_gradientn(
            colors = gradientn1, name = "% Genic",
            breaks = leg.brks, labels = leg.labs,
            limits = c(-1, 1)
          ) +
          ggplot2::scale_color_gradientn(
            colors = gradientn1, name = "% Genic",
            breaks = leg.brks, labels = leg.labs,
            limits = c(-1, 1)
          ) +
          ggplot2::scale_alpha(
            range = c(0.01, 5), name = "Number of bins", guide = F,
            trans = scales::trans_new(
              "square",
              function(x) {
                x^(pow)
              },
              function(x) x^(1 / pow)
            )
          ) +
          ggplot2::coord_cartesian(
            expand = F,
            xlim = lim$x,
            ylim = lim$y
          ) +
          ggplot2::labs(
            x = sprintf("%s \u25ba", xaxis_label),
            y = sprintf("%s \u25ba", yaxis_label), title = title_of_plot
          ) +
          ggplot2::annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = Inf, color = "gray") +
          ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = "white"),
            panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1),
            legend.position = "none",
            panel.grid = ggplot2::element_blank(),
            plot.background = ggplot2::element_blank(),
            legend.background = ggplot2::element_blank(),
            legend.margin = ggplot2::margin(0.015, 0, 0, 0, unit = "npc"),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.line = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5, color = "black", size = 12, family = "Helvetica"),
            strip.text = ggplot2::element_text(color = "black"),
            strip.background = ggplot2::element_rect(fill = "black"),
            axis.title = ggplot2::element_text(family = "Helvetica", color = "black")
          )
      }
      pdat <- MASS::kde2d(x = pdat$x, y = pdat$y, n = 100)
      brks <- pretty(c(pdat$z), 30)
      nbrks <- length(brks)
      fac <- diff(brks[1:2]) * nrow(d) / sum(pdat$z)
      b <- isoband::isobands(x = pdat$x, y = pdat$y, z = t(pdat$z), brks[1:(nbrks - 1)], brks[2:nbrks])
      bands <- isoband::iso_to_sfg(b)
      bdat <- sf::st_sf(level = 1:length(bands), geometry = sf::st_sfc(bands))
      if (!all(sf::st_is_valid(bdat))) {
        bdat <- sf::st_make_valid(bdat)
      }
      bnds <- sf::st_bbox(bdat %>% dplyr::filter(level > 1))

      bdat$ttl <- title_of_plot
      if (show_legend == TRUE) {
        # if no legend position is specified
        if (is.null(legend_pos)) {
          pm + patchwork::inset_element(leg,left = 0,bottom = 0.6,right = 0.4,top = 1)
        }
        # else it is right
        else if (legend_pos == "right") {
          pm + patchwork::inset_element(leg, left = 0.6, bottom = 0.6, right = 1, top = 1)
        }
        # else it is left
        else if (legend_pos == "left") {
          pm + patchwork::inset_element(leg,left = 0,bottom = 0.6,right = 0.4,top = 1)
        } else {
          print("The only two options for 'legend_pos' is 'right' or 'left'")
        }
      } else{
        patchwork::wrap_plots(pm, nrow = 1)
      }
    })
  ps[[cell_line]]
  # save plot
  # ggsave((sprintf("%s/%s.scatterplot.png",out_dir, output_filename)), ps[[cell_line]], height = height_of_plot, width = width_of_plot, device = "png", dpi = 600)
  ggplot2::ggsave((sprintf("%s/%s.scatterplot.pdf",out_dir, output_filename)), ps[[cell_line]], height = height_of_plot, width = width_of_plot, device = cairo_pdf)
  # print completed statement
  print("Generated genic/intergenic scatterplot!")
  })
}
