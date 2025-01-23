#!/usr/bin/env Rscript
#' Title
#'
#' @title Plot scatterplots of annotated genic/intergenic regions and density-based clusters.
#'
#' @description Plots scatterplots of annotated genic/intergenic regions and previously identified clusters of similar-behaving genomic regions.
#'
#' @param out_dir a character string specifying the output directory for the scatterplot.
#' @param genome_assembly a character string specifying the genome assembly. Allowed values include "hg38" or "mm10".
#' @param treated_samp_norm_bw a character string specifying the normalized bigwig file for the treated sample.
#' @param wildtype_samp_norm_bw a character string specifying the normalized bigwig file for the wildtype sample.
#' @param cell_line a character string specifying the cell line of the samples.
#' @param histone_mark a character string specifying the broad histone mark being analyzed.
#' @param number_of_clusters an integer specifying the total number of clusters identified by the HDBSCAN algorithm.
#' @param annotated_clusters a character string specifying the R object containing the annotated clusters generated using annotate_clust().
#' @param are_R_objects a logical indicating whether the inputted bigwig files are R objects. It'll use load() for the reps as opposed to reading them in via rtracklayer::import.bed(). Defaults to FALSE.
#' @param output_filename a character string specifying the file name for the resuting scatterplots to be saved on disk.
#' @param title_of_plot a character string specifying title of the plot.
#' @param pow a numeric specifying the power of. Returns the value of x to the power of y (x^y) and this is used for the scales (i.e. show bins if surpassing a certain intensity). Defaults to 1.25.
#' @param max_x a numeric specifying the maximum value for the x-axis.
#' @param max_y a numeric specifying the maximum value for the y-axis.
#' @param min_x a numeric specifying the minimum value for the x-axis.
#' @param min_y a numeric specifying the minimum value for the y-axis.
#' @param hexbins an integer specifying the number of hexagons across the x axis. The number of hexagons is determined by this integer: a larger value increases the number of hexagons, each containing fewer genomic bins, while a smaller value reduces the number of hexagons, with each containing more genomic bins. The default is 150.
#' @param show_scales a logical indicating whether to show scales or not. Defaults to TRUE.
#' @param xaxis_label a character string specifying the label for the x-axis. This is normally the wildtype sample.
#' @param yaxis_label a character string specifying the label for the y-axis. This is normally the treated sample.
#' @param height_of_figure A numeric specifying the height of the plots. If not supplied, defaults to 5.
#' @param width_of_figure A numeric specifying the width of the plots. If not supplied, defaults to 5.
#' @param plot_title_font_size An integer specifying the font size for each plot title. If not supplied, defaults to 12.
#' @param legend_font_size An integer specifying the font size for the legends. If not supplied, defaults to 7.
#' @param axis_title_font_size An integer specifying the font size for the axis titles. If not supplied, defaults to 10.
#' @param include_additional_density_plot A logical indicating whether to include an additional density plot that displays only the number of bins in non-annotated regions. Defaults to FALSE.
#' @param filter_extreme_bins A logical specifying whether to filter out the top and bottom 1% of bins when plotting. Defaults to TRUE.
#' @return plots of annotated genic/intergenic bins and previously identified density-based clusters.
#' @export
#'
#' @include annotate_clust.R
#' @include norm_bw.R
#'
#' @example inst/examples/example_density_based_scatterplot.R
density_based_scatterplot <- function(out_dir,
                                      genome_assembly,
                                      treated_samp_norm_bw,
                                      wildtype_samp_norm_bw,
                                      are_R_objects = FALSE,
                                      cell_line,
                                      histone_mark,
                                      annotated_clusters,
                                      number_of_clusters,
                                      output_filename,
                                      title_of_plot,
                                      plot_title_font_size = 12,
                                      pow = NULL,
                                      legend_font_size = 7,
                                      min_x,
                                      min_y,
                                      max_x,
                                      max_y,
                                      hexbins = NULL,
                                      show_scales = TRUE,
                                      xaxis_label,
                                      yaxis_label,
                                      axis_title_font_size = 10,
                                      height_of_figure = 6,
                                      width_of_figure = 15,
                                      include_additional_density_plot=FALSE,
                                      filter_extreme_bins = TRUE) {
  suppressWarnings({
    # directory parameters
    out_dir <- paste0(out_dir)
    # gene and intergenic reference files
    if (genome_assembly == "hg38") {
      gene <- hg38_gene
      igr <- hg38_igr
    } else if (genome_assembly == "mm10") {
      gene <- mm10_gene
      igr <- mm10_igr
    }
    # samples info
    cell_line <- paste0(cell_line)
    mark <- paste0(histone_mark)
    # samples labels
    treated_samp_label <- basename(tools::file_path_sans_ext(treated_samp_norm_bw))
    wildtype_samp_label <- basename(tools::file_path_sans_ext(wildtype_samp_norm_bw))
    # number of cluster
    number_of_clusters <- as.integer(paste0(number_of_clusters))
    # parameters for plot
    max_x <- as.numeric(max_x)
    min_x <- as.numeric(min_x)
    max_y <- as.numeric(max_y)
    min_y <- as.numeric(min_y)
    hexbins <- as.numeric(paste0(hexbins))
    xaxis_label <- paste0(xaxis_label)
    yaxis_label <- paste0(yaxis_label)
    height_of_figure <- as.numeric(paste0(height_of_figure))
    width_of_figure <- as.numeric(paste0(width_of_figure))
    title_of_plot <- paste0(title_of_plot)
    # font sizes
    plot_title_font_size <- as.integer(plot_title_font_size)
    legend_font_size <- as.integer(legend_font_size)
    axis_title_font_size <- as.integer(axis_title_font_size)
    # pow
    if (is.null(pow)) {
      pow <- as.numeric(1.25)
    } else {
      pow <- as.numeric(paste0(pow))
    }
    # loading function
    loadRData <- function(fileName) {
      # loads an RData file, and returns it
      load(fileName)
      get(ls()[ls() != "fileName"])
    }
    # check if inputted files are R objects or not
    if (are_R_objects == "FALSE") {
      # load via import.bed
      treated_samp <- rtracklayer::import.bw(paste0(treated_samp_norm_bw))
      wildtype_samp <- rtracklayer::import.bw(paste0(wildtype_samp_norm_bw))
    } else if (are_R_objects == "TRUE") {
      # otherwise load via loading function
      treated_samp <- loadRData(paste0(treated_samp_norm_bw))
      wildtype_samp <- loadRData(paste0(wildtype_samp_norm_bw))
    }
    d <- list(treated_samp, wildtype_samp)
    names(d) <- c(treated_samp_label, wildtype_samp_label)
    # subset by overlaps
    d[[wildtype_samp_label]] <- IRanges::subsetByOverlaps(d[[wildtype_samp_label]], d[[treated_samp_label]])
    d[[treated_samp_label]] <- IRanges::subsetByOverlaps(d[[treated_samp_label]], d[[wildtype_samp_label]])
    lapply(d, length)
    r <- lapply(d, function(y) y[y$score != 0]) %>%
      Reduce(function(a, b) a[IRanges::overlapsAny(a, b)], .) %>%
      GenomicRanges::granges()
    # load annotated clusters from using the annotate_clust function
    cons <- loadRData(paste0(annotated_clusters))
    # Initialize an empty list of cluster labels
    clust_labs <- c()

    # Loop to generate lists and append them to master_list
    for (i in 1:number_of_clusters) {
      # Generate a new list (example)
      new_list <- c(toupper(letters[i]))

      # Append new_list to master_list
      clust_labs <- append(clust_labs, new_list)
    }

    # Initialize clus as NA (or any default value)
    clus <- rep("NA", length(r))

    # Loop through each column and update clus if overlap is found
    for (col in clust_labs) {
      if (col %in% names(cons)) { # Check if column exists in cons
        clus[clus == "NA" & IRanges::overlapsAny(r, cons[[col]])] <- col
      }
    }

    clus[clus == "NA"] <- NA

    # colours for the clusters
    cclr <- setNames(
      pals::alphabet((number_of_clusters))[seq(1:number_of_clusters)],
      # c("C", "B", "A")
      clust_labs
    )
    # genic and intergenic overlap with enriched regions
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

    lineclr <- "black"
    horz <- FALSE
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

    # code for legend
    leg <- ggplot2::ggplot() +
      ggplot2::geom_tile(ggplot2::aes(x = opa, y = clr, fill = value),
        data = cmat
      ) +
      ggplot2::scale_fill_identity() +
      ggplot2::labs(
        x = "# of bins \u25ba",
        y = "% genic \u25ba"
      ) +
      ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = "white", fill = NA, size = 0.5),
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(
          family = "Helvetica",
          color = "white",
          size = legend_font_size, vjust = 1.5, hjust = 0.25
        ),
        axis.title.y = ggplot2::element_text(
          family = "Helvetica",
          color = "white",
          size = legend_font_size, vjust = -1.5, hjust = 0.25
        )
        # axis.title = element_text(
        #   family = "Helvetica",
        #   color = "white",
        #   size = 6,vjust=-3,hjust=0.5
        # )
      )

    ########### start of code for plots #############

    ps <- split(d, c(cell_line, cell_line)) %>%
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
          dplyr::mutate(r = olap)

        if (filter_extreme_bins == TRUE){
          pdat <- pdat %>%
            dplyr::filter(
              x > quantile(x, .01),
              x < quantile(x, .99),
              y > quantile(y, .01),
              y < quantile(y, .99)
            )
        } else if (filter_extreme_bins == FALSE) {
          pdat = pdat
        }

        # colour gradient for bins
        if (is.null(hexbins)) {
          hexbins <- 100
        } else {
          hexbins <- as.integer(hexbins)
        }
        hex <- hexbin::hexbin(pdat$x, pdat$y, xbins = hexbins, IDs = T)
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

        # min and max values for x- and y-axes

        # if (is.null(min_x)) {
        #   min_x <- as.numeric(min(pdat2$x[pdat2$count > 10]))
        # } else {
        #   min_x <- as.numeric(paste0(min_x))
        # }
        # if (is.null(min_y)) {
        #   min_y <- as.numeric(min(pdat2$y[pdat2$count > 10]))
        # } else {
        #   min_y <- as.numeric(paste0(min_y))
        # }
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

        pow <- pow

        # remove and add scales for pm
        pm_with_scales <- function(x){
          x +
            ggplot2::theme(
              panel.background = ggplot2::element_rect(fill = "black"),
              legend.position = "none",
              panel.grid = ggplot2::element_blank(),
              plot.background = ggplot2::element_blank(),
              text = ggplot2::element_text(color = "black", family = "Helvetica", size = plot_title_font_size),
              legend.background = ggplot2::element_blank(),
              legend.margin = ggplot2::margin(0.015, 0, 0, 0, unit = "npc"),
              axis.line = ggplot2::element_blank(),
              plot.title = ggplot2::element_blank(),
              strip.text = ggplot2::element_text(color = "white"),
              strip.background = ggplot2::element_rect(fill = "black"),
              axis.title = ggplot2::element_text(family = "Helvetica", color = "black", size = axis_title_font_size)
            )
        }

        pm_without_scales <- function(x){
          x +
            ggplot2::theme(
              panel.background = ggplot2::element_rect(fill = "black"),
              legend.position = "none",
              panel.grid = ggplot2::element_blank(),
              plot.background = ggplot2::element_blank(),
              text = ggplot2::element_text(color = "black", family = "Helvetica", size = plot_title_font_size),
              legend.background = ggplot2::element_blank(),
              legend.margin = ggplot2::margin(0.015, 0, 0, 0, unit = "npc"),
              axis.line = ggplot2::element_blank(),
              axis.text = ggplot2::element_blank(),
              axis.ticks = ggplot2::element_blank(),
              plot.title = ggplot2::element_blank(),
              strip.text = ggplot2::element_text(color = "white"),
              strip.background = ggplot2::element_rect(fill = "black"),
              axis.title = ggplot2::element_text(family = "Helvetica", color = "black", size = axis_title_font_size)
            )
        }

        # genic/intergenic plot
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
            range = c(0.01, 1), name = "Number of bins", guide = FALSE,
            trans = scales::trans_new(
              "square",
              function(x) {
                x^(pow)
              },
              function(x) x^(1 / pow)
            )
          ) +
          ggplot2::coord_cartesian(
            expand = FALSE,
            xlim = lim$x,
            ylim = lim$y
          ) +
          ggplot2::facet_grid(. ~ ttl) +
          ggplot2::labs(
            x = sprintf("%s \u25ba", xaxis_label),
            y = sprintf("%s \u25ba", yaxis_label)
          ) +
          ggplot2::annotate("segment",
            x = -Inf, xend = Inf, y = Inf,
            yend = Inf, color = "white"
          )

        if (show_scales == TRUE) {
          pm <- pm_with_scales(pm)
        } else if (show_scales == FALSE) {
          pm <- pm_without_scales(pm)
        }

        # for magma plot

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
        bnds <- sf::st_bbox(bdat %>% dplyr::filter(level > 0))

        bdat$ttl <- title_of_plot

        # magma plot with and without scales
        pl_with_scales <- function(x){
          x +
            ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = "black"),
            panel.grid = ggplot2::element_blank(),
            plot.background = ggplot2::element_blank(),
            text = ggplot2::element_text(color = "black", family = "Helvetica", size = plot_title_font_size),
            legend.title = ggplot2::element_text(
              color = "white", family = "Helvetica",
              size = legend_font_size
            ),
            legend.justification = c(0, 1),
            legend.background = ggplot2::element_blank(),
            legend.position = c(0.05, 0.95),
            legend.direction = "horizontal",
            legend.text = ggplot2::element_blank(),
            plot.title = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(color = "white"),
            strip.background = ggplot2::element_rect(fill = "black"),
            axis.title = ggplot2::element_text(family = "Helvetica", size = axis_title_font_size)
          )}
        pl_without_scales <- function(x){
            x +
              ggplot2::theme(
                panel.background = ggplot2::element_rect(fill = "black"),
                panel.grid = ggplot2::element_blank(),
                plot.background = ggplot2::element_blank(),
                text = ggplot2::element_text(color = "black", family = "Helvetica", size = plot_title_font_size),
                legend.title = ggplot2::element_text(
                  color = "white", family = "Helvetica",
                  size = legend_font_size
                ),
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                legend.justification = c(0, 1),
                legend.background = ggplot2::element_blank(),
                legend.position = c(0.05, 0.95),
                legend.direction = "horizontal",
                legend.text = ggplot2::element_blank(),
                plot.title = ggplot2::element_blank(),
                strip.text = ggplot2::element_text(color = "white"),
                strip.background = ggplot2::element_rect(fill = "black"),
                axis.title = ggplot2::element_text(family = "Helvetica", size = axis_title_font_size)
              )}

        pl <- ggplot2::ggplot() +
          ggplot2::geom_sf(data = bdat, ggplot2::aes(fill = level), color = NA) +
          ggplot2::scale_fill_gradientn(
            colors = pals::magma(10),
            name = "# of bins \u25ba"
          ) +
          ggplot2::coord_sf(
            expand = FALSE,
            xlim = lim$x,
            ylim = lim$y
          ) +
          ggplot2::labs(
            x = sprintf("%s \u25ba", xaxis_label),
            y = sprintf("%s \u25ba", yaxis_label)
          ) +
          ggplot2::facet_grid(. ~ ttl) +
          ggplot2::guides(fill = ggplot2::guide_colorbar(
            ticks.colour = NA,
            frame.colour = "white", frame.linewidth = 0.25,
            barwidth = 2.5, barheight = 0.5,
            title.position = "top", title.hjust = 0.5,
            title.vjust = -1,
            alpha = 0.8
          )) +
          ggplot2::annotate("segment",
            x = -Inf, xend = Inf, y = Inf,
            yend = Inf, color = "white"
          )

        if (show_scales == TRUE) {
          pl <- pl_with_scales(pl)
        } else if (show_scales == FALSE) {
          pl <- pl_without_scales(pl)
        }

        pl$coordinates$aspect <- function(f) {
          NULL
        }

        d <- lapply(x[2:1], function(y) {
          IRanges::findOverlaps(r, y) %>%
            IRanges::to() %>%
            {
              y[.]
            } %>%
            GenomicRanges::score()
        }) %>%
          dplyr::bind_cols() %>%
          `names<-`(c("x", "y")) %>%
          dplyr::mutate(clu = clus)

        if (filter_extreme_bins == TRUE){
          d <- d %>%
            dplyr::filter(
              x > quantile(x, .01),
              x < quantile(x, .99),
              y > quantile(y, .01),
              y < quantile(y, .99)
            )
        } else if (filter_extreme_bins == FALSE){
          d <- d
        }

        hull <- na.omit(d) %>%
          dplyr::group_by(clu) %>%
          dplyr::slice(chull(x, y))
        lab <- na.omit(d) %>%
          dplyr::group_by(clu) %>%
          dplyr::summarise(
            x = mean(x),
            y = mean(y)
          )

        bdat2 <- bdat
        bdat2$ttl <- "Density-based clusters"
        d$ttl <- bdat2$ttl[1]
        gradientn <- paste0("#FFFFFF", sprintf("%02X", c(0, round(seq(0, 255, length.out = nbrks - 1)))))

        # density-based cluster plot scaling and no scaling
        with_scaling_text_for_density_cluster <- function(x) {
          x +
            ggplot2::theme(
              panel.background = ggplot2::element_rect(fill = "black"),
              legend.position = "none",
              panel.grid = ggplot2::element_blank(),
              plot.background = ggplot2::element_blank(),
              text = ggplot2::element_text(color = "black", family = "Helvetica", size = plot_title_font_size),
              legend.background = ggplot2::element_blank(),
              plot.title = ggplot2::element_blank(),
              strip.text = ggplot2::element_text(color = "white"),
              strip.background = ggplot2::element_rect(fill = "black"),
              legend.margin = ggplot2::margin(0.015, 0, 0, 0, unit = "npc"),
              axis.title = ggplot2::element_text(family = "Helvetica", size = axis_title_font_size)
            )
        }
        without_scaling_text_for_density_cluster <- function(x) {
          x +
            ggplot2::theme(
              panel.background = ggplot2::element_rect(fill = "black"),
              legend.position = "none",
              panel.grid = ggplot2::element_blank(),
              # to remove scales
              axis.text = ggplot2::element_blank(),
              axis.ticks = ggplot2::element_blank(),
              plot.background = ggplot2::element_blank(),
              text = ggplot2::element_text(color = "black", family = "Helvetica", size = plot_title_font_size),
              legend.background = ggplot2::element_blank(),
              plot.title = ggplot2::element_blank(),
              strip.text = ggplot2::element_text(color = "white"),
              strip.background = ggplot2::element_rect(fill = "black"),
              legend.margin = ggplot2::margin(0.015, 0, 0, 0, unit = "npc"),
              axis.title = ggplot2::element_text(family = "Helvetica", size = axis_title_font_size)
            )
        }

        #### pr = density-based plot ####
        pr <- ggplot2::ggplot() +
          ggplot2::geom_sf(data = bdat2, ggplot2::aes(fill = level), color = NA) +
          ggrastr::geom_point_rast(
            data = d, ggplot2::aes(x = x, y = y, color = clu),
            # alpha = .1, size = .05, raster.dpi = 600,
            alpha = 0.01, size = .05, raster.dpi = 600,
            raster.height = 5, raster.width = 5
          ) +
          ggrepel::geom_label_repel(
            ggplot2::aes(
              label = clu, x = x, y = y,
              color = clu
            ),
            data = lab,
            show.legend = FALSE,
            point.padding = NA, box.padding = 0,
            direction = "x", inherit.aes = FALSE
          ) +
          ggplot2::scale_fill_gradientn(colors = gradientn) +
          ggplot2::scale_color_manual(values = cclr) +
          ggplot2::coord_sf(
            expand = FALSE,
            xlim = lim$x,
            ylim = lim$y
          ) +
          ggplot2::facet_grid(. ~ ttl) +
          ggplot2::annotate("segment",
            x = -Inf, xend = Inf, y = Inf,
            yend = Inf, color = "white"
          ) +
          ggplot2::labs(
            x = sprintf("%s \u25ba", xaxis_label),
            y = sprintf("%s \u25ba", yaxis_label)
          )
        if (show_scales == TRUE) {
          pr <- with_scaling_text_for_density_cluster(pr)
        } else if (show_scales == FALSE) {
          pr <- without_scaling_text_for_density_cluster(pr)
        }
        pr$coordinates$aspect <- function(f) {
          NULL
        } # patchwork::area(t = 1, l = 1, b = 20, r = 20),
        # patchwork::area(t = 2, l = 2, b = 8, r = 8)
        lo <- c(
          patchwork::area(t = 1, l = 1, b = 20, r = 20),
          patchwork::area(t = 2, l = 2, b = 7, r = 7)
        )
        pm_leg <- pm + leg + patchwork::plot_layout(design = lo)
        if(include_additional_density_plot == TRUE){
        patchwork::wrap_plots(pl, pm_leg, pr, nrow = 1)
        } else if(include_additional_density_plot == FALSE){
          patchwork::wrap_plots(pm_leg, pr, nrow = 1)
        }
      })
    ggplot2::ggsave((sprintf("%s/%s_density_scatterplot.pdf", out_dir, output_filename)), ps[[cell_line]], height = height_of_figure, width = width_of_figure, device = cairo_pdf, units = "cm")
    print("Density-based scatterplots generated!")
  })
}
