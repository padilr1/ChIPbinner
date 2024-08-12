#!/usr/bin/env Rscript
#' Title
#'
#' @title Generate density-based clusters.
#'
#' @description This function identifies similarly-behaving genomic compartments, such as regions with the highest loss or gain of a specific histone mark.
#'
#' @param out_dir Output directory for the scatterplot.
#' @param genome_assembly Must of be one of hg38 or mm10.
#' @param treated_samp_norm_bw The normalized bigwig file for the treated sample.
#' @param wildtype_samp_norm_bw The normalized bigwig file for the wildtype sample.
#' @param cell_line The cell line of the samples.
#' @param histone_mark The broad histone mark used for the analysis.
#' @param number_of_clusters The total number of clusters identified by the HDBSCAN algorithm.
#' @param annotated_clusters The annotated clusters generated using annotate_clust(), which is an R object that needs to be loaded.
#' @param are_R_objects Boolean term (true or false) to indicate if the inputted bigwig files are R objects. It'll use load() for the reps as opposed to reading them in via rtracklayer::import.bed(). Default to FALSE.
#' @param output_filename Filename for the resuting scatterplot.
#' @param title_of_plot The title of the plot.
#' @param pow The power of. Returns the value of x to the power of y (x^y) and this is used for the scales (i.e. show bins if surpassing a certain intensity). Defaults to 1.25.
#' @param show_legend Boolean term whether to show legend or not.
#' @param min The min for the x and y axis.
#' @param max The max for the x and y axis
#' @param bin_size The number of bins represented by each hexagonal shape.
#' @param show_scales Boolean term whether to show scales or not.
#' @param xaxis_label Label for the x-axis. This is normally the baseline label (i.e. WT)
#' @param yaxis_label Label for the y-axis. This is normally the treated sample label.
#' @param height_of_figure A numeric specifying the height of the plots.
#' @param width_of_figure A numeric specifying the width of the plots.
#' @param plot_title_font_size An integer specifying the font size for each plot title.
#' @param legend_font_size An integer specifying the font size for the legends.
#' @param axis_title_font_size An integer specifying the font size for the axis titles.
#' @return a plot of density-based clusters.
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
                                      plot_title_font_size=12,
                                      pow = NULL,
                                      show_legend = FALSE,
                                      legend_font_size = 7,
                                      min,
                                      max,
                                      bin_size,
                                      show_scales = TRUE,
                                      xaxis_label,
                                      yaxis_label,
                                      axis_title_font_size=10,
                                      height_of_figure=6,
                                      width_of_figure=15) {
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
    min <- as.numeric(paste0(min))
    max <- as.numeric(paste0(max))
    bin_size <- as.numeric(paste0(bin_size))
    xaxis_label <- paste0(xaxis_label)
    yaxis_label <- paste0(yaxis_label)
    height_of_figure <- as.numeric(paste0(height_of_figure))
    width_of_figure <- as.numeric(paste0(width_of_figure))
    title_of_plot <- paste0(title_of_plot)
    # font sizes
    plot_title_font_size = as.integer(plot_title_font_size)
    legend_font_size = as.integer(legend_font_size)
    axis_title_font_size = as.integer(axis_title_font_size)
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
        clus[clus == "NA" & overlapsAny(r, cons[[col]])] <- col
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
    olap <- tibble(
      gene = overlapsAny(r, gene),
      igr = overlapsAny(r, igr)
    ) %>%
      mutate(out = case_when(
        gene & !igr ~ 1,
        !gene & igr ~ -1,
        TRUE ~ 0
      )) %>%
      pull(out)

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
      mutate(clr = 1:dplyr::n()) %>%
      reshape2::melt(id.vars = "clr", variable.name = "opa") %>%
      mutate(opa = as.integer(opa))
    leg <- ggplot() +
      geom_tile(aes(x = opa, y = clr, fill = value),
        data = cmat
      ) +
      scale_fill_identity() +
      labs(
        x = "# of bins \u25ba",
        y = "% genic \u25ba"
      ) +
      coord_fixed(expand = FALSE) +
      theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "white", fill = NA, size = 0.5),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x=element_text(family = "Helvetica",
                                  color = "white",
                                  size = legend_font_size,vjust=1.5,hjust=0.25),
        axis.title.y=element_text(family = "Helvetica",
                                  color = "white",
                                  size = legend_font_size,vjust=-1.5,hjust=0.25)
        # axis.title = element_text(
        #   family = "Helvetica",
        #   color = "white",
        #   size = 6,vjust=-3,hjust=0.5
        # )
      )

    ps <- split(d, c(cell_line, cell_line)) %>%
      lapply(function(x) {
        pdat <- lapply(x[2:1], function(y) {
          IRanges::findOverlaps(r, y) %>%
            to() %>%
            {
              y[.]
            } %>%
            score()
        }) %>%
          bind_cols() %>%
          `names<-`(c("x", "y")) %>%
          mutate(r = olap) %>%
          dplyr::filter(
            x > quantile(x, .01),
            x < quantile(x, .99),
            y > quantile(y, .01),
            y < quantile(y, .99)
          )


        hex <- hexbin::hexbin(pdat$x, pdat$y, xbins = 75, IDs = T)
        pdat$cell <- hex@cID
        hex <- data.frame(hexbin::hcell2xy(hex),
          cell = hex@cell,
          count = hex@count
        )

        t1 <- "Intergenic vs genic ratio"
        pdat2 <- pdat %>%
          group_by(cell) %>%
          summarise(prop = mean(r, na.rm = T)) %>%
          ungroup() %>%
          right_join(hex, by = "cell") %>%
          mutate(
            logcount = log10(count),
            ttl = t1
          )
        lim <- data.frame(
          x = c(
            min(pdat2$x[pdat2$count > 10]),
            max(pdat2$x[pdat2$count > 10])
          ) %>%
            scales::expand_range(mul = .05),
          y = c(
            min(pdat2$y[pdat2$count > 10]),
            max(pdat2$y[pdat2$count > 10])
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
        if (show_scales == TRUE) {
          pm <- ggplot() +
            geom_point(aes(x = x, y = y, color = c),
              alpha = 0, data = tdat,
              show.legend = T
            ) +
            geom_hex(aes(x = x, y = y, fill = prop, alpha = count, color = prop),
              stat = "identity", color = NA, data = pdat2, size = 5
            ) +
            scale_fill_gradientn(
              colors = gradientn1, name = "% Genic",
              breaks = leg.brks, labels = leg.labs,
              limits = c(-1, 1)
            ) +
            scale_color_gradientn(
              colors = gradientn1, name = "% Genic",
              breaks = leg.brks, labels = leg.labs,
              limits = c(-1, 1)
            ) +
            scale_alpha(
              range = c(0.01, 1), name = "Number of bins", guide = FALSE,
              trans = scales::trans_new(
                "square",
                function(x) {
                  x^(pow)
                },
                function(x) x^(1 / pow)
              )
            ) +
            coord_cartesian(
              expand = FALSE,
              xlim = lim$x,
              ylim = lim$y
            ) +
            facet_grid(. ~ ttl) +
            labs(
              x = sprintf("%s \u25ba", xaxis_label),
              y = sprintf("%s \u25ba", yaxis_label)
            ) +
            annotate("segment",
              x = -Inf, xend = Inf, y = Inf,
              yend = Inf, color = "white"
            ) +
            theme(
              panel.background = element_rect(fill = "black"),
              legend.position = "none",
              panel.grid = element_blank(),
              plot.background = element_blank(),
              text = element_text(color = "black",family = "Helvetica",size = plot_title_font_size),
              legend.background = element_blank(),
              legend.margin = margin(0.015, 0, 0, 0, unit = "npc"),
              axis.line = element_blank(),
              plot.title = element_blank(),
              strip.text = element_text(color = "white"),
              strip.background = element_rect(fill = "black"),
              axis.title = element_text(family = "Helvetica", color = "black",size=axis_title_font_size)
            )

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
          bnds <- sf::st_bbox(bdat %>% filter(level > 1))

          bdat$ttl <- title_of_plot
          if (show_legend == FALSE) {
            pl <- ggplot() +
              geom_sf(data = bdat, aes(fill = level), color = NA) +
              scale_fill_gradientn(
                colors = pals::magma(10),
                name = "# of bins \u25ba"
              ) +
              coord_sf(
                expand = FALSE,
                xlim = lim$x,
                ylim = lim$y
              ) +
              labs(
                x = sprintf("%s \u25ba", xaxis_label),
                y = sprintf("%s \u25ba", yaxis_label)
              ) +
              facet_grid(. ~ ttl) +
              guides(fill = guide_colorbar(
                ticks.colour = NA,
                frame.colour = "white", frame.linewidth = 0.25,
                barwidth = 2.5, barheight = 0.5,
                title.position = "top", title.hjust = 0.5,
                title.vjust = -1,
                alpha = 0.8
              )) +
              annotate("segment",
                x = -Inf, xend = Inf, y = Inf,
                yend = Inf, color = "white"
              ) +
              theme(
                panel.background = element_rect(fill = "black"),
                panel.grid = element_blank(),
                plot.background = element_blank(),
                text = element_text(color = "black",family = "Helvetica",size = plot_title_font_size),
                legend.title = element_text(
                  color = "white", family = "Helvetica",
                  size = legend_font_size
                ),
                legend.justification = c(0, 1),
                legend.background = element_blank(),
                legend.position = "none",
                legend.direction = "horizontal",
                legend.text = element_blank(),
                plot.title = element_blank(),
                strip.text = element_text(color = "white"),
                strip.background = element_rect(fill = "black"),
                axis.title = element_text(family = "Helvetica",size = axis_title_font_size)
              )
          } else {
            pl <- ggplot() +
              geom_sf(data = bdat, aes(fill = level), color = NA) +
              scale_fill_gradientn(
                colors = pals::magma(10),
                name = "# of bins \u25ba"
              ) +
              coord_sf(
                expand = FALSE,
                xlim = lim$x,
                ylim = lim$y
              ) +
              labs(
                x = sprintf("%s \u25ba", xaxis_label),
                y = sprintf("%s \u25ba", yaxis_label)
              ) +
              facet_grid(. ~ ttl) +
              guides(fill = guide_colorbar(
                ticks.colour = NA,
                frame.colour = "white", frame.linewidth = 0.25,
                barwidth = 2.5, barheight = 0.5,
                title.position = "top", title.hjust = 0.5,
                title.vjust = -1,
                alpha = 0.8
              )) +
              annotate("segment",
                x = -Inf, xend = Inf, y = Inf,
                yend = Inf, color = "white"
              ) +
              theme(
                panel.background = element_rect(fill = "black"),
                panel.grid = element_blank(),
                plot.background = element_blank(),
                text = element_text(color = "black",family = "Helvetica",size = plot_title_font_size),
                legend.title = element_text(
                  color = "white", family = "Helvetica",
                  size = legend_font_size
                ),
                legend.justification = c(0, 1),
                legend.background = element_blank(),
                legend.position = c(0.05, 0.95),
                legend.direction = "horizontal",
                legend.text = element_blank(),
                plot.title = element_blank(),
                strip.text = element_text(color = "white"),
                strip.background = element_rect(fill = "black"),
                axis.title = element_text(family = "Helvetica",size = axis_title_font_size)
              )
          }
          pl$coordinates$aspect <- function(f) {
            NULL
          }


          d <- lapply(x[2:1], function(y) {
            IRanges::findOverlaps(r, y) %>%
              to() %>%
              {
                y[.]
              } %>%
              score()
          }) %>%
            bind_cols() %>%
            `names<-`(c("x", "y")) %>%
            mutate(clu = clus) %>%
            dplyr::filter(
              x > quantile(x, .01),
              x < quantile(x, .99),
              y > quantile(y, .01),
              y < quantile(y, .99)
            )

          hull <- na.omit(d) %>%
            group_by(clu) %>%
            dplyr::slice(chull(x, y))
          lab <- na.omit(d) %>%
            group_by(clu) %>%
            summarise(
              x = mean(x),
              y = mean(y)
            )

          bdat2 <- bdat
          bdat2$ttl <- "Density-based clusters"
          d$ttl <- bdat2$ttl[1]
          gradientn <- paste0("#FFFFFF", sprintf("%02X", c(0, round(seq(0, 255, length.out = nbrks - 1)))))
          pr <- ggplot() +
            geom_sf(data = bdat2, aes(fill = level), color = NA) +
            ggrastr::geom_point_rast(
              data = d, aes(x = x, y = y, color = clu),
              alpha = .1, size = .02, raster.dpi = 300,
              raster.height = 5, raster.width = 5
            ) +
            ggrepel::geom_label_repel(
              aes(
                label = clu, x = x, y = y,
                color = clu
              ),
              data = lab,
              show.legend = FALSE,
              point.padding = NA, box.padding = 0,
              direction = "x", inherit.aes = FALSE
            ) +
            scale_fill_gradientn(colors = gradientn) +
            scale_color_manual(values = cclr) +
            coord_sf(
              expand = FALSE,
              xlim = lim$x,
              ylim = lim$y
            ) +
            facet_grid(. ~ ttl) +
            annotate("segment",
              x = -Inf, xend = Inf, y = Inf,
              yend = Inf, color = "white"
            ) +
            labs(
              x = sprintf("%s \u25ba", xaxis_label),
              y = sprintf("%s \u25ba", yaxis_label)
            ) +
            theme(
              panel.background = element_rect(fill = "black"),
              legend.position = "none",
              panel.grid = element_blank(),
              plot.background = element_blank(),
              text = element_text(color = "black",family = "Helvetica",size = plot_title_font_size),
              legend.background = element_blank(),
              plot.title = element_blank(),
              strip.text = element_text(color = "white"),
              strip.background = element_rect(fill = "black"),
              legend.margin = margin(0.015, 0, 0, 0, unit = "npc"),
              axis.title = element_text(family = "Helvetica",size=axis_title_font_size)
            )
          if (show_legend == FALSE) {
            pr$coordinates$aspect <- function(f) {
              NULL
            }
            # lo <- c(
            #   patchwork::area(t = 1, l = 1, b = 20, r = 20),
            #   patchwork::area(t = 2, l = 2, b = 8, r = 8)
            # )
            # pm_leg <- pm + leg + plot_layout(design = lo)
            patchwork::wrap_plots(pl, pm, pr, nrow = 1)
          } else {
            pr$coordinates$aspect <- function(f) {
              NULL
            } # patchwork::area(t = 1, l = 1, b = 20, r = 20),
            # patchwork::area(t = 2, l = 2, b = 8, r = 8)
            lo <- c(
              patchwork::area(t = 1, l = 1, b = 20, r = 20),
              patchwork::area(t = 2, l = 2, b = 7, r = 7)
            )
            pm_leg <- pm + leg + plot_layout(design = lo)
            patchwork::wrap_plots(pl, pm_leg, pr, nrow = 1)
          }
        } else {
          pm <- ggplot() +
            geom_point(aes(x = x, y = y, color = c),
              alpha = 0, data = tdat,
              show.legend = T
            ) +
            geom_hex(aes(x = x, y = y, fill = prop, alpha = count, color = prop),
              stat = "identity", color = NA, data = pdat2, size = 5
            ) +
            scale_fill_gradientn(
              colors = gradientn1, name = "% Genic",
              breaks = leg.brks, labels = leg.labs,
              limits = c(-1, 1)
            ) +
            scale_color_gradientn(
              colors = gradientn1, name = "% Genic",
              breaks = leg.brks, labels = leg.labs,
              limits = c(-1, 1)
            ) +
            scale_alpha(
              range = c(0.01, 1), name = "Number of bins", guide = FALSE,
              trans = scales::trans_new(
                "square",
                function(x) {
                  x^(pow)
                },
                function(x) x^(1 / pow)
              )
            ) +
            coord_cartesian(
              expand = FALSE,
              xlim = lim$x,
              ylim = lim$y
            ) +
            facet_grid(. ~ ttl) +
            labs(
              x = sprintf("%s \u25ba", xaxis_label),
              y = sprintf("%s \u25ba", yaxis_label)
            ) +
            annotate("segment",
              x = -Inf, xend = Inf, y = Inf,
              yend = Inf, color = "white"
            ) +
            theme(
              panel.background = element_rect(fill = "black"),
              legend.position = "none",
              panel.grid = element_blank(),
              plot.background = element_blank(),
              text = element_text(color = "black",family = "Helvetica",size = plot_title_font_size),
              legend.background = element_blank(),
              legend.margin = margin(0.015, 0, 0, 0, unit = "npc"),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_blank(),
              strip.text = element_text(color = "white"),
              strip.background = element_rect(fill = "black"),
              axis.title = element_text(family = "Helvetica", color = "black",size = axis_title_font_size)
            )

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
          bnds <- sf::st_bbox(bdat %>% filter(level > 1))

          bdat$ttl <- title_of_plot
          if (show_legend == FALSE) {
            pl <- ggplot() +
              geom_sf(data = bdat, aes(fill = level), color = NA) +
              scale_fill_gradientn(
                colors = pals::magma(10),
                name = "# of bins \u25ba"
              ) +
              coord_sf(
                expand = FALSE,
                xlim = lim$x,
                ylim = lim$y
              ) +
              labs(
                x = sprintf("%s \u25ba", xaxis_label),
                y = sprintf("%s \u25ba", yaxis_label)
              ) +
              facet_grid(. ~ ttl) +
              guides(fill = guide_colorbar(
                ticks.colour = NA,
                frame.colour = "white", frame.linewidth = 0.25,
                barwidth = 2.5, barheight = 0.5,
                title.position = "top", title.hjust = 0.5,
                title.vjust = -1,
                alpha = 0.8
              )) +
              annotate("segment",
                x = -Inf, xend = Inf, y = Inf,
                yend = Inf, color = "white"
              ) +
              theme(
                panel.background = element_rect(fill = "black"),
                panel.grid = element_blank(),
                plot.background = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                text = element_text(color = "black",family = "Helvetica",size = plot_title_font_size),
                legend.title = element_text(
                  color = "white", family = "Helvetica",
                  size = legend_font_size
                ),
                legend.justification = c(0.1, 1),
                legend.background = element_blank(),
                legend.position = "none",
                legend.direction = "horizontal",
                legend.text = element_blank(),
                plot.title = element_blank(),
                strip.text = element_text(color = "white"),
                strip.background = element_rect(fill = "black"),
                axis.title = element_text(family = "Helvetica",size = axis_title_font_size)
              )
          } else {
            pl <- ggplot() +
              geom_sf(data = bdat, aes(fill = level), color = NA) +
              ggplot2::scale_fill_gradientn(
                colors = pals::magma(10),
                name = "# of bins \u25ba"
              ) +
              coord_sf(
                expand = FALSE,
                xlim = lim$x,
                ylim = lim$y
              ) +
              labs(
                x = sprintf("%s \u25ba", xaxis_label),
                y = sprintf("%s \u25ba", yaxis_label)
              ) +
              facet_grid(. ~ ttl) +
              guides(fill = guide_colorbar(
                ticks.colour = NA,
                frame.colour = "white", frame.linewidth = 0.25,
                barwidth = 2.5,
                barheight = 0.5,
                title.position = "top", title.hjust = 0.5,
                title.vjust = -1,
                alpha = 0.8
              )) +
              annotate("segment",
                x = -Inf, xend = Inf, y = Inf,
                yend = Inf, color = "white"
              ) +
              theme(
                panel.background = element_rect(fill = "black"),
                panel.grid = element_blank(),
                plot.background = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                text = element_text(color = "black",family = "Helvetica",size = plot_title_font_size),
                legend.title = element_text(
                  color = "white", family = "Helvetica",
                  size = legend_font_size
                ),
                legend.justification = c(0.1, 1),
                legend.background = element_blank(),
                legend.position = c(0.05, 0.95),
                legend.direction = "horizontal",
                legend.text = element_blank(),
                plot.title = element_blank(),
                strip.text = element_text(color = "white"),
                strip.background = element_rect(fill = "black"),
                axis.title = element_text(family = "Helvetica",size=axis_title_font_size)
              )
          }
          pl$coordinates$aspect <- function(f) {
            NULL
          }


          d <- lapply(x[2:1], function(y) {
            IRanges::findOverlaps(r, y) %>%
              to() %>%
              {
                y[.]
              } %>%
              score()
          }) %>%
            bind_cols() %>%
            `names<-`(c("x", "y")) %>%
            mutate(clu = clus) %>%
            dplyr::filter(
              x > quantile(x, .01),
              x < quantile(x, .99),
              y > quantile(y, .01),
              y < quantile(y, .99)
            )

          hull <- na.omit(d) %>%
            group_by(clu) %>%
            dplyr::slice(chull(x, y))
          lab <- na.omit(d) %>%
            group_by(clu) %>%
            summarise(
              x = mean(x),
              y = mean(y)
            )

          bdat2 <- bdat
          bdat2$ttl <- "Density-based clusters"
          d$ttl <- bdat2$ttl[1]
          gradientn <- paste0("#FFFFFF", sprintf("%02X", c(0, round(seq(0, 255, length.out = nbrks - 1)))))
          pr <- ggplot() +
            geom_sf(data = bdat2, aes(fill = level), color = NA) +
            ggrastr::geom_point_rast(
              data = d, aes(x = x, y = y, color = clu),
              alpha = .1, size = .02, raster.dpi = 300,
              raster.height = 5, raster.width = 5
            ) +
            ggrepel::geom_label_repel(
              aes(
                label = clu, x = x, y = y,
                color = clu
              ),
              data = lab,
              show.legend = FALSE,
              point.padding = NA, box.padding = 0,
              direction = "x", inherit.aes = FALSE
            ) +
            scale_fill_gradientn(colors = gradientn) +
            scale_color_manual(values = cclr) +
            coord_sf(
              expand = FALSE,
              xlim = lim$x,
              ylim = lim$y
            ) +
            facet_grid(. ~ ttl) +
            annotate("segment",
              x = -Inf, xend = Inf, y = Inf,
              yend = Inf, color = "white"
            ) +
            labs(
              x = sprintf("%s \u25ba", xaxis_label),
              y = sprintf("%s \u25ba", yaxis_label)
            ) +
            theme(
              panel.background = element_rect(fill = "black"),
              legend.position = "none",
              panel.grid = element_blank(),
              plot.background = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              text = element_text(color = "black",family = "Helvetica",size = plot_title_font_size),
              legend.background = element_blank(),
              plot.title = element_blank(),
              strip.text = element_text(color = "white"),
              strip.background = element_rect(fill = "black"),
              legend.margin = margin(0.015, 0, 0, 0, unit = "npc"),
              axis.title = element_text(family = "Helvetica",size=axis_title_font_size)
            )
          if (show_legend == FALSE) {
            pr$coordinates$aspect <- function(f) {
              NULL
            }
            # lo <- c(
            #   patchwork::area(t = 1, l = 1, b = 20, r = 20),
            #   patchwork::area(t = 2, l = 2, b = 8, r = 8)
            # )
            # pm_leg <- pm + leg + plot_layout(design = lo)
            patchwork::wrap_plots(pl, pm, pr, nrow = 1)
          } else {
            pr$coordinates$aspect <- function(f) {
              NULL
            }
            lo <- c(
              # patchwork::area(t = 1, l = 1, b = 20, r = 20)
              # patchwork::area(t = 2, l = 2, b = 8, r = 8)
              # https://patchwork.data-imaginist.com/reference/area.html
              # This means that t and l should always be less or equal to b and r respectively
              patchwork::area(t = 1, l = 1, b = 20, r = 20),
              patchwork::area(t = 2, l = 2, b = 7, r = 7)
            )
            pm_leg <- pm + leg + plot_layout(design = lo)
            patchwork::wrap_plots(pl, pm_leg, pr, nrow = 1)
          }
        }
      })
    # ggplot2::ggsave((sprintf("%s/%s_density_scatterplot.pdf", out_dir, output_filename)), ps[[cell_line]], height = 2.4, width = 6.2, device = cairo_pdf)
    ggplot2::ggsave((sprintf("%s/%s_density_scatterplot.pdf", out_dir, output_filename)), ps[[cell_line]], height = height_of_figure, width = width_of_figure, device = cairo_pdf, units = "cm")
    # ggplot2::ggsave((sprintf("%s/%s_density_scatterplot.png", out_dir, output_filename)), ps[[cell_line]], height = height_of_figure, width = width_of_figure, device = "png", dpi = 600, units = "cm")
    print("Density-based scatterplots generated!")
  })
}
