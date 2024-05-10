#!/usr/bin/env Rscript
#+ message = FALSE, warning = FALSE
#+ #!/usr/bin/env Rscript
#+ message = FALSE, warning = FALSE
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(lattice)
library(gridExtra)
library(DiffBind)
library(patchwork)
library(pacman)
library(sf)
library(MASS)
library(lwgeom)
library(ggrepel)
library(hexbin)
library(ggrastr)
library(viridis)
library(pals)
library(isoband)
# start of function
#' Title
#'
#' @param wd
#' @param path_for_norm_bw
#' @param treatment
#' @param control
#' @param cell_line
#' @param histone_mark
#' @param gene
#' @param intergenic
#' @param title_of_plot
#' @param window_size
#' @param number_of_clusters
#' @param pow
#' @param show_legend
#' @param min
#' @param max
#' @param bin_size
#' @param show_scales
#' @param xaxis_label
#' @param yaxis_label
#' @param height_of_figure
#' @param width_of_figure
#'
#' @return
#' @export
#'
#' @examples
density_based_scatterplots <- function(wd,path_for_norm_bw,treatment,control,cell_line,histone_mark,gene,intergenic,title_of_plot,window_size,number_of_clusters,pow=NULL,show_legend=FALSE,min,max,bin_size,show_scales=TRUE,xaxis_label,yaxis_label,height_of_figure,width_of_figure){
# parameters
setwd(wd)
getwd()
control=paste0(control)
test=paste0(treatment)
cell_line=paste0(cell_line)
mark=paste0(histone_mark)
gene <- import.bed(gene)
igr <- import.bed(intergenic)
min=as.numeric(paste0(min))
max=as.numeric(paste0(max))
bin_size=as.numeric(paste0(bin_size))
xaxis_label=paste0(xaxis_label)
yaxis_label=paste0(yaxis_label)
height_of_figure=as.numeric(paste0(height_of_figure))
width_of_figure=as.numeric(paste0(width_of_figure))
# cluster_size <- paste0(cluster_size)
# number_of_clusters <- as.integer(paste0(number_of_clusters))
path <- paste0(path_for_norm_bw)
window_size=paste0(".",window_size,"kb.")
title_of_plot <- paste0(title_of_plot)
#
s <- list.files(path = path, pattern = window_size,full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c("line","samp","mark",NA,NA,NA),'\\.',F) %>%
  mutate(f = file.path(path, f)) %>%
  dplyr::filter(line == cell_line) %>%
  dplyr::filter(samp == control | samp == test)
# cell_line <- as.character(s$line[[1]])
# NEED TO ORDER SAMPLE FIRST THEN PARENTAL #
odr <- c(test,control)
s <- s %>%
  dplyr::slice(match(odr,samp))
#
d <- deframe(s[,c('samp', 'f')]) %>%
  lapply(import.bw)
#
#
d[[control]] <- subsetByOverlaps(d[[control]],d[[test]])
d[[test]]<- subsetByOverlaps(d[[test]],d[[control]])
#
lapply(d,length)
#
#
r <- lapply(d, function(y) y[y$score != 0]) %>%
  Reduce(function(a, b) a[overlapsAny(a, b)], .) %>%
  granges()
#
load((sprintf('data/cons/cons.%s.%s.%s.%s%srda',cell_line,control,test,mark,window_size)))
#
number_of_clusters=as.integer(paste0(number_of_clusters))
# pow
if(is.null(pow)){
  pow = as.numeric(1.25)} else {
    pow = as.numeric(paste0(pow))
  }

if ( number_of_clusters == 2 ){
clus <- case_when(
  overlapsAny(r, cons$A) ~ 'A',
  overlapsAny(r, cons$B) ~ 'B',
  T ~ 'NA'
)
clus[clus == 'NA'] <- NA

olap <- tibble(gene = overlapsAny(r, gene),
               igr = overlapsAny(r, igr)) %>%
  mutate(out = case_when(
    gene & !igr ~ 1,
    !gene & igr ~ -1,
    TRUE ~ 0
  )) %>%
  pull(out)

lineclr <- "black"
horz <- F
gradientn1 <- brewer.rdylbu(50)
cramp <- colorRampPalette(c("#000000ff","#ffffff00"), alpha = T)(5)

leg.brks <- seq(-1, 1, length.out = 19)[seq(2, 18, by = 2)]
leg.labs <- c(sprintf('Genic\u25bc'), rep('', 3), '50%',
              rep('', 3), sprintf('Genic\u25b2'))
len <- 9
pal <- brewer.rdylbu(len)
cmat <- seq(0, 255, length.out = len + 1) %>%
  {.[-1]} %>%
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
            data = cmat) +
  scale_fill_identity() +
  labs(x = "# of bins \u25ba",
       y = "% genic \u25ba") +
  coord_fixed(expand = F) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "white", fill = NA, size = 0.5),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Arial',
                                  color = "white",
                                  size = 7))
#####################
ps <- split(d, s$line) %>%
  lapply(function(x) {
    pdat <- lapply(x[2:1], function(y) {
      findOverlaps(r, y) %>%
        to() %>%
        {y[.]} %>%
        score()
    }) %>%
      bind_cols() %>%
      `names<-`(c('x', 'y')) %>%
      mutate(r = olap) %>%
      dplyr::filter(x > quantile(x, .01),
                    x < quantile(x, .99),
                    y > quantile(y, .01),
                    y < quantile(y, .99))

    # size of bins
    # hex <- hexbin(pdat$x, pdat$y, xbins = 75, IDs = T)
    hex <- hexbin(pdat$x, pdat$y, xbins = bin_size, IDs = T)
    pdat$cell <- hex@cID
    hex <- data.frame(hcell2xy(hex),
                      cell = hex@cell,
                      count = hex@count)

    t1 <- 'Intergenic vs genic ratio'
    pdat2 <- pdat %>%
      group_by(cell) %>%
      summarise(prop = mean(r, na.rm = T)) %>%
      ungroup %>%
      right_join(hex, by = "cell") %>%
      mutate(logcount = log10(count),
             ttl = t1)
    # lim <- data.frame(x = c(min(pdat2$x[pdat2$count > 10]),
    #                         max(pdat2$x[pdat2$count > 10])) %>%
    #                     scales::expand_range(mul = .05),
    #                   y = c(min(pdat2$y[pdat2$count > 10]),
    #                         max(pdat2$y[pdat2$count > 10])) %>%
    #                     scales::expand_range(mul = .05))
    lim <- data.frame(x = c(min,max) %>%
                        scales::expand_range(mul = .05),
                      y = c(min,max) %>%
                        scales::expand_range(mul = .05))
    # lim <- data.frame(x = c(min(pdat2$x[pdat2$count > 10]),
    #                         max) %>%
    #                     scales::expand_range(mul = .05),
    #                   y = c(min(pdat2$y[pdat2$count > 10]),
    #                         max) %>%
    #                     scales::expand_range(mul = .05))
    yrang <- diff(lim$y)
    tdat <- data.frame(x = rep(median(pdat$x), 2),
                       y = rep(median(pdat$y), 2),
                       c = c(0,1),
                       ttl = t1)
    pow <- pow
    # 1.25 ORIGINALLY
    if (show_scales == TRUE){
    pm <- ggplot() +
      geom_point(aes(x = x, y = y, color = c), alpha = 0, data = tdat,
                 show.legend = T) +
      geom_hex(aes(x = x, y = y, fill = prop, alpha = count, color = prop),
               stat = "identity", color = NA, data = pdat2, size = 5) +
      scale_fill_gradientn(colors = gradientn1, name = "% Genic",
                           breaks = leg.brks, labels = leg.labs,
                           limits = c(-1, 1)) +
      scale_color_gradientn(colors = gradientn1, name = "% Genic",
                            breaks = leg.brks, labels = leg.labs,
                            limits = c(-1, 1)) +
      scale_alpha(range = c(0.01, 1), name = "Number of bins", guide = F,
                  trans = scales::trans_new("square",
                                            function(x) {
                                              x^(pow)
                                            },
                                            function(x) x^(1/pow))) +
      coord_cartesian(expand = F,
                      xlim = lim$x,
                      ylim = lim$y) +
      facet_grid(.~ttl) +
      labs(x = sprintf('%s \u25ba',xaxis_label),
           y = sprintf('%s \u25ba',yaxis_label)) +
      annotate('segment', x = -Inf, xend =Inf, y = Inf,
               yend = Inf, color = 'white') +
      annotate("segment", x=-Inf, xend=Inf,y=-Inf, yend=Inf,color="white") +
      theme(panel.background = element_rect(fill = "black"),
            legend.position = "none",
            panel.grid = element_blank(),
            plot.background = element_blank(),
            legend.background = element_blank(),
            legend.margin = margin(0.015, 0, 0, 0, unit="npc"),
            # axis.text = element_blank(),
            # axis.ticks = element_blank(),
            axis.line = element_blank(),
            plot.title = element_blank(),
            strip.text = element_text(color = "white"),
            # strip.background = element_rect(fill = "#726866"),
            strip.background = element_rect(fill = 'black'),
            axis.title = element_text(family = "Arial", color = "black"))

    pdat <- kde2d(x = pdat$x, y = pdat$y, n = 100)
    brks <- pretty(c(pdat$z), 30)
    nbrks <- length(brks)
    fac <- diff(brks[1:2]) * nrow(d)/sum(pdat$z)
    b <- isobands(x = pdat$x, y = pdat$y, z = t(pdat$z), brks[1:(nbrks - 1)], brks[2:nbrks])
    bands <- iso_to_sfg(b)
    bdat <- st_sf(level = 1:length(bands), geometry = st_sfc(bands))
    if(!all(st_is_valid(bdat))) {
      bdat <- lwgeom::st_make_valid(bdat)
    }
    bnds <- st_bbox(bdat %>% filter(level > 1))

    bdat$ttl <- title_of_plot
    if(show_legend == FALSE){
    pl <- ggplot() +
      geom_sf(data = bdat, aes(fill = level), color = NA) +
      scale_fill_gradientn(colors = pals::magma(10),
                           name = "# of bins \u25ba") +
      coord_sf(expand = F,
               xlim = lim$x,
               ylim = lim$y) +
      labs(x = sprintf('%s \u25ba', xaxis_label),
           y = sprintf('%s \u25ba', yaxis_label)) +
      facet_grid(. ~ ttl) +
      guides(fill = guide_colorbar(barwidth = 3, barheight = 0.5,
                                   title.position = "top", title.hjust = 0.5,
                                   frame.colour = "white", frame.linewidth = 1,
                                   ticks = F)) +
      annotate('segment', x = -Inf, xend =Inf, y = Inf,
               yend = Inf, color = 'white') +
      annotate("segment", x=-Inf, xend=Inf,y=-Inf, yend=Inf,color="white")+
      theme(panel.background = element_rect(fill = "black"),
            panel.grid = element_blank(),
            plot.background = element_blank(),
            # axis.text = element_blank(),
            # axis.ticks = element_blank(),
            text = element_text(color = "black"),
            legend.title = element_text(color = "white", family = "Arial",
                                        size = 7),
            legend.justification = c(0, 1),
            legend.background = element_blank(),
            # legend.position = c(0.05, 0.95),
            legend.position="none",
            legend.direction = "horizontal",
            legend.text = element_blank(),
            #axis.ticks = element_line(color = "black"),
            #axis.text = element_text(color = "black"),
            plot.title = element_blank(),
            strip.text = element_text(color = "white"),
            #strip.background = element_rect(fill = "#726866"),
            strip.background = element_rect(fill = 'black'),
            axis.title = element_text(family = "Arial"))
    } else{
      pl <- ggplot() +
        geom_sf(data = bdat, aes(fill = level), color = NA) +
        scale_fill_gradientn(colors = pals::magma(10),
                             name = "# of bins \u25ba") +
        coord_sf(expand = F,
                 xlim = lim$x,
                 ylim = lim$y) +
        labs(x = sprintf('%s \u25ba', xaxis_label),
             y = sprintf('%s \u25ba', yaxis_label)) +
        facet_grid(. ~ ttl) +
        guides(fill = guide_colorbar(barwidth = 3, barheight = 0.5,
                                     title.position = "top", title.hjust = 0.5,
                                     frame.colour = "white", frame.linewidth = 1,
                                     ticks = F)) +
        annotate('segment', x = -Inf, xend =Inf, y = Inf,
                 yend = Inf, color = 'white') +
        annotate("segment", x=-Inf, xend=Inf,y=-Inf, yend=Inf,color="white")+
        theme(panel.background = element_rect(fill = "black"),
              panel.grid = element_blank(),
              plot.background = element_blank(),
              # axis.text = element_blank(),
              # axis.ticks = element_blank(),
              text = element_text(color = "black"),
              legend.title = element_text(color = "white", family = "Arial",
                                          size = 7),
              legend.justification = c(0, 1),
              legend.background = element_blank(),
              legend.position = c(0.05, 0.95),
              # legend.position="none",
              legend.direction = "horizontal",
              legend.text = element_blank(),
              #axis.ticks = element_line(color = "black"),
              #axis.text = element_text(color = "black"),
              plot.title = element_blank(),
              strip.text = element_text(color = "white"),
              #strip.background = element_rect(fill = "#726866"),
              strip.background = element_rect(fill = 'black'),
              axis.title = element_text(family = "Arial"))
    }} else{
      pm <- ggplot() +
        geom_point(aes(x = x, y = y, color = c), alpha = 0, data = tdat,
                   show.legend = T) +
        geom_hex(aes(x = x, y = y, fill = prop, alpha = count, color = prop),
                 stat = "identity", color = NA, data = pdat2, size = 5) +
        scale_fill_gradientn(colors = gradientn1, name = "% Genic",
                             breaks = leg.brks, labels = leg.labs,
                             limits = c(-1, 1)) +
        scale_color_gradientn(colors = gradientn1, name = "% Genic",
                              breaks = leg.brks, labels = leg.labs,
                              limits = c(-1, 1)) +
        scale_alpha(range = c(0.01, 1), name = "Number of bins", guide = F,
                    trans = scales::trans_new("square",
                                              function(x) {
                                                x^(pow)
                                              },
                                              function(x) x^(1/pow))) +
        coord_cartesian(expand = F,
                        xlim = lim$x,
                        ylim = lim$y) +
        facet_grid(.~ttl) +
        labs(x = sprintf('%s \u25ba',xaxis_label),
             y = sprintf('%s \u25ba',yaxis_label)) +
        annotate('segment', x = -Inf, xend =Inf, y = Inf,
                 yend = Inf, color = 'white') +
        annotate("segment", x=-Inf, xend=Inf,y=-Inf, yend=Inf,color="white") +
        theme(panel.background = element_rect(fill = "black"),
              legend.position = "none",
              panel.grid = element_blank(),
              plot.background = element_blank(),
              legend.background = element_blank(),
              legend.margin = margin(0.015, 0, 0, 0, unit="npc"),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_blank(),
              strip.text = element_text(color = "white"),
              # strip.background = element_rect(fill = "#726866"),
              strip.background = element_rect(fill = 'black'),
              axis.title = element_text(family = "Arial", color = "black"))

      pdat <- kde2d(x = pdat$x, y = pdat$y, n = 100)
      brks <- pretty(c(pdat$z), 30)
      nbrks <- length(brks)
      fac <- diff(brks[1:2]) * nrow(d)/sum(pdat$z)
      b <- isobands(x = pdat$x, y = pdat$y, z = t(pdat$z), brks[1:(nbrks - 1)], brks[2:nbrks])
      bands <- iso_to_sfg(b)
      bdat <- st_sf(level = 1:length(bands), geometry = st_sfc(bands))
      if(!all(st_is_valid(bdat))) {
        bdat <- lwgeom::st_make_valid(bdat)
      }
      bnds <- st_bbox(bdat %>% filter(level > 1))

      bdat$ttl <- title_of_plot
      if(show_legend == FALSE){
        pl <- ggplot() +
          geom_sf(data = bdat, aes(fill = level), color = NA) +
          scale_fill_gradientn(colors = pals::magma(10),
                               name = "# of bins \u25ba") +
          coord_sf(expand = F,
                   xlim = lim$x,
                   ylim = lim$y) +
          labs(x = sprintf('%s \u25ba', xaxis_label),
               y = sprintf('%s \u25ba', yaxis_label)) +
          facet_grid(. ~ ttl) +
          guides(fill = guide_colorbar(barwidth = 3, barheight = 0.5,
                                       title.position = "top", title.hjust = 0.5,
                                       frame.colour = "white", frame.linewidth = 1,
                                       ticks = F)) +
          annotate('segment', x = -Inf, xend =Inf, y = Inf,
                   yend = Inf, color = 'white') +
          annotate("segment", x=-Inf, xend=Inf,y=-Inf, yend=Inf,color="white")+
          theme(panel.background = element_rect(fill = "black"),
                panel.grid = element_blank(),
                plot.background = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                text = element_text(color = "black"),
                legend.title = element_text(color = "white", family = "Arial",
                                            size = 7),
                legend.justification = c(0, 1),
                legend.background = element_blank(),
                # legend.position = c(0.05, 0.95),
                legend.position="none",
                legend.direction = "horizontal",
                legend.text = element_blank(),
                #axis.ticks = element_line(color = "black"),
                #axis.text = element_text(color = "black"),
                plot.title = element_blank(),
                strip.text = element_text(color = "white"),
                #strip.background = element_rect(fill = "#726866"),
                strip.background = element_rect(fill = 'black'),
                axis.title = element_text(family = "Arial"))
      } else{
        pl <- ggplot() +
          geom_sf(data = bdat, aes(fill = level), color = NA) +
          scale_fill_gradientn(colors = pals::magma(10),
                               name = "# of bins \u25ba") +
          coord_sf(expand = F,
                   xlim = lim$x,
                   ylim = lim$y) +
          labs(x = sprintf('%s \u25ba', xaxis_label),
               y = sprintf('%s \u25ba', yaxis_label)) +
          facet_grid(. ~ ttl) +
          guides(fill = guide_colorbar(barwidth = 3, barheight = 0.5,
                                       title.position = "top", title.hjust = 0.5,
                                       frame.colour = "white", frame.linewidth = 1,
                                       ticks = F)) +
          annotate('segment', x = -Inf, xend =Inf, y = Inf,
                   yend = Inf, color = 'white') +
          annotate("segment", x=-Inf, xend=Inf,y=-Inf, yend=Inf,color="white")+
          theme(panel.background = element_rect(fill = "black"),
                panel.grid = element_blank(),
                plot.background = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                text = element_text(color = "black"),
                legend.title = element_text(color = "white", family = "Arial",
                                            size = 7),
                legend.justification = c(0, 1),
                legend.background = element_blank(),
                legend.position = c(0.05, 0.95),
                # legend.position="none",
                legend.direction = "horizontal",
                legend.text = element_blank(),
                #axis.ticks = element_line(color = "black"),
                #axis.text = element_text(color = "black"),
                plot.title = element_blank(),
                strip.text = element_text(color = "white"),
                #strip.background = element_rect(fill = "#726866"),
                strip.background = element_rect(fill = 'black'),
                axis.title = element_text(family = "Arial"))
    }}
    pl$coordinates$aspect <- function(f) { NULL }


    d <- lapply(x[2:1], function(y) {
      findOverlaps(r, y) %>%
        to() %>%
        {y[.]} %>%
        score()
    }) %>%
      bind_cols() %>%
      `names<-`(c('x', 'y')) %>%
      mutate(clu = clus) %>%
      dplyr::filter(x > quantile(x, .01),
                    x < quantile(x, .99),
                    y > quantile(y, .01),
                    y < quantile(y, .99))

    hull <- na.omit(d) %>%
      group_by(clu) %>%
      dplyr::slice(chull(x, y))
    lab <- na.omit(d) %>%
      group_by(clu) %>%
      summarise(x = mean(x),
                y = mean(y))
    cclr <- setNames(tableau20(4)[seq(1,6,2)],
                     c('B', 'A'))

    bdat2 <- bdat
    bdat2$ttl <- 'Density-based clusters'
    d$ttl <- bdat2$ttl[1]
    # gradientn <- paste0('#FFFFFF', as.hexmode(round(seq(0,255, length.out = nbrks))))
    gradientn <- paste0('#FFFFFF', sprintf("%02X", c(0, round(seq(0, 255, length.out = nbrks - 1)))))
    if (show_scales == TRUE){
    pr <- ggplot() +
      geom_sf(data = bdat2, aes(fill = level), color = NA) +
      geom_point_rast(data = d, aes(x = x, y = y, color = clu),
                      alpha = .1, size = .02, raster.dpi = 300,
                      raster.height = 5, raster.width = 5) +
      geom_label_repel(aes(label = clu, x = x, y = y,
                           color = clu), data = lab,
                       show.legend = F,
                       point.padding = NA, box.padding = 0,
                       direction = "x", inherit.aes = F) +
      scale_fill_gradientn(colors = gradientn) +
      scale_color_manual(values = cclr) +
      coord_sf(expand = F,
               xlim = lim$x,
               ylim = lim$y) +
      facet_grid(.~ttl) +
      annotate('segment', x = -Inf, xend =Inf, y = Inf,
               yend = Inf, color = 'white') +
      annotate("segment", x=-Inf, xend=Inf,y=-Inf, yend=Inf,color="white")+
      labs(x = sprintf('%s \u25ba', xaxis_label),
           y = sprintf('%s \u25ba', yaxis_label)) +
      theme(panel.background = element_rect(fill = "black"),
            legend.position = "none",
            panel.grid = element_blank(),
            plot.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(color = "black"),
            legend.background = element_blank(),
            #axis.ticks = element_line(color = "black"),
            #axis.text = element_text(color = "black"),
            plot.title = element_blank(),
            strip.text = element_text(color = "white"),
            #strip.background = element_rect(fill = '#726866'),
            strip.background = element_rect(fill = 'black'),
            legend.margin = margin(0.015, 0, 0, 0, unit="npc"),
            axis.title = element_text(family = "Arial"))
    } else {
      pr <- ggplot() +
        geom_sf(data = bdat2, aes(fill = level), color = NA) +
        geom_point_rast(data = d, aes(x = x, y = y, color = clu),
                        alpha = .1, size = .02, raster.dpi = 300,
                        raster.height = 5, raster.width = 5) +
        geom_label_repel(aes(label = clu, x = x, y = y,
                             color = clu), data = lab,
                         show.legend = F,
                         point.padding = NA, box.padding = 0,
                         direction = "x", inherit.aes = F) +
        scale_fill_gradientn(colors = gradientn) +
        scale_color_manual(values = cclr) +
        coord_sf(expand = F,
                 xlim = lim$x,
                 ylim = lim$y) +
        facet_grid(.~ttl) +
        annotate('segment', x = -Inf, xend =Inf, y = Inf,
                 yend = Inf, color = 'white') +
        annotate("segment", x=-Inf, xend=Inf,y=-Inf, yend=Inf,color="white")+
        labs(x = sprintf('%s \u25ba', xaxis_label),
             y = sprintf('%s \u25ba', yaxis_label)) +
        theme(panel.background = element_rect(fill = "black"),
              legend.position = "none",
              panel.grid = element_blank(),
              plot.background = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              text = element_text(color = "black"),
              legend.background = element_blank(),
              #axis.ticks = element_line(color = "black"),
              #axis.text = element_text(color = "black"),
              plot.title = element_blank(),
              strip.text = element_text(color = "white"),
              #strip.background = element_rect(fill = '#726866'),
              strip.background = element_rect(fill = 'black'),
              legend.margin = margin(0.015, 0, 0, 0, unit="npc"),
              axis.title = element_text(family = "Arial"))
    }
    if(show_legend == FALSE) {
    pr$coordinates$aspect <- function(f) { NULL }
    lo <- c(patchwork::area(t = 1, l = 1, b = 20, r = 20),
            #original parameters:t = 1, l = 1, b = 20, r = 20
            patchwork::area(t = 2, l = 2, b = 8, r = 8))
    #original params: t = 2, l = 2, b = 8, r = 8
    pm_leg <- pm + leg + plot_layout(design = lo)
    wrap_plots(pl, pm, pr, nrow = 1)
    } else{
      pr$coordinates$aspect <- function(f) { NULL }
      lo <- c(patchwork::area(t = 1, l = 1, b = 20, r = 20),
              #original parameters:t = 1, l = 1, b = 20, r = 20
              patchwork::area(t = 2, l = 2, b = 8, r = 8))
      #original params: t = 2, l = 2, b = 8, r = 8
      pm_leg <- pm + leg + plot_layout(design = lo)
      wrap_plots(pl, pm_leg, pr, nrow = 1)
    }
  })

#assign((sprintf('%s.%s.%s.%s.clus.fig',cell_line,control,test,mark)), ps)
save(ps,file=(sprintf('figs/%s.%s.%s.%s%sclus.fig.rda',cell_line,control,test,mark,window_size)))
# ggsave((sprintf('figs/%s.%s.%s.%s%sclus.pdf',cell_line,control,test,mark,window_size)), ps[[cell_line]], height = 2.4, width = 6.2, device = cairo_pdf)
# ggsave((sprintf('figs/%s.%s.%s.%s%sclus.png',cell_line,control,test,mark,window_size)), ps[[cell_line]], height = 2.4, width = 6.2, device = "png",dpi = 600,units = "cm")
ggsave((sprintf('figs/%s.%s.%s.%s%sclus.pdf',cell_line,control,test,mark,window_size)), ps[[cell_line]], height = height_of_figure, width = width_of_figure,dpi = 600,units = "cm", device = cairo_pdf)
ggsave((sprintf('figs/%s.%s.%s.%s%sclus.png',cell_line,control,test,mark,window_size)), ps[[cell_line]], height = height_of_figure, width = width_of_figure, device = "png",dpi = 600,units = "cm")
  } else {
    clus <- case_when(
      overlapsAny(r, cons$A) ~ 'A',
      overlapsAny(r, cons$B) ~ 'B',
      overlapsAny(r, cons$C) ~ 'C',
      T ~ 'NA'
    )
    clus[clus == 'NA'] <- NA

    olap <- tibble(gene = overlapsAny(r, gene),
                   igr = overlapsAny(r, igr)) %>%
      mutate(out = case_when(
        gene & !igr ~ 1,
        !gene & igr ~ -1,
        TRUE ~ 0
      )) %>%
      pull(out)

    lineclr <- "black"
    horz <- F
    gradientn1 <- brewer.rdylbu(50)
    cramp <- colorRampPalette(c("#000000ff","#ffffff00"), alpha = T)(5)

    leg.brks <- seq(-1, 1, length.out = 19)[seq(2, 18, by = 2)]
    leg.labs <- c(sprintf('Genic\u25bc'), rep('', 3), '50%',
                  rep('', 3), sprintf('Genic\u25b2'))
    len <- 9
    pal <- brewer.rdylbu(len)
    cmat <- seq(0, 255, length.out = len + 1) %>%
      {.[-1]} %>%
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
                data = cmat) +
      scale_fill_identity() +
      labs(x = "# of bins \u25ba",
           y = "% genic \u25ba") +
      coord_fixed(expand = F) +
      theme(panel.background = element_blank(),
            plot.background = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_blank(),
            panel.border = element_rect(colour = "white", fill = NA, size = 0.5),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_text(family = 'Arial',
                                      color = "white",
                                      size = 7))

    ps <- split(d, s$line) %>%
      lapply(function(x) {
        pdat <- lapply(x[2:1], function(y) {
          findOverlaps(r, y) %>%
            to() %>%
            {y[.]} %>%
            score()
        }) %>%
          bind_cols() %>%
          `names<-`(c('x', 'y')) %>%
          mutate(r = olap) %>%
          dplyr::filter(x > quantile(x, .01),
                        x < quantile(x, .99),
                        y > quantile(y, .01),
                        y < quantile(y, .99))


        hex <- hexbin(pdat$x, pdat$y, xbins = 75, IDs = T)
        pdat$cell <- hex@cID
        hex <- data.frame(hcell2xy(hex),
                          cell = hex@cell,
                          count = hex@count)

        t1 <- 'Intergenic vs genic ratio'
        pdat2 <- pdat %>%
          group_by(cell) %>%
          summarise(prop = mean(r, na.rm = T)) %>%
          ungroup %>%
          right_join(hex, by = "cell") %>%
          mutate(logcount = log10(count),
                 ttl = t1)
        lim <- data.frame(x = c(min(pdat2$x[pdat2$count > 10]),
                                max(pdat2$x[pdat2$count > 10])) %>%
                            scales::expand_range(mul = .05),
                          y = c(min(pdat2$y[pdat2$count > 10]),
                                max(pdat2$y[pdat2$count > 10])) %>%
                            scales::expand_range(mul = .05))
        yrang <- diff(lim$y)
        tdat <- data.frame(x = rep(median(pdat$x), 2),
                           y = rep(median(pdat$y), 2),
                           c = c(0,1),
                           ttl = t1)
        pow <- pow
        if (show_scales == TRUE){
        pm <- ggplot() +
          geom_point(aes(x = x, y = y, color = c), alpha = 0, data = tdat,
                     show.legend = T) +
          geom_hex(aes(x = x, y = y, fill = prop, alpha = count, color = prop),
                   stat = "identity", color = NA, data = pdat2,size = 5) +
          scale_fill_gradientn(colors = gradientn1, name = "% Genic",
                               breaks = leg.brks, labels = leg.labs,
                               limits = c(-1, 1)) +
          scale_color_gradientn(colors = gradientn1, name = "% Genic",
                                breaks = leg.brks, labels = leg.labs,
                                limits = c(-1, 1)) +
          scale_alpha(range = c(0.01, 1), name = "Number of bins", guide = F,
                      trans = scales::trans_new("square",
                                                function(x) {
                                                  x^(pow)
                                                },
                                                function(x) x^(1/pow))) +
          coord_cartesian(expand = F,
                          xlim = lim$x,
                          ylim = lim$y) +
          facet_grid(.~ttl) +
          labs(x = sprintf("%s \u25ba",xaxis_label),
               y = sprintf("%s \u25ba",yaxis_label)) +
          annotate('segment', x = -Inf, xend =Inf, y = Inf,
                   yend = Inf, color = 'white') +
          theme(panel.background = element_rect(fill = "black"),
                legend.position = "none",
                panel.grid = element_blank(),
                plot.background = element_blank(),
                legend.background = element_blank(),
                legend.margin = margin(0.015, 0, 0, 0, unit="npc"),
                # axis.text = element_blank(),
                # axis.ticks = element_blank(),
                axis.line = element_blank(),
                plot.title = element_blank(),
                strip.text = element_text(color = "white"),
                # strip.background = element_rect(fill = "#726866"),
                strip.background = element_rect(fill = 'black'),
                axis.title = element_text(family = "Arial", color = "black"))

        pdat <- kde2d(x = pdat$x, y = pdat$y, n = 100)
        brks <- pretty(c(pdat$z), 30)
        nbrks <- length(brks)
        fac <- diff(brks[1:2]) * nrow(d)/sum(pdat$z)
        b <- isobands(x = pdat$x, y = pdat$y, z = t(pdat$z), brks[1:(nbrks - 1)], brks[2:nbrks])
        bands <- iso_to_sfg(b)
        bdat <- st_sf(level = 1:length(bands), geometry = st_sfc(bands))
        if(!all(st_is_valid(bdat))) {
          bdat <- lwgeom::st_make_valid(bdat)
        }
        bnds <- st_bbox(bdat %>% filter(level > 1))

        bdat$ttl <- title_of_plot
        if(show_legend == FALSE){
        pl <- ggplot() +
          geom_sf(data = bdat, aes(fill = level), color = NA) +
          scale_fill_gradientn(colors = pals::magma(10),
                               name = "# of bins \u25ba") +
          coord_sf(expand = F,
                   xlim = lim$x,
                   ylim = lim$y) +
          labs(x = sprintf("%s \u25ba",xaxis_label),
               y = sprintf("%s \u25ba",yaxis_label)) +
          facet_grid(. ~ ttl) +
          guides(fill = guide_colorbar(barwidth = 3, barheight = 0.5,
                                       title.position = "top", title.hjust = 0.5,
                                       frame.colour = "white", frame.linewidth = 1,
                                       ticks = F)) +
          annotate('segment', x = -Inf, xend =Inf, y = Inf,
                   yend = Inf, color = 'white') +
          theme(panel.background = element_rect(fill = "black"),
                panel.grid = element_blank(),
                plot.background = element_blank(),
                # axis.text = element_blank(),
                # axis.ticks = element_blank(),
                text = element_text(color = "black"),
                legend.title = element_text(color = "white", family = "Arial",
                                            size = 7),
                legend.justification = c(0, 1),
                legend.background = element_blank(),
                # legend.position = c(0.05, 0.95),
                legend.position = "none",
                legend.direction = "horizontal",
                legend.text = element_blank(),
                #axis.ticks = element_line(color = "black"),
                #axis.text = element_text(color = "black"),
                plot.title = element_blank(),
                strip.text = element_text(color = "white"),
                #strip.background = element_rect(fill = "#726866"),
                strip.background = element_rect(fill = 'black'),
                axis.title = element_text(family = "Arial"))
        } else {
          pl <- ggplot() +
            geom_sf(data = bdat, aes(fill = level), color = NA) +
            scale_fill_gradientn(colors = pals::magma(10),
                                 name = "# of bins \u25ba") +
            coord_sf(expand = F,
                     xlim = lim$x,
                     ylim = lim$y) +
            labs(x = sprintf("%s \u25ba",xaxis_label),
                 y = sprintf("%s \u25ba",yaxis_label)) +
            facet_grid(. ~ ttl) +
            guides(fill = guide_colorbar(barwidth = 3, barheight = 0.5,
                                         title.position = "top", title.hjust = 0.5,
                                         frame.colour = "white", frame.linewidth = 1,
                                         ticks = F)) +
            annotate('segment', x = -Inf, xend =Inf, y = Inf,
                     yend = Inf, color = 'white') +
            theme(panel.background = element_rect(fill = "black"),
                  panel.grid = element_blank(),
                  plot.background = element_blank(),
                  # axis.text = element_blank(),
                  # axis.ticks = element_blank(),
                  text = element_text(color = "black"),
                  legend.title = element_text(color = "white", family = "Arial",
                                              size = 7),
                  legend.justification = c(0, 1),
                  legend.background = element_blank(),
                  legend.position = c(0.05, 0.95),
                  # legend.position = "none",
                  legend.direction = "horizontal",
                  legend.text = element_blank(),
                  #axis.ticks = element_line(color = "black"),
                  #axis.text = element_text(color = "black"),
                  plot.title = element_blank(),
                  strip.text = element_text(color = "white"),
                  #strip.background = element_rect(fill = "#726866"),
                  strip.background = element_rect(fill = 'black'),
                  axis.title = element_text(family = "Arial"))
        }
        pl$coordinates$aspect <- function(f) { NULL }


        d <- lapply(x[2:1], function(y) {
          findOverlaps(r, y) %>%
            to() %>%
            {y[.]} %>%
            score()
        }) %>%
          bind_cols() %>%
          `names<-`(c('x', 'y')) %>%
          mutate(clu = clus) %>%
          dplyr::filter(x > quantile(x, .01),
                        x < quantile(x, .99),
                        y > quantile(y, .01),
                        y < quantile(y, .99))

        hull <- na.omit(d) %>%
          group_by(clu) %>%
          dplyr::slice(chull(x, y))
        lab <- na.omit(d) %>%
          group_by(clu) %>%
          summarise(x = mean(x),
                    y = mean(y))
        cclr <- setNames(tableau20(6)[seq(1,6,2)],
                         c('C', 'B', 'A'))

        bdat2 <- bdat
        bdat2$ttl <- 'Density-based clusters'
        d$ttl <- bdat2$ttl[1]
        # gradientn <- paste0('#ffffff', as.hexmode(round(seq(0,255, length.out = nbrks))))
        gradientn <- paste0('#FFFFFF', sprintf("%02X", c(0, round(seq(0, 255, length.out = nbrks - 1)))))
        pr <- ggplot() +
          geom_sf(data = bdat2, aes(fill = level), color = NA) +
          geom_point_rast(data = d, aes(x = x, y = y, color = clu),
                          alpha = .1, size = .02, raster.dpi = 300,
                          raster.height = 5, raster.width = 5) +
          geom_label_repel(aes(label = clu, x = x, y = y,
                               color = clu), data = lab,
                           show.legend = F,
                           point.padding = NA, box.padding = 0,
                           direction = "x", inherit.aes = F) +
          scale_fill_gradientn(colors = gradientn) +
          scale_color_manual(values = cclr) +
          coord_sf(expand = F,
                   xlim = lim$x,
                   ylim = lim$y) +
          facet_grid(.~ttl) +
          annotate('segment', x = -Inf, xend =Inf, y = Inf,
                   yend = Inf, color = 'white') +
          labs(x = sprintf("%s \u25ba", xaxis_label),
               y = sprintf("%s \u25ba",yaxis_label)) +
          theme(panel.background = element_rect(fill = "black"),
                legend.position = "none",
                panel.grid = element_blank(),
                plot.background = element_blank(),
                # axis.text = element_blank(),
                # axis.ticks = element_blank(),
                text = element_text(color = "black"),
                legend.background = element_blank(),
                #axis.ticks = element_line(color = "black"),
                #axis.text = element_text(color = "black"),
                plot.title = element_blank(),
                strip.text = element_text(color = "white"),
                #strip.background = element_rect(fill = '#726866'),
                strip.background = element_rect(fill = 'black'),
                legend.margin = margin(0.015, 0, 0, 0, unit="npc"),
                axis.title = element_text(family = "Arial"))
        if(show_legend == FALSE) {
          pr$coordinates$aspect <- function(f) { NULL }
          lo <- c(patchwork::area(t = 1, l = 1, b = 20, r = 20),
                  #original parameters:t = 1, l = 1, b = 20, r = 20
                  patchwork::area(t = 2, l = 2, b = 8, r = 8))
          #original params: t = 2, l = 2, b = 8, r = 8
          pm_leg <- pm + leg + plot_layout(design = lo)
          wrap_plots(pl, pm, pr, nrow = 1)
        } else{
          pr$coordinates$aspect <- function(f) { NULL }
          lo <- c(patchwork::area(t = 1, l = 1, b = 20, r = 20),
                  #original parameters:t = 1, l = 1, b = 20, r = 20
                  patchwork::area(t = 2, l = 2, b = 8, r = 8))
          #original params: t = 2, l = 2, b = 8, r = 8
          pm_leg <- pm + leg + plot_layout(design = lo)
          wrap_plots(pl, pm_leg, pr, nrow = 1)
        }
        } else{
          pm <- ggplot() +
            geom_point(aes(x = x, y = y, color = c), alpha = 0, data = tdat,
                       show.legend = T) +
            geom_hex(aes(x = x, y = y, fill = prop, alpha = count, color = prop),
                     stat = "identity", color = NA, data = pdat2,size = 5) +
            scale_fill_gradientn(colors = gradientn1, name = "% Genic",
                                 breaks = leg.brks, labels = leg.labs,
                                 limits = c(-1, 1)) +
            scale_color_gradientn(colors = gradientn1, name = "% Genic",
                                  breaks = leg.brks, labels = leg.labs,
                                  limits = c(-1, 1)) +
            scale_alpha(range = c(0.01, 1), name = "Number of bins", guide = F,
                        trans = scales::trans_new("square",
                                                  function(x) {
                                                    x^(pow)
                                                  },
                                                  function(x) x^(1/pow))) +
            coord_cartesian(expand = F,
                            xlim = lim$x,
                            ylim = lim$y) +
            facet_grid(.~ttl) +
            labs(x = sprintf("%s \u25ba",xaxis_label),
                 y = sprintf("%s \u25ba",yaxis_label)) +
            annotate('segment', x = -Inf, xend =Inf, y = Inf,
                     yend = Inf, color = 'white') +
            theme(panel.background = element_rect(fill = "black"),
                  legend.position = "none",
                  panel.grid = element_blank(),
                  plot.background = element_blank(),
                  legend.background = element_blank(),
                  legend.margin = margin(0.015, 0, 0, 0, unit="npc"),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.line = element_blank(),
                  plot.title = element_blank(),
                  strip.text = element_text(color = "white"),
                  # strip.background = element_rect(fill = "#726866"),
                  strip.background = element_rect(fill = 'black'),
                  axis.title = element_text(family = "Arial", color = "black"))

          pdat <- kde2d(x = pdat$x, y = pdat$y, n = 100)
          brks <- pretty(c(pdat$z), 30)
          nbrks <- length(brks)
          fac <- diff(brks[1:2]) * nrow(d)/sum(pdat$z)
          b <- isobands(x = pdat$x, y = pdat$y, z = t(pdat$z), brks[1:(nbrks - 1)], brks[2:nbrks])
          bands <- iso_to_sfg(b)
          bdat <- st_sf(level = 1:length(bands), geometry = st_sfc(bands))
          if(!all(st_is_valid(bdat))) {
            bdat <- lwgeom::st_make_valid(bdat)
          }
          bnds <- st_bbox(bdat %>% filter(level > 1))

          bdat$ttl <- title_of_plot
          if(show_legend == FALSE){
            pl <- ggplot() +
              geom_sf(data = bdat, aes(fill = level), color = NA) +
              scale_fill_gradientn(colors = pals::magma(10),
                                   name = "# of bins \u25ba") +
              coord_sf(expand = F,
                       xlim = lim$x,
                       ylim = lim$y) +
              labs(x = sprintf("%s \u25ba",xaxis_label),
                   y = sprintf("%s \u25ba",yaxis_label)) +
              facet_grid(. ~ ttl) +
              guides(fill = guide_colorbar(barwidth = 3, barheight = 0.5,
                                           title.position = "top", title.hjust = 0.5,
                                           frame.colour = "white", frame.linewidth = 1,
                                           ticks = F)) +
              annotate('segment', x = -Inf, xend =Inf, y = Inf,
                       yend = Inf, color = 'white') +
              theme(panel.background = element_rect(fill = "black"),
                    panel.grid = element_blank(),
                    plot.background = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    text = element_text(color = "black"),
                    legend.title = element_text(color = "white", family = "Arial",
                                                size = 7),
                    legend.justification = c(0, 1),
                    legend.background = element_blank(),
                    # legend.position = c(0.05, 0.95),
                    legend.position = "none",
                    legend.direction = "horizontal",
                    legend.text = element_blank(),
                    #axis.ticks = element_line(color = "black"),
                    #axis.text = element_text(color = "black"),
                    plot.title = element_blank(),
                    strip.text = element_text(color = "white"),
                    #strip.background = element_rect(fill = "#726866"),
                    strip.background = element_rect(fill = 'black'),
                    axis.title = element_text(family = "Arial"))
          } else {
            pl <- ggplot() +
              geom_sf(data = bdat, aes(fill = level), color = NA) +
              scale_fill_gradientn(colors = pals::magma(10),
                                   name = "# of bins \u25ba") +
              coord_sf(expand = F,
                       xlim = lim$x,
                       ylim = lim$y) +
              labs(x = sprintf("%s \u25ba",xaxis_label),
                   y = sprintf("%s \u25ba",yaxis_label)) +
              facet_grid(. ~ ttl) +
              guides(fill = guide_colorbar(barwidth = 3, barheight = 0.5,
                                           title.position = "top", title.hjust = 0.5,
                                           frame.colour = "white", frame.linewidth = 1,
                                           ticks = F)) +
              annotate('segment', x = -Inf, xend =Inf, y = Inf,
                       yend = Inf, color = 'white') +
              theme(panel.background = element_rect(fill = "black"),
                    panel.grid = element_blank(),
                    plot.background = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    text = element_text(color = "black"),
                    legend.title = element_text(color = "white", family = "Arial",
                                                size = 7),
                    legend.justification = c(0, 1),
                    legend.background = element_blank(),
                    legend.position = c(0.05, 0.95),
                    # legend.position = "none",
                    legend.direction = "horizontal",
                    legend.text = element_blank(),
                    #axis.ticks = element_line(color = "black"),
                    #axis.text = element_text(color = "black"),
                    plot.title = element_blank(),
                    strip.text = element_text(color = "white"),
                    #strip.background = element_rect(fill = "#726866"),
                    strip.background = element_rect(fill = 'black'),
                    axis.title = element_text(family = "Arial"))
          }
          pl$coordinates$aspect <- function(f) { NULL }


          d <- lapply(x[2:1], function(y) {
            findOverlaps(r, y) %>%
              to() %>%
              {y[.]} %>%
              score()
          }) %>%
            bind_cols() %>%
            `names<-`(c('x', 'y')) %>%
            mutate(clu = clus) %>%
            dplyr::filter(x > quantile(x, .01),
                          x < quantile(x, .99),
                          y > quantile(y, .01),
                          y < quantile(y, .99))

          hull <- na.omit(d) %>%
            group_by(clu) %>%
            dplyr::slice(chull(x, y))
          lab <- na.omit(d) %>%
            group_by(clu) %>%
            summarise(x = mean(x),
                      y = mean(y))
          cclr <- setNames(tableau20(6)[seq(1,6,2)],
                           c('C', 'B', 'A'))

          bdat2 <- bdat
          bdat2$ttl <- 'Density-based clusters'
          d$ttl <- bdat2$ttl[1]
          # gradientn <- paste0('#ffffff', as.hexmode(round(seq(0,255, length.out = nbrks))))
          gradientn <- paste0('#FFFFFF', sprintf("%02X", c(0, round(seq(0, 255, length.out = nbrks - 1)))))
          pr <- ggplot() +
            geom_sf(data = bdat2, aes(fill = level), color = NA) +
            geom_point_rast(data = d, aes(x = x, y = y, color = clu),
                            alpha = .1, size = .02, raster.dpi = 300,
                            raster.height = 5, raster.width = 5) +
            geom_label_repel(aes(label = clu, x = x, y = y,
                                 color = clu), data = lab,
                             show.legend = F,
                             point.padding = NA, box.padding = 0,
                             direction = "x", inherit.aes = F) +
            scale_fill_gradientn(colors = gradientn) +
            scale_color_manual(values = cclr) +
            coord_sf(expand = F,
                     xlim = lim$x,
                     ylim = lim$y) +
            facet_grid(.~ttl) +
            annotate('segment', x = -Inf, xend =Inf, y = Inf,
                     yend = Inf, color = 'white') +
            labs(x = sprintf("%s \u25ba", xaxis_label),
                 y = sprintf("%s \u25ba",yaxis_label)) +
            theme(panel.background = element_rect(fill = "black"),
                  legend.position = "none",
                  panel.grid = element_blank(),
                  plot.background = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  text = element_text(color = "black"),
                  legend.background = element_blank(),
                  #axis.ticks = element_line(color = "black"),
                  #axis.text = element_text(color = "black"),
                  plot.title = element_blank(),
                  strip.text = element_text(color = "white"),
                  #strip.background = element_rect(fill = '#726866'),
                  strip.background = element_rect(fill = 'black'),
                  legend.margin = margin(0.015, 0, 0, 0, unit="npc"),
                  axis.title = element_text(family = "Arial"))
          if(show_legend == FALSE) {
            pr$coordinates$aspect <- function(f) { NULL }
            lo <- c(patchwork::area(t = 1, l = 1, b = 20, r = 20),
                    #original parameters:t = 1, l = 1, b = 20, r = 20
                    patchwork::area(t = 2, l = 2, b = 8, r = 8))
            #original params: t = 2, l = 2, b = 8, r = 8
            pm_leg <- pm + leg + plot_layout(design = lo)
            wrap_plots(pl, pm, pr, nrow = 1)
          } else{
            pr$coordinates$aspect <- function(f) { NULL }
            lo <- c(patchwork::area(t = 1, l = 1, b = 20, r = 20),
                    #original parameters:t = 1, l = 1, b = 20, r = 20
                    patchwork::area(t = 2, l = 2, b = 8, r = 8))
            #original params: t = 2, l = 2, b = 8, r = 8
            pm_leg <- pm + leg + plot_layout(design = lo)
            wrap_plots(pl, pm_leg, pr, nrow = 1)
          }
        }
      })
    #assign((sprintf('%s.%s.%s.%s.clus.fig',cell_line,control,test,mark)), ps)
    save(ps,file=(sprintf('plots/%s.%s.%s.%s%sclus.fig.rda',cell_line,control,test,mark,window_size)))
#     ggsave((sprintf('figs/%s.%s.%s.%s%sclus.pdf',cell_line,control,test,mark,window_size)), ps[[cell_line]], height = 2.4, width = 6.2, device = cairo_pdf)}
# ggsave((sprintf('figs/%s.%s.%s.%s%sclus.png',cell_line,control,test,mark,window_size)), ps[[cell_line]], height = 2.4, width = 6.2, device = "png")
    ggsave((sprintf('figs/%s.%s.%s.%s%sclus.pdf',cell_line,control,test,mark,window_size)), ps[[cell_line]], height = 2.4, width = 6.2, device = cairo_pdf)}
ggsave((sprintf('figs/%s.%s.%s.%s%sclus.png',cell_line,control,test,mark,window_size)), ps[[cell_line]], height = height_of_figure, width = width_of_figure, device = "png",dpi = 600,units = "cm")
}
