#!/usr/bin/env Rscript
#+ message = FALSE, warning = FALSE
# library(data.table)
# library(tidyverse)
# library(rtracklayer)
# library(matrixStats)
# library(isoband)
# library(sf)
# library(ggrepel)
# library(viridis)
# library(MASS)
# library(lwgeom)
# library(hexbin)
# library(pals)
# library(patchwork)
# library(gdata)
#' Title
#'
#' @param path_for_norm_bw Directory
#' @param out_dir
#' @param gene
#' @param intergenic
#' @param cell_line
#' @param histone_mark
#' @param treated_samp_label
#' @param control_label
#' @param window_size
#' @param title_of_plot
#' @param max_x
#' @param max_y
#' @param pow
#' @param xaxis_label
#' @param yaxis_label
#' @param min_x
#' @param min_y
#' @param show_scales
#'
#' @return
#' @export
#'
#' @examples
genic_intergenic_scatterplot <- function(path_for_norm_bw,
                                          out_dir,
                                          gene,
                                          intergenic,
                                          cell_line,
                                          histone_mark,
                                          treated_samp_label,
                                          control_label,
                                          window_size,
                                          title_of_plot,
                                          max_x,
                                          max_y,
                                          pow=NULL,
                                          xaxis_label,
                                          yaxis_label,
                                          min_x=NULL,
                                          min_y=NULL,
                                          show_scales=TRUE){
# setwd(wd)
# getwd()
path <- paste0(path_for_norm_bw)
out_dir <- paste0(out_dir)
gene <- import.bed(gene)
igr <- import.bed(intergenic)
cell_line=paste0(cell_line)
control_label=paste0(control_label)
treated_samp_label=paste0(treated_samp_label)
mark=paste0(histone_mark)
window_size=paste0(".",window_size,"kb.")
# x and y axis labels
xaxis_label=paste0(xaxis_label)
yaxis_label=paste0(yaxis_label)
# parameters for max x and y
max_x = as.numeric(max_x)
max_y = as.numeric(max_y)
# new titles
title_of_plot=paste0(title_of_plot)
# list files
s <- list.files(path = path, pattern = window_size,full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c("line","samp","mark",NA,NA),'\\.',F) %>%
  mutate(f = file.path(path, f)) %>%
  dplyr::filter(line == cell_line) %>%
  dplyr::filter(samp == control_label | samp == treated_samp_label)
# cell_line <- as.character(s$line[[1]])
# NEED TO ORDER SAMPLE FIRST THEN PARENTAL #
odr <- c(treated_samp_label,control_label)
s <- s %>%
  dplyr::slice(match(odr,samp))
#
d <- deframe(s[,c('samp', 'f')]) %>%
  lapply(import.bw)
#
#
d[[control_label]] <- subsetByOverlaps(d[[control_label]],d[[treated_samp_label]])
d[[treated_samp_label]]<- subsetByOverlaps(d[[treated_samp_label]],d[[control_label]])
#
lapply(d,length)
#
#
r <- lapply(d, function(y) y[y$score != 0]) %>%
  Reduce(function(a, b) a[overlapsAny(a, b)], .) %>%
  granges()
# overlaps
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
        panel.border = element_rect(colour = "white", fill = NA, size = 0.1),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Arial',
                                  color = "black",
                                  size = 10))
leg
setwd(out_dir)
ggsave(filename = "legend.png", plot=leg, height = 1, width = 1,device = "png",dpi=600,bg = "white")
# pow
if(is.null(pow)){
  pow = as.numeric(1.25)} else {
    pow = as.numeric(paste0(pow))
  }
# make legend
# cfunc <- colorRampPalette(c("#93003a","#ad1042","#c32a49",
#                             "#d4444d","#e35e50","#ee7851",
#                             "#f59350","#faad4c","#fbc846",
#                             "#f9e43c","#f3ff2c","#cbed92",
#                             "#b2daa9","#9cc6b3","#88b2b8",
#                             "#759fb8","#638bb6","#5178b2",
#                             "#3e66ac","#2854a5","#00429d"))
# lineclr <- "black"
# horz <- F
# gradientn1 <- cfunc(21)
# cramp <- colorRampPalette(c("#000000ff","#ffffff00"), alpha = T)(5)
# leg.brks <- seq(-1, 1, length.out = 19)[seq(2, 18, by = 2)]
# leg.labs <- c(sprintf('Genic\u25bc'), rep('', 3), '50%',
#               rep('', 3), sprintf('Genic\u25b2'))
# len <- 9
# pal <- cfunc(len)
# cmat <- seq(0, 255, length.out = len + 1) %>%
# {.[-1]} %>%
#   round() %>%
#   as.hexmode() %>%
#   format(width = 2, upper.case = T) %>%
#   lapply(function(x) {
#     paste0(pal, x)
#   }) %>%
#   do.call(cbind, .) %>%
#   data.frame() %>%
#   `colnames<-`(1:len) %>%
#   mutate(clr = 1:dplyr::n()) %>%
#   reshape2::melt(id.vars = "clr", variable.name = "opa") %>%
#   mutate(opa = as.integer(opa))
# leg <- ggplot() +
#   geom_tile(aes(x = opa, y = clr, fill = value),
#             data = cmat) +
#   scale_fill_identity() +
#   labs(x = "# of bins \u25ba",
#        y = "% genic \u25ba") +
#   coord_fixed(expand = F) +
#   theme(panel.background = element_blank(),
#         plot.background = element_blank(),
#         panel.grid = element_blank(),
#         axis.line = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_text(family = 'Arial',
#                                   color = "black"))
# leg

#### plot ####
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


    hex <- hexbin(pdat$x, pdat$y, xbins = 100, IDs = T)
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
    # lim <- data.frame(x = c(min(pdat2$x[pdat2$count > 5]),
    #                         max(pdat2$x[pdat2$count > 5])) %>%
    #                     scales::expand_range(mul = .05),
    #                   y = c(min(pdat2$y[pdat2$count > 5]),
    #                         max(pdat2$y[pdat2$count > 5])) %>%
    #                     scales::expand_range(mul = .05))
    if(is.null(min_x)){
      min_x = as.numeric(min(pdat2$x[pdat2$count > 10]))} else {
        min_x = as.numeric(paste0(min_x))
      }
    if(is.null(min_y)){
      min_y = as.numeric(min(pdat2$y[pdat2$count > 10]))} else {
        min_y = as.numeric(paste0(min_y))
      }
    lim <- data.frame(x = c(min_x,
                            max_x) %>%
                        scales::expand_range(mul = .05),
                      y = c(min_y,
                            max_y) %>%
                        scales::expand_range(mul = .05))
    yrang <- diff(lim$y)
    tdat <- data.frame(x = rep(median(pdat$x), 2),
                       y = rep(median(pdat$y), 2),
                       c = c(0,1),
                       ttl = t1)
    # what is pow?
    pow <- pow
    if (show_scales == TRUE) {
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
      scale_alpha(range = c(0.01, 5), name = "Number of bins", guide = F,
                  trans = scales::trans_new("square",
                                            function(x) {
                                              x^(pow)
                                            },
                                            function(x) x^(1/pow))) +
      coord_cartesian(expand = F,
                      xlim = lim$x,
                      ylim = lim$y) +
      # facet_grid(.~ttl) +
      labs(x = sprintf("%s \u25ba",xaxis_label),
           y = sprintf("%s \u25ba",yaxis_label),title=title_of_plot) +
      # annotate('segment', x = -Inf, xend =Inf, y = Inf,
      #          yend = Inf, color = 'white') +
      annotate("segment", x=-Inf, xend=Inf,y=-Inf, yend=Inf,color="gray") +
      theme(panel.background = element_rect(fill = "white"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            legend.position = "none",
            panel.grid = element_blank(),
            plot.background = element_blank(),
            legend.background = element_blank(),
            legend.margin = margin(0.015, 0, 0, 0, unit="npc"),
            # axis.text = element_blank(),
            # axis.ticks = element_blank(),
            # axis.line = element_blank(),
            plot.title = element_text(hjust = 0.5,color = "black",size=12,family="Arial"),
            strip.text = element_text(color = "black"),
            # strip.background = element_rect(fill = "#726866"),
            strip.background = element_rect(fill = 'black'),
            axis.title = element_text(family = "Arial", color = "black"))
    } else {
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
        scale_alpha(range = c(0.01, 5), name = "Number of bins", guide = F,
                    trans = scales::trans_new("square",
                                              function(x) {
                                                x^(pow)
                                              },
                                              function(x) x^(1/pow))) +
        coord_cartesian(expand = F,
                        xlim = lim$x,
                        ylim = lim$y) +
        # facet_grid(.~ttl) +
        labs(x = sprintf("%s \u25ba",xaxis_label),
             y = sprintf("%s \u25ba",yaxis_label),title=title_of_plot) +
        # annotate('segment', x = -Inf, xend =Inf, y = Inf,
        #          yend = Inf, color = 'white') +
        annotate("segment", x=-Inf, xend=Inf,y=-Inf, yend=Inf,color="gray") +
        theme(panel.background = element_rect(fill = "white"),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              legend.position = "none",
              panel.grid = element_blank(),
              plot.background = element_blank(),
              legend.background = element_blank(),
              legend.margin = margin(0.015, 0, 0, 0, unit="npc"),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(hjust = 0.5,color = "black",size=12,family="Arial"),
              strip.text = element_text(color = "black"),
              # strip.background = element_rect(fill = "#726866"),
              strip.background = element_rect(fill = 'black'),
              axis.title = element_text(family = "Arial", color = "black"))
    }
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

    # pl <- ggplot() +
    #   geom_sf(data = bdat, aes(fill = level), color = NA) +
    #   scale_fill_gradientn(colors = pals::magma(10),
    #                        name = "# of bins \u25ba") +
    #   coord_sf(expand = F,
    #            xlim = lim$x,
    #            ylim = lim$y) +
    #   labs(x = sprintf("%s \u25ba",control_label),
    #        y = sprintf("%s \u25ba",treated_samp_label)) +
    #   facet_grid(. ~ ttl) +
    #   guides(fill = guide_colorbar(barwidth = 3, barheight = 0.5,
    #                                title.position = "top", title.hjust = 0.5,
    #                                frame.colour = "white", frame.linewidth = 1,
    #                                ticks = F)) +
    #   # annotate('segment', x = -Inf, xend =Inf, y = Inf,
    #   #          yend = Inf, color = 'white') +
    #   theme(panel.background = element_rect(fill = "black"),
    #         panel.grid = element_blank(),
    #         plot.background = element_blank(),
    #         axis.text = element_blank(),
    #         axis.ticks = element_blank(),
    #         text = element_text(color = "black"),
    #         legend.title = element_text(color = "white", family = "Arial",
    #                                     size = 7),
    #         legend.justification = c(0, 1),
    #         legend.background = element_blank(),
    #         legend.position = c(0.05, 0.95),
    #         legend.direction = "horizontal",
    #         legend.text = element_blank(),
    #         #axis.ticks = element_line(color = "black"),
    #         #axis.text = element_text(color = "black"),
    #         plot.title = element_blank(),
    #         strip.text = element_text(color = "white"),
    #         #strip.background = element_rect(fill = "#726866"),
    #         strip.background = element_rect(fill = 'black'),
    #         axis.title = element_text(family = "Arial"))
    # pl$coordinates$aspect <- function(f) { NULL }
    # lo <- c(area(t = 1, l = 1, b = 20, r = 20),
    #         area(t = 2, l = 2, b = 8, r = 8))
    # pm_leg <- pm + leg + plot_layout(design = lo)
    # wrap_plots(pl, pm_leg, nrow = 1)
    wrap_plots(pm,nrow = 1)
  })
ps[[cell_line]]
ggsave((sprintf('%s.%s.%s.%s%sscatterplot.png',cell_line,control_label,treated_samp_label,mark,window_size)), ps[[cell_line]], height = 4, width = 4, device = "png",dpi=600)
# library(gdata)
# # rename object using library gdata
# mv(from = "ps", to = paste0(sprintf('%s.%s.%s.%s%sscatterplot',cell_line,control_label,treated_samp_label,mark,window_size)))
# #assign((sprintf('%s.%s.%s.%s.clus.fig',cell_line,control_label,treated_samp_label,mark)), ps)
# save(list=paste0(sprintf('%s.%s.%s.%s%sscatterplot',cell_line,control_label,treated_samp_label,mark,window_size)),file=(sprintf('%s.%s.%s.%s%sscatterplot_data.rda',cell_line,control_label,treated_samp_label,mark,window_size)))
# # save legend as data
# save(leg,file = "legend.RData")
print("Generated genic/intergenic scatterplot!")
}
