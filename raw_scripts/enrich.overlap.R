#!/usr/bin/env Rscript
#+ message = FALSE, warning = FALSE
library(pacman)
pacman::p_load(argparse,tidyverse,rtracklayer,GenomicRanges,LOLA,highcharter)
#+
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-wd",help="Working directory")
#
parser$add_argument("-c","--control",type="character",help="Control sample")
parser$add_argument("-t","--test",type="character",help="Test sample")
parser$add_argument("-l","--line",type="character",help="Cell line")
parser$add_argument("-m","--mark",type="character",help="Histone mark")
#
parser$add_argument("-cl","--cluster",type="character",help="Cluster for enrichment")
#
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
args <- parser$parse_args()
#
if ( args$verbose ) { 
  write("Start generating overlap plots...\n", stderr()) 
}
#need to set wd 
setwd(args$wd)
getwd()
# load params
cell_line=paste0(args$line)
control=paste0(args$control)
test=paste0(args$test)
mark=paste0(args$mark)
cluster=paste0(args$cluster)
#load gene and igr bed
g <- import.bed('regionDB/hg38/ensembl/regions/gene.bed')
ig <- import.bed('regionDB/hg38/ensembl/regions/intergenic.bed')
#load cons
load((sprintf('data/cons/cons.%s.%s.%s.%s.rda',cell_line,control,test,mark)))
# background
uni <- import.bed(sprintf('data/pooled.bed/%s.%s.%s.%s.pooled.bed',cell_line,control,test,mark)) 
uni_igr <- uni[overlapsAny(uni, ig) & !overlapsAny(uni, g)]
uni_g <- uni[overlapsAny(uni, g) & !overlapsAny(uni, ig)]
# user set
qSet <- cons[[cluster]]
qSet.igr <- qSet[overlapsAny(qSet, ig) & !overlapsAny(qSet, g)]
qSet.g <- qSet[overlapsAny(qSet, g) & !overlapsAny(qSet, ig)]
# load functions
signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
         symbols = c("****", "***", "**", "*", "ns"))
}
# collection of functional annotations
ensemblDB <- loadRegionDB('regionDB/mm10',collections = "ensembl")
ccreDB <- loadRegionDB('regionDB/mm10',collections = "ccre")
repeatDB <- loadRegionDB('regionDB/mm10',collections = "repeats")
#repeatmaskerDB <- loadRegionDB('repeatmasker')
# plots
# ensembl - enrichment - genome-wide
res <- runLOLA(cons[[cluster]],uni,ensemblDB,cores = 6)
ensembl.ol.plot.enrich <- res %>% mutate(reg = sub('.bed', '', filename),
                        sig = as.character(signif.num(qValue)),
                        q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                        q = -log10(q),
                    info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                        oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(ensembl.ol.plot.enrich,file=(sprintf('plots/%s.%s.%s.%s.ensembl.ol.plot.enrich.rda',cell_line,control,test,mark)))
# ensembl - enrichment - igr 
res <- runLOLA(qSet.igr,uni_igr,ensemblDB,cores = 6)
ensembl.ol.plot.enrich <- res %>% mutate(reg = sub('.bed', '', filename),
                                         sig = as.character(signif.num(qValue)),
                                         q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                         q = -log10(q),
                                         info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                         oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(ensembl.ol.plot.enrich,file=(sprintf('plots/%s.%s.%s.%s.ensembl.ol.plot.enrich.igr.rda',cell_line,control,test,mark)))
# ensembl - enrichment - genic
res <- runLOLA(qSet.g,uni_g,ensemblDB,cores = 6)
ensembl.ol.plot.enrich <- res %>% mutate(reg = sub('.bed', '', filename),
                                         sig = as.character(signif.num(qValue)),
                                         q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                         q = -log10(q),
                                         info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                         oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(ensembl.ol.plot.enrich,file=(sprintf('plots/%s.%s.%s.%s.ensembl.ol.plot.enrich.g.rda',cell_line,control,test,mark)))
# ensembl - genomewide - depletion
res <- runLOLA(cons[[cluster]],uni,ensemblDB,cores = 6,direction = "depletion")
ensembl.ol.plot.dep <- res %>% mutate(reg = sub('.bed', '', filename),
                                  sig = as.character(signif.num(qValue)),
                                  q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                  q = -log10(q),
                                  info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                  oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(ensembl.ol.plot.dep,file=(sprintf('plots/%s.%s.%s.%s.ensembl.ol.plot.dep.rda',cell_line,control,test,mark)))
# ensembl - igr - depletion
res <- runLOLA(qSet.igr,uni_igr,ensemblDB,cores = 6,direction = "depletion")
ensembl.ol.plot.dep <- res %>% mutate(reg = sub('.bed', '', filename),
                                      sig = as.character(signif.num(qValue)),
                                      q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                      q = -log10(q),
                                      info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                      oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(ensembl.ol.plot.dep,file=(sprintf('plots/%s.%s.%s.%s.ensembl.ol.plot.dep.igr.rda',cell_line,control,test,mark)))
#ensembl - genic - depletion
res <- runLOLA(qSet.g,uni_g,ensemblDB,cores = 6,direction = "depletion")
ensembl.ol.plot.dep <- res %>% mutate(reg = sub('.bed', '', filename),
                                      sig = as.character(signif.num(qValue)),
                                      q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                      q = -log10(q),
                                      info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                      oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(ensembl.ol.plot.dep,file=(sprintf('plots/%s.%s.%s.%s.ensembl.ol.plot.dep.g.rda',cell_line,control,test,mark)))
# ccreDB enrichment - genomewide
res <- runLOLA(cons[[cluster]],uni,ccreDB,cores = 6)
ccre.ol.plot.enrich <- res %>% mutate(reg = sub('.bed', '', filename),
                              sig = as.character(signif.num(qValue)),
                              q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                              q = -log10(q),
                              info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                              oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(ccre.ol.plot.enrich,file=(sprintf('plots/%s.%s.%s.%s.ccre.ol.plot.enrich.rda',cell_line,control,test,mark)))
# igr - enrichment -ccreDB
res <- runLOLA(qSet.igr,uni_igr,ccreDB,cores = 6)
ccre.ol.plot.enrich <- res %>% mutate(reg = sub('.bed', '', filename),
                                      sig = as.character(signif.num(qValue)),
                                      q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                      q = -log10(q),
                                      info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                      oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(ccre.ol.plot.enrich,file=(sprintf('plots/%s.%s.%s.%s.ccre.ol.plot.enrich.igr.rda',cell_line,control,test,mark)))
# genic - enrichment - ccre
res <- runLOLA(qSet.g,uni_g,ccreDB,cores = 6)
ccre.ol.plot.enrich <- res %>% mutate(reg = sub('.bed', '', filename),
                                      sig = as.character(signif.num(qValue)),
                                      q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                      q = -log10(q),
                                      info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                      oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(ccre.ol.plot.enrich,file=(sprintf('plots/%s.%s.%s.%s.ccre.ol.plot.enrich.g.rda',cell_line,control,test,mark)))
# ccreDB depletion - genomewide
res <- runLOLA(cons[[cluster]],uni,ccreDB,cores = 6,direction = "depletion")
ccre.ol.plot.dep <- res %>% mutate(reg = sub('.bed', '', filename),
                                      sig = as.character(signif.num(qValue)),
                                      q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                      q = -log10(q),
                                      info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                      oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(ccre.ol.plot.dep,file=(sprintf('plots/%s.%s.%s.%s.ccre.ol.plot.dep.rda',cell_line,control,test,mark)))
# igr - ccre - depletion
res <- runLOLA(qSet.igr,uni_igr,ccreDB,cores = 6,direction = "depletion")
ccre.ol.plot.dep <- res %>% mutate(reg = sub('.bed', '', filename),
                                   sig = as.character(signif.num(qValue)),
                                   q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                   q = -log10(q),
                                   info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                   oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(ccre.ol.plot.dep,file=(sprintf('plots/%s.%s.%s.%s.ccre.ol.plot.dep.igr.rda',cell_line,control,test,mark)))
# genic - ccre - depletion
res <- runLOLA(qSet.g,uni_g,ccreDB,cores = 6,direction = "depletion")
ccre.ol.plot.dep <- res %>% mutate(reg = sub('.bed', '', filename),
                                   sig = as.character(signif.num(qValue)),
                                   q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                   q = -log10(q),
                                   info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                   oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(ccre.ol.plot.dep,file=(sprintf('plots/%s.%s.%s.%s.ccre.ol.plot.dep.g.rda',cell_line,control,test,mark)))
#
# repeatDB - genomewide - enrichment
res <- runLOLA(cons[[cluster]],uni,repeatDB,cores = 6)
repeat.ol.plot.enrich <- res %>% mutate(reg = sub('.bed', '', filename),
                                      sig = as.character(signif.num(qValue)),
                                      q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                      q = -log10(q),
                                      info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                      oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(repeat.ol.plot.enrich,file=(sprintf('plots/%s.%s.%s.%s.repeat.ol.plot.enrich.rda',cell_line,control,test,mark)))
# igr enrichment repeatDB
res <- runLOLA(qSet.igr,uni_igr,repeatDB,cores = 6)
repeat.ol.plot.enrich <- res %>% mutate(reg = sub('.bed', '', filename),
                                      sig = as.character(signif.num(qValue)),
                                      q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                      q = -log10(q),
                                      info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                      oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(repeat.ol.plot.enrich,file=(sprintf('plots/%s.%s.%s.%s.repeat.ol.plot.enrich.igr.rda',cell_line,control,test,mark)))
# genic - enrichment - repeatDB
res <- runLOLA(qSet.g,uni_g,repeatDB,cores = 6)
repeat.ol.plot.enrich <- res %>% mutate(reg = sub('.bed', '', filename),
                                        sig = as.character(signif.num(qValue)),
                                        q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                        q = -log10(q),
                                        info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                        oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(repeat.ol.plot.enrich,file=(sprintf('plots/%s.%s.%s.%s.repeat.ol.plot.enrich.g.rda',cell_line,control,test,mark)))
# repeatDB - genomewide - depletion
res <- runLOLA(cons[[cluster]],uni,repeatDB,cores = 6,direction = "depletion")
repeat.ol.plot.dep <- res %>% mutate(reg = sub('.bed', '', filename),
                                   sig = as.character(signif.num(qValue)),
                                   q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                   q = -log10(q),
                                   info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                   oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(repeat.ol.plot.dep,file=(sprintf('plots/%s.%s.%s.%s.repeat.ol.plot.dep.rda',cell_line,control,test,mark)))
# repeatDB - igr - depletion
res <- runLOLA(qSet.igr,uni_igr,repeatDB,cores = 6,direction = "depletion")
repeat.ol.plot.dep <- res %>% mutate(reg = sub('.bed', '', filename),
                                   sig = as.character(signif.num(qValue)),
                                   q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                   q = -log10(q),
                                   info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                   oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(repeat.ol.plot.dep,file=(sprintf('plots/%s.%s.%s.%s.repeat.ol.plot.dep.igr.rda',cell_line,control,test,mark)))
# repeat - genic - depletion
res <- runLOLA(qSet.g,uni_g,repeatDB,cores = 6,direction = "depletion")
repeat.ol.plot.dep <- res %>% mutate(reg = sub('.bed', '', filename),
                                     sig = as.character(signif.num(qValue)),
                                     q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                     q = -log10(q),
                                     info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                     oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))) %>%
  hchart('scatter', hcaes(x = oddsRatio, y = q, size = support,
                          group = reg), maxSize = '10%') %>%
  hc_xAxis( title = list(text = 'Odds ratio'), type = 'logarithmic') %>%
  hc_yAxis(title = list(text = '-log10(padj)')) %>%
  hc_tooltip(pointFormatter = JS("function() {
      return this.options.info;
  }")) %>%
  hc_chart(zoomType = 'xy')
save(repeat.ol.plot.dep,file=(sprintf('plots/%s.%s.%s.%s.repeat.ol.plot.dep.g.rda',cell_line,control,test,mark)))