#!/usr/bin/env Rscript
#' Title
#' @title Enrichment & depletion analysis for clusters of bins.
#' @description Performs overlap enrichment & depletion analysis between bins found in a specific cluster and a class of annotated regions from a curated database. Uses Fisher's exact test to assess statistical significance of the overlap against a background of all bins found in both cell lines being compared.
#' @param genome_assembly a character string specifying the genome assembly. Allowed values include "hg38" or "mm10".
#' @param annotated_clusters an R object containing annotated clusters, generated from the annotate_clust() function.
#' @param query_cluster a character string specifying the cluster being to be tested. Normally, the clusters are label alphabetically from A to Z. For example, "B" to specify cluster B.
#' @param pooled_bed_file a BED file generated from 'pre_clus()' consisting of pooled genomic coordinates found in both samples being compared. Used for the background.
#' @param functional_db an R object consisting of a curated database of functional annotations.
#' @param region a string specifying the region to be tested. Must be one of "genome_wide", "genic" or "intergenic". Genome_wide indicates global testing of bins, including those found in both intergenic and genic regions. By specifying the region to be either genic or intergenic, the user can evaluate exclusively genic or intergenic bins overlapping a specific class of annotated regions. In these cases, the background is stratified to only genic or intergenic regions to avoid spurious associations to annotations confounded by their predominantly genic or intergenic localization.
#' @param cores an integer specifying the number of parallel cores to be used.
#' @param n_elements an integer specifying the number of elements to be plotted. This will not affect the output table.
#' @param cutoff_for_overlap an integer specifying the minimum number of overlap between the bins in the specified cluster and the class of annotated regions.
#' @param file_plot_name a character string specifying the filename of the output plot.
#' @param output_table_name a character string specifying the filename for the output table.
#' @param width_of_plot a numeric specifying the width of the plot.
#' @param height_of_plot a numeric specifying the height of the plot.
#' @param out_dir a character string specifying the output directory.
#'
#' @return a table and plot of the results of the enrichment & depletion analysis.
#' @export
#'
#' @include annotate_clust.R
#'
#' @example inst/examples/example_enrich_clust.R
enrich_clust <- function(genome_assembly,
                         annotated_clusters,
                         query_cluster,
                         pooled_bed_file,
                         functional_db,
                         region,
                         cores,
                         n_elements,
                         cutoff_for_overlap,
                         file_plot_name,
                         output_table_name,
                         width_of_plot,
                         height_of_plot,
                         out_dir) {
  suppressWarnings({
    # function to filter
    '%!in%' <- function(x,y)!('%in%'(x,y))
    # out dir
    out_dir <- paste0(out_dir)
    # gene and intergenic reference files
    if (genome_assembly == "hg38") {
      gene <- hg38_gene
      igr <- hg38_igr
    } else if (genome_assembly == "mm10") {
      gene <- mm10_gene
      igr <- mm10_igr
    }
    # loading function
    loadRData <- function(fileName){
      #loads an RData file, and returns it
      load(fileName)
      get(ls()[ls() != "fileName"])
    }
    # load annotated clusters
    cons <- loadRData(paste0(annotated_clusters))
    # the cluster being analyzed, for example cluster A, B or C...
    cluster <- as.character(paste0(query_cluster))
    # load pooled BED file from pre_clus(), comprising of genomic coordinates found in the two samples being compared.
    pooled_BED <- rtracklayer::import.bed(pooled_bed_file)
    # number of parallel cores to use
    cores <- as.integer(paste0(cores))
    # load signif annotation
    signif.num <- function(x) {
      symnum(x,
        corr = FALSE, na = FALSE, legend = FALSE,
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        symbols = c("****", "***", "**", "*", "ns")
      )
    }
    # collection of functional annotations
    ## database of functional annotations to use
    DB <- loadRData(paste0(functional_db))
    # background
    uni <- pooled_BED
    uni_igr <- uni[IRanges::overlapsAny(uni, igr) & !IRanges::overlapsAny(uni, gene)]
    uni_g <- uni[IRanges::overlapsAny(uni, gene) & !IRanges::overlapsAny(uni, igr)]
    # user set
    qSet <- cons[[cluster]]
    qSet_igr <- qSet[IRanges::overlapsAny(qSet, igr) & !IRanges::overlapsAny(qSet, gene)]
    qSet_g <- qSet[IRanges::overlapsAny(qSet, gene) & !IRanges::overlapsAny(qSet, igr)]
    # run LOLA
    ## genome-wide, genic or intergenic can be set for the region
    region <- paste0(region)
    ## take the top most significantly enriched or depleted elements
    n_elements <- as.integer(paste0(n_elements))
    ## cut off for overlaps
    cutoff_for_overlap <- as.integer(paste0(cutoff_for_overlap))
    ## plot attributes
    file_plot_name <- paste0(file_plot_name)
    width_of_plot <- as.integer(paste0(width_of_plot))
    height_of_plot <- as.integer(paste0(height_of_plot))
    output_table_name <- paste0(output_table_name)
    ## filter out
    "%!in%" <- function(x, y) !("%in%"(x, y))
    if (region == "genome_wide") {
      ## GENOME_WIDE ##
      ### enrichment
      enrichment <- LOLA::runLOLA(qSet, uni, DB, cores = cores) %>% as.data.frame()
      if (as.integer(nrow(enrichment)) == 0 | !("qValue" %in% names(enrichment))) {
        print("No overlap enrichment with specified database.")
        enrichment <- NULL
      } else {
        enrichment <- enrichment %>%
          dplyr::mutate(
            reg = sub(".bed.gz", "", filename),
            sig = as.character(signif.num(qValue)),
            q = dplyr::case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
            q = -log10(q),
            info = sprintf("Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d", oddsRatio, qValue, support, size),
            oddsRatio = dplyr::case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))
          ) %>%
          dplyr::mutate(overlap = .$support) %>%
          dplyr::select(c("oddsRatio", "reg", "qValue", "overlap", "sig")) %>%
          dplyr::mutate(type = "enrichment")
        if (nrow(enrichment %>% dplyr::filter(sig %!in% "ns")) == 0){
          enrichment <- NULL
        } else {
          enrichment <- enrichment
        }
      }
      ### depletion
      depletion <- LOLA::runLOLA(qSet, uni, DB, cores = cores, direction = "depletion")
      if (as.integer(nrow(depletion)) == 0 | !("qValue" %in% names(depletion))) {
        print("No depletion with specified database.")
        depletion <- NULL
      } else {
        depletion <- depletion %>%
          dplyr::mutate(
            reg = sub(".bed.gz", "", filename),
            sig = as.character(signif.num(qValue)),
            q = dplyr::case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
            q = -log10(q),
            info = sprintf("Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d", oddsRatio, qValue, support, size),
            oddsRatio = dplyr::case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))
          ) %>%
          dplyr::mutate(overlap = .$support) %>%
          dplyr::select(c("oddsRatio", "reg", "qValue", "overlap", "sig")) %>%
          dplyr::mutate(type = "depletion")
        if (nrow(depletion %>% dplyr::filter(sig %!in% "ns")) == 0){
          depletion <- NULL
        } else {
          depletion <- depletion
        }
      }
    } else if (region == "genic") {
      ## GENIC ##
      ### enrichment
      enrichment <- LOLA::runLOLA(qSet_g, uni_g, DB, cores = cores) %>% as.data.frame() %>%
        dplyr::filter(filename %!in% c("gene.bed.gz"))
      if (as.integer(nrow(enrichment)) == 0 | !("qValue" %in% names(enrichment))) {
        print("No overlap enrichment with specified database.")
        enrichment <- NULL
      } else {
        enrichment <- enrichment %>%
          dplyr::mutate(
            reg = sub(".bed.gz", "", filename),
            sig = as.character(signif.num(qValue)),
            q = dplyr::case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
            q = -log10(q),
            info = sprintf("Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d", oddsRatio, qValue, support, size),
            oddsRatio = dplyr::case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))
          ) %>%
          dplyr::mutate(overlap = .$support) %>%
          dplyr::select(c("oddsRatio", "reg", "qValue", "overlap", "sig")) %>%
          dplyr::mutate(type = "enrichment")
        if (nrow(enrichment %>% dplyr::filter(sig %!in% "ns")) == 0){
          enrichment <- NULL
        } else {
          enrichment <- enrichment
        }
      }
      ### depletion
      depletion <- LOLA::runLOLA(qSet_g, uni_g, DB, cores = cores, direction = "depletion") %>%
        dplyr::filter(filename %!in% c("gene.bed.gz"))
      if (as.integer(nrow(depletion)) == 0 | !("qValue" %in% names(depletion))) {
        print("No depletion with specified database.")
        depletion <- NULL
      } else {
        depletion <- depletion %>%
          dplyr::mutate(
            reg = sub(".bed.gz", "", filename),
            sig = as.character(signif.num(qValue)),
            q = dplyr::case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
            q = -log10(q),
            info = sprintf("Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d", oddsRatio, qValue, support, size),
            oddsRatio = dplyr::case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))
          ) %>%
          dplyr::mutate(overlap = .$support) %>%
          dplyr::select(c("oddsRatio", "reg", "qValue", "overlap", "sig")) %>%
          dplyr::mutate(type = "depletion")
        if (nrow(depletion %>% dplyr::filter(sig %!in% "ns")) == 0){
          depletion <- NULL
        } else {
          depletion <- depletion
        }
      }
    } else if (region == "intergenic") {
      ## INTERGENIC ##
      ### enrichment
      enrichment <- LOLA::runLOLA(qSet_igr, uni_igr, DB, cores = cores) %>% as.data.frame() %>%
        dplyr::filter(filename %!in% c("intergenic.bed.gz"))
      if (as.integer(nrow(enrichment)) == 0 | !("qValue" %in% names(enrichment))) {
        print("No overlap enrichment with specified database.")
        enrichment <- NULL
      } else {
        enrichment <- enrichment %>%
          dplyr::mutate(
            reg = sub(".bed.gz", "", filename),
            sig = as.character(signif.num(qValue)),
            q = dplyr::case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
            q = -log10(q),
            info = sprintf("Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d", oddsRatio, qValue, support, size),
            oddsRatio = dplyr::case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))
          ) %>%
          dplyr::mutate(overlap = .$support) %>%
          dplyr::select(c("oddsRatio", "reg", "qValue", "overlap", "sig")) %>%
          dplyr::mutate(type = "enrichment")
        if (nrow(enrichment %>% dplyr::filter(sig %!in% "ns")) == 0){
          enrichment <- NULL
        } else {
          enrichment <- enrichment
        }
      }
      ### depletion
      depletion <- LOLA::runLOLA(qSet_igr, uni_igr, DB, cores = cores, direction = "depletion") %>%
        dplyr::filter(filename %!in% c("intergenic.bed.gz"))
      if (as.integer(nrow(depletion)) == 0 | !("qValue" %in% names(depletion))) {
        print("No depletion with specified database.")
        depletion <- NULL
      } else {
        depletion <- depletion %>%
          dplyr::mutate(
            reg = sub(".bed.gz", "", filename),
            sig = as.character(signif.num(qValue)),
            q = dplyr::case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
            q = -log10(q),
            info = sprintf("Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d", oddsRatio, qValue, support, size),
            oddsRatio = dplyr::case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))
          ) %>%
          dplyr::mutate(overlap = .$support) %>%
          dplyr::select(c("oddsRatio", "reg", "qValue", "overlap", "sig")) %>%
          dplyr::mutate(type = "depletion")
        if (nrow(depletion %>% dplyr::filter(sig %!in% "ns")) == 0){
          depletion <- NULL
        } else {
          depletion <- depletion
        }
      }
    } else {
      print("Please use genome_wide, genic or intergenic region.")
      break
    }
    if (is.null(nrow(enrichment)) == FALSE & is.null(nrow(depletion)) == FALSE) {
      ### aggregate
      agg <- rbind(enrichment, depletion)
      ### plot only the top selected elements
      rbind(enrichment %>%
        dplyr::arrange(dplyr::desc(oddsRatio)) %>% dplyr::filter(overlap > cutoff_for_overlap) %>% dplyr::filter(sig %!in% "ns") %>% dplyr::slice(1:n_elements), depletion %>%
        dplyr::arrange(dplyr::desc(oddsRatio)) %>% dplyr::filter(overlap > cutoff_for_overlap) %>% dplyr::filter(sig %!in% "ns") %>% dplyr::slice(1:n_elements)) %>%
        dplyr::mutate(reg = forcats::fct_inorder(reg)) %>%
        ggplot2::ggplot(ggplot2::aes(x = oddsRatio, y = reg)) +
        ggplot2::geom_segment(ggplot2::aes(
          x = oddsRatio, y = reg, xend = oddsRatio, yend = reg,
          color = overlap
        )) +
        ggplot2::geom_point(ggplot2::aes(size = overlap), color = "#ff8c00", stat = "identity") +
        ggplot2::geom_text(ggplot2::aes(label = sig), vjust = -0.4, size = 4) +
        ggplot2::scale_size(name = "# overlaps") +
        ggplot2::labs(y = "Region", x = "Odds ratio") +
        ggplot2::scale_color_viridis_c(
          name = "# overlaps",
          begin = 0.2, end = 0.8,
          option = "A", guide = "legend",
          direction = 1
        ) +
        ggplot2::facet_wrap(. ~ type, scales = "free") +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill = "white"),
          axis.line.y = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_line(size = 0.5, color = "black"),
          axis.text = ggplot2::element_text(color = "black", size = 9),
          axis.text.x = ggplot2::element_text(color = "black", size = 9),
          axis.text.y = ggplot2::element_text(color = "black", size = 9),
          axis.title = ggplot2::element_text(color = "black", size = 9),
          panel.grid.major.y = ggplot2::element_line(color = "grey", linetype = "solid"),
          axis.ticks.y = ggplot2::element_line(color = "grey", linetype = "solid"),
          panel.grid.major.x = ggplot2::element_line(color = "grey", linetype = "dashed"),
          legend.key = ggplot2::element_rect(fill = "white"),
          legend.background = ggplot2::element_rect(fill = "white"),
          strip.text.x = ggplot2::element_text(color = "white", size = 9),
          strip.background.x = ggplot2::element_rect(fill = "black"),
          legend.position = "none"
        ) +
        ggplot2::scale_size_continuous(range = c(4, 8))
      # save plot
      ggplot2::ggsave(dpi = 600, filename = paste0(file_plot_name,".pdf"), path = out_dir, units = "in", width = width_of_plot, height = height_of_plot,device = cairo_pdf)
      # write csv of matrix
      readr::write_csv(agg, file = sprintf("%s/%s.csv", out_dir, output_table_name))
      print("Enrichment/depletion output successfully generated!")
    } else if (is.null(nrow(enrichment)) == FALSE & is.null(nrow(depletion)) == TRUE) {
      enrichment %>%
        dplyr::arrange(dplyr::desc(oddsRatio)) %>%
        dplyr::filter(overlap > cutoff_for_overlap) %>%
        dplyr::filter(sig %!in% "ns") %>%
        dplyr::slice(1:n_elements) %>%
        dplyr::mutate(reg = forcats::fct_inorder(reg)) %>%
        ggplot2::ggplot(ggplot2::aes(x = oddsRatio, y = reg)) +
        ggplot2::geom_segment(ggplot2::aes(
          x = oddsRatio, y = reg, xend = oddsRatio, yend = reg,
          color = overlap
        )) +
        ggplot2::geom_point(ggplot2::aes(size = overlap), color = "#ff8c00", stat = "identity") +
        ggplot2::geom_text(ggplot2::aes(label = sig), vjust = -0.4, size = 4) +
        ggplot2::scale_size(name = "# overlaps") +
        ggplot2::labs(y = "Region", x = "Odds ratio",title="Enrichment") +
        ggplot2::scale_color_viridis_c(
          name = "# overlaps",
          begin = 0.2, end = 0.8,
          option = "A", guide = "legend",
          direction = 1
        ) +
        # facet_wrap(. ~ type, scales = "free") +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill = "white"),
          plot.title = ggplot2::element_text(color="black",size=10,hjust=0.5),
          axis.line.y = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_line(size = 0.5, color = "black"),
          axis.text = ggplot2::element_text(color = "black", size = 9),
          axis.text.x = ggplot2::element_text(color = "black", size = 9,hjust=1),
          axis.text.y = ggplot2::element_text(color = "black", size = 9,hjust = 1),
          axis.title = ggplot2::element_text(color = "black", size = 9),
          panel.grid.major.y = ggplot2::element_line(color = "grey", linetype = "solid"),
          axis.ticks.y = ggplot2::element_line(color = "grey", linetype = "solid"),
          panel.grid.major.x = ggplot2::element_line(color = "grey", linetype = "dashed"),
          legend.key = ggplot2::element_rect(fill = "white"),
          legend.background = ggplot2::element_rect(fill = "white"),
          strip.text.x = ggplot2::element_text(color = "white", size = 9),
          strip.background.x = ggplot2::element_rect(fill = "black"),
          legend.position = "none"
        ) +
        ggplot2::scale_size_continuous(range = c(4, 8))
      # save plot
      ggplot2::ggsave(dpi = 600, filename = paste0(file_plot_name,".pdf"), path = out_dir, units = "in", width = width_of_plot, height = height_of_plot,device = cairo_pdf)
      # write csv of matrix
      readr::write_csv(enrichment, file = sprintf("%s/%s.csv", out_dir, output_table_name))
      print("Enrichment output successfully generated! No significant depletion found.")
    } else if (is.null(nrow(enrichment)) == TRUE & is.null(nrow(depletion)) == FALSE) {
      depletion %>%
        dplyr::arrange(dplyr::desc(oddsRatio)) %>%
        dplyr::filter(overlap > cutoff_for_overlap) %>%
        dplyr::filter(sig %!in% "ns") %>%
        dplyr::slice(1:n_elements) %>%
        dplyr::mutate(reg = forcats::fct_inorder(reg)) %>%
        ggplot2::ggplot(ggplot2::aes(x = oddsRatio, y = reg)) +
        ggplot2::geom_segment(ggplot2::aes(
          x = oddsRatio, y = reg, xend = oddsRatio, yend = reg,
          color = overlap
        )) +
        ggplot2::geom_point(ggplot2::aes(size = overlap), color = "#ff8c00", stat = "identity") +
        ggplot2::geom_text(ggplot2::aes(label = sig), vjust = -0.4, size = 4) +
        ggplot2::scale_size(name = "# overlaps") +
        ggplot2::labs(y = "Region", x = "Odds ratio",title = "Depletion") +
        ggplot2::scale_color_viridis_c(
          name = "# overlaps",
          begin = 0.2, end = 0.8,
          option = "A", guide = "legend",
          direction = 1
        ) +
        ggplot2::facet_wrap(. ~ type, scales = "free") +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill = "white"),
          plot.title = ggplot2::element_text(color="black",size=10,hjust=0.5),
          axis.line.y = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_line(size = 0.5, color = "black"),
          axis.text = ggplot2::element_text(color = "black", size = 9),
          axis.text.x = ggplot2::element_text(color = "black", size = 9,hjust=1),
          axis.text.y = ggplot2::element_text(color = "black", size = 9,hjust=1),
          axis.title = ggplot2::element_text(color = "black", size = 9),
          panel.grid.major.y = ggplot2::element_line(color = "grey", linetype = "solid"),
          axis.ticks.y = ggplot2::element_line(color = "grey", linetype = "solid"),
          panel.grid.major.x = ggplot2::element_line(color = "grey", linetype = "dashed"),
          legend.key = ggplot2::element_rect(fill = "white"),
          legend.background = ggplot2::element_rect(fill = "white"),
          strip.text.x = ggplot2::element_text(color = "white", size = 9),
          strip.background.x = ggplot2::element_rect(fill = "black"),
          legend.position = "none"
        ) +
        ggplot2::scale_size_continuous(range = c(4, 8))
      # save plot
      ggplot2::ggsave(dpi = 600, filename = paste0(file_plot_name,".pdf"), path = out_dir, units = "in", width = width_of_plot, height = height_of_plot,device = cairo_pdf)
      # write csv of matrix
      readr::write_csv(depletion, file = sprintf("%s/%s.csv", out_dir, output_table_name))
      print("Depletion output successfully generated! No significant enrichment found.")
    } else {
      print("No significant enrichment and depletion found with specified database!")
    }
  })
  print("Enrichment/depletion analysis finished!")
}
