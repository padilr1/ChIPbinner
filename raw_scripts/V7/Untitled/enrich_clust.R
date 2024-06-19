#!/usr/bin/env Rscript
#' Title
#' @title Enrichment & depletion analysis for clusters of bins.
#' @description Using LOLA, performs overlap enrichment & depletion analysis between bins found in a specific cluster and a class of annotated regions from a curated database. Uses Fisher's exact test to assess statistical significance of the overlap against a background of all bins found in both cell lines being compared.
#' @param gene Annotated genic regions in BED format. hg38 and mm10 are included in the package and can be accessed in 'extdata'.
#' @param intergenic Annotated intergenic regions in BED format. hg38 and mm10 are included in the package and can be accessed in 'extdata'.
#' @param assembly Must of be one of hg38 or mm10.
#' @param annotated_clusters The R object containing annotated clusters, generated from the annotate_clust() function.
#' @param query_cluster The specific cluster being tested - normally, the clusters are label alphabetically from A to Z.
#' @param pooled_bed_file The pooled BED file generated from 'pre_clus()' consisting of genomic coordinates found in both samples being compared. Used for the background.
#' @param functional_db Must be one of ensembl, ccre or repeat. These are the curated databases of functional annotations available for testing.
#' @param region Must be one of genome_wide, genic or intergenic. Genome_wide indicates global testing of bins, including those found in both intergenic and genic regions. By specifying the region to be either genic or intergenic, the user can evaluate exclusively genic or intergenic bins overlapping a specific class of annotated regions. In these cases, the background is stratified to only genic or intergenic regions to avoid spurious associations to annotations confounded by their predominantly genic or intergenic localization.
#' @param cores The number of parallel cores to be used.
#' @param n_elements The number of elements to be plotted. This will not affect the output table.
#' @param cutoff_for_overlap An integer indicating the minimum number of overlap between the bins in the specified cluster and the class of annotated regions.
#' @param file_plot_name The name/label for the output plot.
#' @param output_table_name The name/label for the output table.
#' @param width_of_plot The width of the plot.
#' @param height_of_plot The height of the plot.
#' @param out_dir Output directory.
#'
#' @return A table and plot of the results of the enrichment & depletion analysis.
#' @export
#'
#' @example inst/examples/example_enrich_clust.R
enrich_clust <- function(gene,
                         intergenic,
                         assembly,
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
    # out dir
    out_dir <- paste0(out_dir)
    # load reference files
    gene <- rtracklayer::import.bed(gene)
    igr <- rtracklayer::import.bed(intergenic)
    # load annotated clusters
    load(paste0(annotated_clusters))
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
    functional_db <- paste0(functional_db, "DB")
    if (assembly == "hg38") {
      if (functional_db == "ensemblDB") {
        ensemblDB <- LOLA::loadRegionDB(system.file("extdata/regionDB/hg38", package = "ChIPbinner"), collections = "ensembl")
        DB <- ensemblDB
      } else if (functional_db == "ccreDB") {
        ccreDB <- LOLA::loadRegionDB(system.file("extdata/regionDB/hg38", package = "ChIPbinner"), collections = "ccre")
        DB <- ccreDB
      } else if (functional_db == "repeatDB") {
        repeatDB <- LOLA::loadRegionDB(system.file("extdata/regionDB/hg38", package = "ChIPbinner"), collections = "repeats")
        DB <- repeatDB
      }
      print("Using hg38 assembly for functional annotations.")
    } else if (assembly == "mm10") {
      if (functional_db == "ensemblDB") {
        ensemblDB <- LOLA::loadRegionDB(system.file("extdata/regionDB/mm10", package = "ChIPbinner"), collections = "ensembl")
        DB <- ensemblDB
      } else if (functional_db == "ccreDB") {
        ccreDB <- LOLA::loadRegionDB(system.file("extdata/regionDB/mm10", package = "ChIPbinner"), collections = "ccre")
        DB <- ccreDB
      } else if (functional_db == "repeatDB") {
        repeatDB <- LOLA::loadRegionDB(system.file("extdata/regionDB/mm10", package = "ChIPbinner"), collections = "repeats")
        DB <- repeatDB
      }
      print("Using mm10 assembly for functional annotations.")
    } else {
      print("Please use hg38 or mm10 for genome assembly.")
      break
    }
    # background
    uni <- pooled_BED
    uni_igr <- uni[overlapsAny(uni, igr) & !overlapsAny(uni, gene)]
    uni_g <- uni[overlapsAny(uni, gene) & !overlapsAny(uni, igr)]
    # user set
    qSet <- cons[[cluster]]
    qSet_igr <- qSet[overlapsAny(qSet, igr) & !overlapsAny(qSet, gene)]
    qSet_g <- qSet[overlapsAny(qSet, gene) & !overlapsAny(qSet, igr)]
    # run LOLA
    ## genome-wide, genic or intergenic can be set for the region
    region <- paste0(region)
    ## take the top most significantly enriched or depleted elements
    n_elements <- as.integer(paste0(n_elements))
    ## cut off for overlaps
    cutoff_for_overlap <- as.integer(paste0(cutoff_for_overlap))
    ## plot attributes
    file_plot_name <- paste0(file_plot_name, ".png")
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
          mutate(
            reg = sub(".bed", "", filename),
            sig = as.character(signif.num(qValue)),
            q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
            q = -log10(q),
            info = sprintf("Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d", oddsRatio, qValue, support, size),
            oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))
          ) %>%
          mutate(overlap = .$support) %>%
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
          mutate(
            reg = sub(".bed", "", filename),
            sig = as.character(signif.num(qValue)),
            q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
            q = -log10(q),
            info = sprintf("Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d", oddsRatio, qValue, support, size),
            oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))
          ) %>%
          mutate(overlap = .$support) %>%
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
      enrichment <- LOLA::runLOLA(qSet_g, uni_g, DB, cores = cores) %>% as.data.frame()
      if (as.integer(nrow(enrichment)) == 0 | !("qValue" %in% names(enrichment))) {
        print("No overlap enrichment with specified database.")
        enrichment <- NULL
      } else {
        enrichment <- enrichment %>%
          mutate(
            reg = sub(".bed", "", filename),
            sig = as.character(signif.num(qValue)),
            q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
            q = -log10(q),
            info = sprintf("Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d", oddsRatio, qValue, support, size),
            oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))
          ) %>%
          mutate(overlap = .$support) %>%
          dplyr::select(c("oddsRatio", "reg", "qValue", "overlap", "sig")) %>%
          dplyr::mutate(type = "enrichment")
        if (nrow(enrichment %>% dplyr::filter(sig %!in% "ns")) == 0){
          enrichment <- NULL
        } else {
          enrichment <- enrichment
        }
      }
      ### depletion
      depletion <- LOLA::runLOLA(qSet_g, uni_g, DB, cores = cores, direction = "depletion")
      if (as.integer(nrow(depletion)) == 0 | !("qValue" %in% names(depletion))) {
        print("No depletion with specified database.")
        depletion <- NULL
      } else {
        depletion <- depletion %>%
          mutate(
            reg = sub(".bed", "", filename),
            sig = as.character(signif.num(qValue)),
            q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
            q = -log10(q),
            info = sprintf("Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d", oddsRatio, qValue, support, size),
            oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))
          ) %>%
          mutate(overlap = .$support) %>%
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
      enrichment <- LOLA::runLOLA(qSet_igr, uni_igr, DB, cores = cores) %>% as.data.frame()
      if (as.integer(nrow(enrichment)) == 0 | !("qValue" %in% names(enrichment))) {
        print("No overlap enrichment with specified database.")
        enrichment <- NULL
      } else {
        enrichment <- enrichment %>%
          mutate(
            reg = sub(".bed", "", filename),
            sig = as.character(signif.num(qValue)),
            q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
            q = -log10(q),
            info = sprintf("Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d", oddsRatio, qValue, support, size),
            oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))
          ) %>%
          mutate(overlap = .$support) %>%
          dplyr::select(c("oddsRatio", "reg", "qValue", "overlap", "sig")) %>%
          dplyr::mutate(type = "enrichment")
        if (nrow(enrichment %>% dplyr::filter(sig %!in% "ns")) == 0){
          enrichment <- NULL
        } else {
          enrichment <- enrichment
        }
      }
      ### depletion
      depletion <- LOLA::runLOLA(qSet_igr, uni_igr, DB, cores = cores, direction = "depletion")
      if (as.integer(nrow(depletion)) == 0 | !("qValue" %in% names(depletion))) {
        print("No depletion with specified database.")
        depletion <- NULL
      } else {
        depletion <- depletion %>%
          mutate(
            reg = sub(".bed", "", filename),
            sig = as.character(signif.num(qValue)),
            q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
            q = -log10(q),
            info = sprintf("Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d", oddsRatio, qValue, support, size),
            oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0]))
          ) %>%
          mutate(overlap = .$support) %>%
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
        arrange(desc(oddsRatio)) %>% dplyr::filter(overlap > cutoff_for_overlap) %>% dplyr::filter(sig %!in% "ns") %>% dplyr::slice(1:n_elements), depletion %>%
        arrange(desc(oddsRatio)) %>% dplyr::filter(overlap > cutoff_for_overlap) %>% dplyr::filter(sig %!in% "ns") %>% dplyr::slice(1:n_elements)) %>%
        mutate(reg = fct_inorder(reg)) %>%
        ggplot(aes(x = oddsRatio, y = reg)) +
        geom_segment(aes(
          x = oddsRatio, y = reg, xend = oddsRatio, yend = reg,
          color = overlap
        )) +
        geom_point(aes(size = overlap), color = "#ff8c00", stat = "identity") +
        geom_text(aes(label = sig), vjust = -0.4, size = 4) +
        scale_size(name = "# overlaps") +
        labs(y = "Region", x = "Odds ratio") +
        scale_color_viridis_c(
          name = "# overlaps",
          begin = 0.2, end = 0.8,
          option = "A", guide = "legend",
          direction = 1
        ) +
        facet_wrap(. ~ type, scales = "free") +
        coord_cartesian(clip = "off") +
        theme(
          panel.background = element_rect(fill = "white"),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.x = element_line(size = 0.5, color = "black"),
          axis.text = element_text(color = "black", size = 9),
          axis.text.x = element_text(color = "black", size = 9),
          axis.text.y = element_text(color = "black", size = 9),
          axis.title = element_text(color = "black", size = 9),
          panel.grid.major.y = element_line(color = "grey", linetype = "solid"),
          axis.ticks.y = element_line(color = "grey", linetype = "solid"),
          panel.grid.major.x = element_line(color = "grey", linetype = "dashed"),
          legend.key = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"),
          strip.text.x = element_text(color = "white", size = 9),
          strip.background.x = element_rect(fill = "black"),
          legend.position = "none"
        ) +
        scale_size_continuous(range = c(4, 8))
      # save plot
      ggsave(dpi = 600, filename = file_plot_name, path = out_dir, units = "in", width = width_of_plot, height = height_of_plot)
      # write csv of matrix
      readr::write_csv(agg, file = sprintf("%s/%s.csv", out_dir, output_table_name))
      print("Enrichment/depletion output successfully generated!")
    } else if (is.null(nrow(enrichment)) == FALSE & is.null(nrow(depletion)) == TRUE) {
      enrichment %>%
        arrange(desc(oddsRatio)) %>%
        dplyr::filter(overlap > cutoff_for_overlap) %>%
        dplyr::filter(sig %!in% "ns") %>%
        dplyr::slice(1:n_elements) %>%
        mutate(reg = fct_inorder(reg)) %>%
        ggplot(aes(x = oddsRatio, y = reg)) +
        geom_segment(aes(
          x = oddsRatio, y = reg, xend = oddsRatio, yend = reg,
          color = overlap
        )) +
        geom_point(aes(size = overlap), color = "#ff8c00", stat = "identity") +
        geom_text(aes(label = sig), vjust = -0.4, size = 4) +
        scale_size(name = "# overlaps") +
        labs(y = "Region", x = "Odds ratio",title="Enrichment") +
        scale_color_viridis_c(
          name = "# overlaps",
          begin = 0.2, end = 0.8,
          option = "A", guide = "legend",
          direction = 1
        ) +
        # facet_wrap(. ~ type, scales = "free") +
        coord_cartesian(clip = "off") +
        theme(
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(color="black",size=10,hjust=0.5),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.x = element_line(size = 0.5, color = "black"),
          axis.text = element_text(color = "black", size = 9),
          axis.text.x = element_text(color = "black", size = 9,hjust=1),
          axis.text.y = element_text(color = "black", size = 9,hjust = 1),
          axis.title = element_text(color = "black", size = 9),
          panel.grid.major.y = element_line(color = "grey", linetype = "solid"),
          axis.ticks.y = element_line(color = "grey", linetype = "solid"),
          panel.grid.major.x = element_line(color = "grey", linetype = "dashed"),
          legend.key = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"),
          strip.text.x = element_text(color = "white", size = 9),
          strip.background.x = element_rect(fill = "black"),
          legend.position = "none"
        ) +
        scale_size_continuous(range = c(4, 8))
      # save plot
      ggsave(dpi = 600, filename = file_plot_name, path = out_dir, units = "in", width = width_of_plot, height = height_of_plot)
      # write csv of matrix
      readr::write_csv(enrichment, file = sprintf("%s/%s.csv", out_dir, output_table_name))
      print("Enrichment output successfully generated! No significant depletion found.")
    } else if (is.null(nrow(enrichment)) == TRUE & is.null(nrow(depletion)) == FALSE) {
      depletion %>%
        arrange(desc(oddsRatio)) %>%
        dplyr::filter(overlap > cutoff_for_overlap) %>%
        dplyr::filter(sig %!in% "ns") %>%
        dplyr::slice(1:n_elements) %>%
        mutate(reg = fct_inorder(reg)) %>%
        ggplot(aes(x = oddsRatio, y = reg)) +
        geom_segment(aes(
          x = oddsRatio, y = reg, xend = oddsRatio, yend = reg,
          color = overlap
        )) +
        geom_point(aes(size = overlap), color = "#ff8c00", stat = "identity") +
        geom_text(aes(label = sig), vjust = -0.4, size = 4) +
        scale_size(name = "# overlaps") +
        labs(y = "Region", x = "Odds ratio",title = "Depletion") +
        scale_color_viridis_c(
          name = "# overlaps",
          begin = 0.2, end = 0.8,
          option = "A", guide = "legend",
          direction = 1
        ) +
        facet_wrap(. ~ type, scales = "free") +
        coord_cartesian(clip = "off") +
        theme(
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(color="black",size=10,hjust=0.5),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.x = element_line(size = 0.5, color = "black"),
          axis.text = element_text(color = "black", size = 9),
          axis.text.x = element_text(color = "black", size = 9,hjust=1),
          axis.text.y = element_text(color = "black", size = 9,hjust=1),
          axis.title = element_text(color = "black", size = 9),
          panel.grid.major.y = element_line(color = "grey", linetype = "solid"),
          axis.ticks.y = element_line(color = "grey", linetype = "solid"),
          panel.grid.major.x = element_line(color = "grey", linetype = "dashed"),
          legend.key = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"),
          strip.text.x = element_text(color = "white", size = 9),
          strip.background.x = element_rect(fill = "black"),
          legend.position = "none"
        ) +
        scale_size_continuous(range = c(4, 8))
      # save plot
      ggsave(dpi = 600, filename = file_plot_name, path = out_dir, units = "in", width = width_of_plot, height = height_of_plot)
      # write csv of matrix
      readr::write_csv(depletion, file = sprintf("%s/%s.csv", out_dir, output_table_name))
      print("Depletion output successfully generated! No significant enrichment found.")
    } else {
      print("No significant enrichment and depletion found with specified database!")
    }
  })
  print("Enrichment/depletion analysis finished!")
}
