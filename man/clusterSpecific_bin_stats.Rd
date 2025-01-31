% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterSpecific_bin_stats.R
\name{clusterSpecific_bin_stats}
\alias{clusterSpecific_bin_stats}
\title{Extract results per bin within each cluster from a DESeq analysis}
\usage{
clusterSpecific_bin_stats(
  out_dir,
  genome_assembly,
  treated_sample_bedfiles,
  wildtype_sample_bedfiles,
  treated_condition_label,
  wildtype_condition_label,
  annotated_clusters = NULL,
  output_filename,
  return_results_for_all_bins = FALSE
)
}
\arguments{
\item{out_dir}{a character string specifying the output directory for the sizeFactor-normalized BED files.}

\item{genome_assembly}{a character string specifying the genome assembly. Allowed values include "hg38" or "mm10".}

\item{treated_sample_bedfiles}{a vector specifying the BED file(s) for the treated sample(s) with raw counts.}

\item{wildtype_sample_bedfiles}{a vector specifying the BED file(s) for the wildtype sample(s) with raw counts.}

\item{treated_condition_label}{a character string specifying the condition for the treated sample(s).}

\item{wildtype_condition_label}{a character string specifying the condition for the wildtype sample(s).}

\item{annotated_clusters}{a character string specifying the R object containing the annotated clusters generated using annotate_clust(). If no clusters are provided, then the function proceeds without associating any bins to a cluster.}

\item{output_filename}{a character string specifying the file name for the resulting table to be saved on disk.}

\item{return_results_for_all_bins}{a logical indicating whether to also return a table with all results, including bins not associated with any cluster. Defaults to FALSE.}
}
\value{
a csv file
}
\description{
Performs differential bin analysis on raw counts using DESeq2. Log2 fold changes are shrunk using the "apeglm" method. Results include base means, log2 fold changes, standard errors, p-values, adjusted p-values, genomic coordinates, and cluster assignments for each bin.
}
\details{
Title
}
\examples{
\dontrun{
clusterSpecific_bin_stats(out_dir, mm10, treated_sample_bedfiles = c("NSD1KO_rep1.bed", "NSD1KO_rep2.bed"), wildtype_sample_bedfiles = c("WT_rep1.bed", "WT_rep2.bed"), treated_condition_label = "NSD1KO", wildtype_condition_label = "WT", annotated_clusters = "input_directory/annotated_clusters.rda", output_filename = "results_per_bin_per_cluster")
}
}
