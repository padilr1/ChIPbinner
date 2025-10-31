# ChIPbinner

## Installation

Users can install the package using the two methods indicated below:
```{r, eval=FALSE}
# install the 'remotes' package first and run the installation code below:
remotes::install_github("padilr1/ChIPbinner")

# download a '.tar.gz' file of ChIPbinner and run the installation code below:
install.packages("ChIPbinner_0.99.0.tar.gz",repos=NULL,type="source")

# the package vignette can be viewed using:
browseVignettes("ChIPbinner")
# if the vignette is still not visible, you can try forcing the installation of the package vignette directly:
remotes::install_github("padilr1/ChIPbinner",build_vignettes=TRUE)
```

_The authors are preparing to submit this package to Bioconductor._

## Introduction

ChIPbinner is an open-source R package designed to facilitate genome-wide analysis of broad histone modifications. This tool addresses limitations in existing peak-calling software, which often struggles to accurately detect diffuse and broad signals in ChIP-Seq, CUT&RUN/TAG or related datasets. ChIPbinner instead divides the genome into uniform bins, offering an unbiased, reference-agnostic method to identify and explore differential enrichment across genomic regions.

ChIPbinner provides users with functionalities for signal normalization and combining replicates, alongside visualization tools for exploratory analysis, including genic/intergenic scatterplots, principal component analysis and correlation plots. Finally, it allows users to identify clusters of similarly-behaving bins and statistically assess their overlap with specific classes of annotated regions. 

ChIPbinner improves on previously published software by offering a clustering approach that is independent of differential binding and including additional features for downstream analysis.

Curated databases and complete published datasets for use with ChIPbinner can be found in this <a href="https://github.com/padilr1/ChIPbinner_database.git">GitHub repository</a>.

## Workflow

For a step-by-step tutorial on how to run ChIPbinner please see the <a href="https://padilr1.github.io/ChIPbinnerWorkflow/"> workflow documentation </a>
