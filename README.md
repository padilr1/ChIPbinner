# ChIPbinner

ChIPbinner is an open-source R package designed to facilitate genome-wide analysis of broad histone modifications. This tool addresses limitations in existing peak-calling software, which often struggles to accurately detect diffuse and broad signals in ChIP-Seq, CUT&RUN/TAG or related datasets. ChIPbinner instead divides the genome into uniform bins, offering an unbiased, reference-agnostic method to identify and explore differential enrichment across genomic regions.

ChIPbinner provides users with functionalities for signal normalization and combining replicates, alongside visualization tools for exploratory analysis, including genic/intergenic scatterplots, principal component analysis and correlation plots. Finally, it allows users to identify clusters of similarly-behaving bins and statistically assess their overlap with specific classes of annotated regions. 
