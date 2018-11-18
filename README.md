# sctransform
## R package for modeling single cell UMI expression data using regularized negative binomial regression

This packaged was developed by Christoph Hafemeister in [Rahul Satija's lab](https://satijalab.org/) at the New York Genome Center. A previous version of this work was used in the paper [Developmental diversification of cortical inhibitory interneurons, Nature 555, 2018](https://github.com/ChristophH/in-lineage). We are currently working on integrating the functionality of this package into [Seurat](https://satijalab.org/seurat/), an R package designed for QC, analysis, and exploration of single cell RNA-seq data.

This package is in beta status, please sanity check any results, and kindly notify me of any issues you find.

## Quick start
`devtools::install_github(repo = 'ChristophH/sctransform')`  
`normalized_data <- sctransform::vst(umi_count_matrix)$y`
