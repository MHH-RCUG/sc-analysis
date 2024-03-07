# Setting up download_test_datasets renv 
# in /modules/download_test_datasets
renv::init(bioconductor = TRUE, bare = TRUE) # bioconductor version 3.16 for R 4.2

renv::install("remotes")
renv::install("tidyverse") # version 2.0.0; dbplyr version < "2.3.4" needed
remotes::install_version("dbplyr", version='2.2.1')
renv::install("openxlsx") # version 4.2.5.2 
renv::install("bioc::biomaRt") # version 2.0.0 
renv::install("Seurat") # version 5.0.2
renv::install("R.utils") # version 2.12.3

renv::install("kableExtra") # version 1.4.0
remotes::install_version("matrixStats", version="1.1.0")
renv::install("bioc::SingleR") # version 2.0.0
renv::install("ggsci") # version 3.0.1
renv::install("bioc::SingleCellExperiment") # version 1.20.1
renv::install("bioc::scran") # version 1.26.2
renv::install("viridis") # version 0.6.5
renv::install("pheatmap") # version 1.0.12
renv::install("sessioninfo") # version 1.2.2




