# Setting up download_test_datasets renv 
# in /modules/download_test_datasets
renv::init(bioconductor = TRUE, bare = TRUE) # bioconductor version 3.16 for R 4.2

renv::install("remotes")
renv::install("tidyverse") # version 2.0.0; dbplyr version < "2.3.4" needed
remotes::install_version("dbplyr", version='2.2.1')
renv::install("openxlsx") # version 4.2.5.2 
renv::install("bioc::biomaRt") # version 2.54.1 
renv::install("Seurat") # version 5.0.2
