# Setting up download_test_datasets renv 
# in /modules/download_test_datasets
renv::init(bioconductor = TRUE, bare = TRUE) # bioconductor version 3.16 for R 4.2
renv::install("Seurat") # version 5.0.1
renv::install("openxlsx") # version 4.2.5.2
renv::install("R.utils") # version 2.12.3
