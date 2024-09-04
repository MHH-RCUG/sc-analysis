# Set renv
renv::init(bioconductor = TRUE, bare = TRUE) # bioconductor version 3.16 for R 4.2
# Project '/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/basic' [renv 0.16.0]
renv::settings$snapshot.type("all")

# basic
renv::install("remotes")
renv::install("tidyverse") # version 2.0.0; dbplyr version < "2.3.4" needed
remotes::install_version("dbplyr", version='2.2.1')
renv::install("patchwork") # version 1.2.0
renv::install("openxlsx") # version 4.2.5.2 
renv::install("bioc::biomaRt") # version 2.0.0 
renv::install("Seurat") # version 5.0.2
renv::install("R.utils") # version 2.12.3

# pre-processing
renv::install("bioc::glmGamPoi") # version 1.10.2
renv::install("ape") # version 5.8
renv::install("clustertree") # version 0.5.1
renv::install("cellgeni/sceasy") # version 0.0.7
renv::install("ggsci") # version 3.0.1
renv::install("scCustomize") # version 2.1.2

# rmarkdown
renv::install("kableExtra") # version 1.4.0
renv::install("sessioninfo") # version 1.2.2
renv::install("markdown") # version 1.12
remotes::install_version("xfun", version="0.41")
renv::install("knitcitations") # version 1.0.12
renv::install("DT") # version 0.33

# DEG
renv::install("immunogenomics/presto") # version 1.0.0
renv::install("RGLab/MAST") # version 1.27.1
renv::install("enrichR") # version 3.2 

# dataset mapping
remotes::install_version("matrixStats", version="1.1.0")
renv::install("bioc::SingleR") # version 2.0.0
renv::install("bioc::celldex") # version 1.8.0
renv::install("bioc::SingleCellExperiment") # version 1.20.1
renv::install("bioc::scran") # version 1.26.2
renv::install("viridis") # version 0.6.5
renv::install("pheatmap") # version 1.0.12

# ccc analysis
renv::install("saezlab/liana") # version 0.1.14

# Further packages:
# reticulate 1.35.0
# knitr 1.45
# rmarkdown 2.25
# htmltools 0.5.7
# httpuv 1.6.14
# httr 1.4.7
# bslib 0.6.1
# shiny 1.8.0
# gridGraphics 0.5.1

# Set virtualenv
reticulate::virtualenv_create("/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/basic/virtualenvs/r-reticulate")
renv::use_python(type = "virtualenv", name = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/basic/virtualenvs/r-reticulate")
# Python 3.8.10 in virtualenv; ./virtualenvs/r-reticulate
# In Project '/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/basic' [renv 0.16.0]

#use_virtualenv("/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/basic/virtualenvs/r-reticulate")

reticulate::py_install("pandas", envname = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/basic/virtualenvs/r-reticulate")
reticulate::py_install("leidenalg", envname = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/basic/virtualenvs/r-reticulate")
reticulate::py_install("anndata", envname = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/basic/virtualenvs/r-reticulate")


