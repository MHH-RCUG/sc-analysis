# Set renv
renv::init(bioconductor = TRUE, bare = TRUE) # bioconductor version 3.16 for R 4.2
# Project '/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/ccc' [renv 1.0.5]
renv::settings$snapshot.type("all")

# basic
renv::install("remotes")
renv::install("tidyverse")
renv::install("patchwork") # version 1.2.0
renv::install("ggsci") # version 3.0.2
renv::install("gridGraphics")

# ccc_analysis
renv::install("bioc::scater") # version 1.26.1
renv::install("bioc::ComplexHeatmap") # version 2.14.0
renv::install("bioc::dir.expiry") # version 1.6.0
renv::install("bioc::basilisk.utils") # version 1.10.0
renv::install("saezlab/liana")
remotes::install_version("Seurat", version="5.0.2")

# rmarkdown
renv::install("sessioninfo") # version 1.2.2
renv::install("markdown") # version 1.12
remotes::install_version("xfun", version="0.41")
renv::install("knitcitations") # version 1.0.12
renv::install("kableExtra") # version 1.4.0

# Further packages:
# scater 1.26.1
# scran 1.26.2
# knitr 1.45
# rmarkdown 2.26 !!!
# htmltools 0.5.7
# httpuv 1.6.14
# httr 1.4.7
# bslib 0.6.1
# shiny 1.8.0
# SeuratObject 5.0.1
# xfun 0.42 !!!
# OmnipathR

# Set virtualenv
reticulate::virtualenv_create("/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/ccc/virtualenvs/r-reticulate")
renv::use_python(type = "virtualenv", name = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/ccc/virtualenvs/r-reticulate")
# Python 3.8.10 in virtualenv; ./virtualenvs/r-reticulate
# In Project '/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/ccc' [renv 1.0.5]

reticulate::py_install("pandas", envname = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/ccc/virtualenvs/r-reticulate")


