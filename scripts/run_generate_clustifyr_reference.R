### Requirements
################################################################################
# Download the expression matrix and the meta data
# https://cellbrowser.readthedocs.io/en/master/load.html
# <Cell Browser dataset ID> can be found under 'Info & Download' - 'Dataset Information' - 'Abstract'

# Download manually in Linux: 
# wget https://cells.ucsc.edu/<Cell Browser dataset ID>/exprMatrix.tsv.gz
# wget https://cells.ucsc.edu/<Cell Browser dataset ID>/meta.tsv

# Save under a ref_data_name in 'references' subfolder



### Configuration
################################################################################
param=list()

# Set paths
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
param$path_to_basic_settings=file.path(param$path_to_git,"basic_settings.R")
param$path_to_advanced_settings=file.path(param$path_to_git,"advanced_settings.R")
setwd(param$path_to_git)
param$scriptname = "scripts/dataset_mapping/dataset_mapping_seurat.Rmd"

# Set environment
renv::use_python(type = "virtualenv", name = file.path(param$path_to_git,"env/basic/virtualenvs/r-reticulate"))
renv::load(file.path(param$path_to_git,"env/basic"))

# Set parameter
source(file.path(param$path_to_git,"config/set_parameter.R"))



### Run script
################################################################################
source(file.path(param$path_to_git,"scripts/download_references/generate_clustifyr_reference.R"))

