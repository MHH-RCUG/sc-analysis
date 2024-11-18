### Configuration
################################################################################
param=list()

# set paths
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
param$path_to_basic_settings=file.path(param$path_to_git,"basic_settings.R")
param$path_to_advanced_settings=file.path(param$path_to_git,"advanced_settings.R")
setwd(param$path_to_git)
