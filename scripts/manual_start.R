### Configuration
################################################################################
param=list()

# set paths
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
param$path_to_basic_settings=file.path(param$path_to_git,"basic_settings.R")
param$path_to_advanced_settings=file.path(param$path_to_git,"advanced_settings.R")
setwd(param$path_to_git)

# set environment (one of both)
# basic renv
renv::activate(file.path(param$path_to_git,"env/basic"))
# ccc renv
renv::activate(file.path(param$path_to_git,"env/ccc"))

# set configuration
source(file.path(param$path_to_git,'config/configuration.R'))
