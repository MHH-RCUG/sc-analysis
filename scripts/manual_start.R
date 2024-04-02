### Configuration
################################################################################
param=list()
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
setwd(param$path_to_git)


### Environment
################################################################################
# basic renv
renv::activate(file.path(param$path_to_git,"env/basic"))

# ccc renv
renv::activate(file.path(param$path_to_git,"env/ccc"))
