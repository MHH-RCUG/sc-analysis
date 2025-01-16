### Configuration
################################################################################
#param=list()

# Set paths
#param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
#param$path_to_basic_settings=file.path(param$path_to_git,"basic_settings.R")
#param$path_to_advanced_settings=file.path(param$path_to_git,"advanced_settings.R")
#setwd(param$path_to_git)

source("config/start_settings.R")

# Set environment
#renv::use_python(type = "virtualenv", name = file.path(param$path_to_git,"env/basic/virtualenvs/r-reticulate"))
#renv::load(file.path(param$path_to_git,"env/basic"))

# Set parameter
source(file.path(param$path_to_git,"config/set_parameter.R"))



### Run script
################################################################################
source(file.path(param$path_to_git,"scripts/download_test_datasets/test_dataset_download.R"))
