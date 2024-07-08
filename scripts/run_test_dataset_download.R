### Configuration
################################################################################
param=list()

# set paths
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
setwd(param$path_to_git)

# set environment
renv::load(file.path(param$path_to_git,"env/basic"))
#source(file.path(param$path_to_git,'/config/configuration.R'))



### Parameter
################################################################################

########## Basic settings ########## 

# Dataset
# One of the datasets included in ./modules/download_test_datasets/
param$download_test_datasets="download_10x_pbmc_small_split2samples"



### Run script
################################################################################
if (!is.null(param$download_test_datasets)) {
  # Download test dataset 
  param$path_test_dataset=paste0(param$path_to_git, "/scripts/download_test_datasets/", param$download_test_datasets, ".R")
  if (file.exists(param$path_test_dataset)) {
    message(paste0("Using test dataset '", gsub('download_','', param$download_test_datasets), "'."))
    # Data output in a data subfolder of the directory where it is run 
    # Create output directories
    if (!file.exists("data")) dir.create("data", recursive=TRUE, showWarnings=FALSE)
    setwd(file.path(param$path_to_git,"data"))
    source(param$path_test_dataset)
  } else {
    message("Test dataset does not exist.")
  }
}