### Requirements
################################################################################
# Both datasets (reference and query) must be Seurat objects (rds files) with a @meta.data$annotation annotation column with a cell type annotation or at least a cluster number!



### Configuration
################################################################################
param=list()

# set paths
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
setwd(param$path_to_git)
param$scriptname = "modules/dataset_mapping/dataset_mapping_seurat.Rmd"

# set environment
renv::load(file.path(param$path_to_git,"env/basic"))



### Parameter
################################################################################

########## Basic settings ########## 

# Dataset
param$download_test_datasets="download_10x_pbmc_small_split2samples"
#param$download_test_datasets=NULL



### Run script
################################################################################
if (!is.null(param$download_test_datasets)) {
  # Download test dataset 
  param$path_test_dataset=paste0(param$path_to_git, "/modules/download_test_datasets/", param$download_test_datasets, ".R")
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