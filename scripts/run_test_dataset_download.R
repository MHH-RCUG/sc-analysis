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