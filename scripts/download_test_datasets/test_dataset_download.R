### Test dataset download

if (!is.null(param$download_test_datasets)) {
  param[["test"]] = "Works"
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