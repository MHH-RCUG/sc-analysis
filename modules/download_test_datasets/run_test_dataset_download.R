# Download test dataset 
param$path_test_dataset=paste0(param$path_to_git, "/modules/download_test_datasets/", param$download_test_datasets, ".R")
if (file.exists(param$path_test_dataset)) {
  message(paste0("Using test dataset '", gsub('download_','', param$download_test_datasets), "'."))
  source(param$path_test_dataset)
} else {
  message("Test dataset does not exist.")
}