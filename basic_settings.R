### Basic settings
################################################################################
# Needs to be set



### Project specification
# Set project name
param$project_id = "Testdata"



### Dataset
# Either load existing object
param$data = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/output/Testdata/cluster_analysis/data/sc.rds"

# Or use test dataset included in ./modules/download_test_datasets/
# Does only work for 10x data automatically at the moment
param$download_test_datasets="download_10x_pbmc_1k_healthyDonor_v3Chemistry"

# Or set data path
param$path_data = NULL
  
# Reference
param$species="human"



### For ccc
# Set sender and receiver cell types
param$sender = NULL
param$receiver = NULL

### For dataset mapping
# Set reference object
param$refdata = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/output/Testdata_2/cluster_analysis/data/sc.rds"