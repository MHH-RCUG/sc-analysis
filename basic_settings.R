### Basic settings
################################################################################
# Needs to be set



### Project specification
# Set project name
param$project_id = "Testdata"



### Dataset
# Either load existing object
#param$data = NULL
param$data = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/output/Testdata/pre-processing/data/sc.rds"

# Or use test dataset included in ./modules/download_test_datasets/
# Does only work for 10x data automatically at the moment
param$download_test_datasets=NULL

# Or set data path
param$path_data = data.frame(name=c("Sample1", "Sample2"),  
                                 type=c("10x"), 
                                 path=c("/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/data/counts/sample1/",
                                        "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/data/counts/sample2/"))

# Reference
param$species="human"

