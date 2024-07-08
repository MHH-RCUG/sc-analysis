### Configuration
################################################################################
param=list()

# Set paths
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
setwd(param$path_to_git)
param$scriptname = "scripts/pre-processing/qc.Rmd"

# Set environment
renv::load(file.path(param$path_to_git,"env/basic"))

# Set standard parameter
source(file.path(param$path_to_git,"config/set_standard_parameter.R"))



### Parameter
################################################################################

########## Basic settings ########## 

# Set project name
param$project_id = "Testdata"

# Select data  
# Either load existing object
param$data = NULL
#param$data = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/data/BGHFO/sc.rds"
# Or use test dataset
# One of the datasets included in ./modules/download_test_datasets/
# Does only work for 10x data automatically at the moment
param$download_test_datasets=NULL
#param_advset$download_test_datasets="download_10x_pbmc_small_split2samples"
# Or set data path
#param$path_data = NULL
param$path_data = data.frame(name=c("sample1", "sample2"),  
                             type=c("10x"), 
                             path=c("/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/data/counts/sample1/",
                                    "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/data/counts/sample2/"))

# Reference
param$species="human"



########## Advanced settings ########## 

# Set advanced parameter
source(file.path(param$path_to_git,"/advanced_settings.R"))

# Filter for cells
#param_advset$cell_filter = list(nFeature_RNA=c(500, 2400), nCount_RNA=c(NA, NA), percent_mt=c(NA, 20))

# Normalization
#param_advset$norm="SCT"



### Run markdown
################################################################################

# Set output folder
param$path_out = file.path(param$path_to_git, "output", param$project_id, gsub(".Rmd", "", basename(param$scriptname)))

# Record advanced parameter settings
param = modifyList(x = param, val = param_advset)

# Create output directories
if (!file.exists(param$path_out)) dir.create(param$path_out, recursive=TRUE, showWarnings=FALSE)
# Run markdown
rmarkdown::render(file.path(param$path_to_git,param$scriptname), param = param, output_file = file.path(param$path_out, paste0(param$project_id,"_", gsub("Rmd", "html", basename(param$scriptname)))))
