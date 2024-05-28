### Configuration
################################################################################
param=list()

# Set paths
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
setwd(param$path_to_git)
param$scriptname = "modules/pre-processing/qc.Rmd"

# Set environment
renv::load(file.path(param$path_to_git,"env/basic"))

# Set standard parameter
source(file.path(param$path_to_git,"scripts/set_standard_parameter.R"))



### Parameter
################################################################################

########## Basic settings ########## 

# Set project name
param$project_id = "Testdata"

# Select data:  
  # Either load existing object
  param$data = NULL
  #param$data = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/data/BGHFO/sc.rds"
  # Or set data path
  #param$path_data = NULL
  param$path_data = data.frame(name=c("sample1", "sample2"),  
                               type=c("10x"), 
                               path=c("/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/data/counts/sample1/",
                                      "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/data/counts/sample2/"))

# Set output folder
param$path_out = file.path(param$path_to_git, "output", param$project_id, "qc")

# Reference
param$species="human"



########## Advanced settings ########## 

### Dataset settings
# Dataset
# One of the datasets included in ./modules/download_test_datasets/
# Does only work for 10x data automatically at the moment
#param_advset$download_test_datasets="download_10x_pbmc_small_split2samples"
param$download_test_datasets=NULL

### Reference
# Default is Ensembl release 98 which corresponds to 2020-A reference package of 10x Genomics Cell Ranger
# Ensembl release 110 which corresponds to 2024-A reference package of 10x Genomics Cell Ranger
# ATTENTION: Translation cc genes between human and mouse does not work; Error in getLDS() since version 105 (https://github.com/grimbough/biomaRt/issues/66)
# Means versions > 105 do not work for mouse
param_advset$annot_version=NULL

### Normalisation
# Which normalisation should be used for analysis? ("RNA", "SCT")
param_advset$norm="SCT"

### Set colors
# Colour palette used for samples
param_advset$col_palette_samples = NULL
# Defined colours for samples
param_advset$col_samples = NULL
# Feature Plot colors - Highlights
param_advset$col = "#0086b3"
# Feature Plot colors - Background
param_advset$col_bg = "#D3D3D3"

### Record advanced parameter settings
param = modifyList(x = param, val = param_advset)


### Run markdown
################################################################################
# Create output directories
if (!file.exists(param$path_out)) dir.create(param$path_out, recursive=TRUE, showWarnings=FALSE)
# Run markdown
rmarkdown::render(file.path(param$path_to_git,param$scriptname), param = param, output_file = file.path(param$path_out, paste0(param$project_id,"_", gsub("Rmd", "html", basename(param$scriptname)))))
