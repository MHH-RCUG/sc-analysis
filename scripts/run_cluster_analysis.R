### Configuration
################################################################################
param=list()

# Set paths
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
setwd(param$path_to_git)
param$scriptname = "scripts/pre-processing/cluster_analysis.Rmd"

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
# Load existing object
param$data = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/output/Testdata/pre-processing/data/sc.rds"

# Reference
param$species="human"



########## Advanced settings ########## 

# Set advanced parameter
source(file.path(param$path_to_git,"/advanced_settings.R"))




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
