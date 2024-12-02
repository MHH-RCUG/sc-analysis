### Requirements
################################################################################
# Must be a Seurat object (rds file) with a @meta.data$annotation annotation column with a cell type annotation or at least a cluster number!
# Annotations must be known and some specified in param$sender and param$receiver



### Configuration
################################################################################
#param=list()

# Set paths
#param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
#param$path_to_basic_settings=file.path(param$path_to_git,"basic_settings.R")
#param$path_to_advanced_settings=file.path(param$path_to_git,"advanced_settings.R")
#setwd(param$path_to_git)



source("config/start_settings.R")

param$scriptname = "scripts/ccc_analysis/ccc_analysis.Rmd"

# set environment
#renv::use_python(type = "virtualenv", name = file.path(param$path_to_git,"env/basic/virtualenvs/r-reticulate"))
#renv::load(file.path(param$path_to_git,"env/basic"))

# Set parameter
source(file.path(param$path_to_git,"config/set_parameter.R"))



### Run markdown
################################################################################
rmarkdown::render(file.path(param$path_to_git,param$scriptname), param = param, output_file = file.path(param$path_out, paste0(param$project_id,"_", gsub("Rmd", "html", basename(param$scriptname)))))
