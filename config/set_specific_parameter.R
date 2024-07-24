### Specific parameter 
################################################################################

# Read advanced parameter
source(file.path(param$path_to_git,"/advanced_settings.R"))

# Set output folder
param$path_out = file.path(param$path_to_git, "output", param$project_id, gsub(".Rmd", "", basename(param$scriptname)))

# Set advanced parameter
param = modifyList(x = param, val = param_advset)

# Create output directories
if (!file.exists(param$path_out)) dir.create(param$path_out, recursive=TRUE, showWarnings=FALSE)