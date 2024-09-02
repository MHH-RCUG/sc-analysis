### Parameter
################################################################################
# Read standard parameter
source(file.path(param$path_to_git,"config/standard_parameter.R"))

# Read basic settings
source(param$path_to_basic_settings)

# Read advanced parameter
source(param$path_to_advanced_settings)

# Set output folder
param$path_out = file.path(param$path_to_git, "output", param$project_id, gsub(".Rmd", "", basename(param$scriptname)))

# Set advanced parameter
param = modifyList(x = param, val = param_advset)

# Create output directories
if (!file.exists(param$path_out)) dir.create(param$path_out, recursive=TRUE, showWarnings=FALSE)

