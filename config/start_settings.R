### Configuration
################################################################################
param=list()

# Set paths
source("config/tool_path_configuration.R")
param$path_to_git=path_to_scanalysis
if (is.null(path_to_basic_settings)) {param$path_to_basic_settings=file.path(param$path_to_git,"basic_settings.R")}
if (is.null(path_to_advanced_settings)) {param$path_to_basic_settings=file.path(param$path_to_git,"advanced_settings.R")}
if (getwd()!=param$path_to_git) {
  setwd(param$path_to_git)
}

# Set environment
renv::use_python(type = "virtualenv", name = file.path(param$path_to_git,"env/basic/virtualenvs/r-reticulate"))
renv::load(file.path(param$path_to_git,"env/basic"))

# Set R configuration
source(file.path(param$path_to_git,'config/configuration.R'))
