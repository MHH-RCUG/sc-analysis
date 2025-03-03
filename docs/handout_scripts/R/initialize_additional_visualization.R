### Initialize additional visualization
################################################################################

### Required libraries
library(Seurat) # main
library(ggplot2) # plots
library(magrittr) # %>% operator



### Parameter setting 
param=list()
# Location of rds object
param$data = DataPath
# Output path (only needed for saving plots)
param$path_out = OutputPath
# Location of inspect_rds.R script
# Needed packages: Seurat, ggplot2, magrittr
param$scriptpath =  ScriptPath

# Source additional functions
source(file.path(param$scriptpath, "R/functions_plotting_style.R"))



### Load Seurat S4 objects 
# Test if file is defined
if (is.null(param$data)) {
  message("Dataset is not specified")
} else {
  # Test if file exists
  if (file.exists(param$data)) {
    # Read object
    message(paste0("Load dataset:", param$data))
    sc = base::readRDS(param$data)

  } else {
    message("Dataset does not exist")
  }
}


# Transfer original params to loaded object
if ("parameters" %in% names(sc@misc)) {
  # Retrieve previous parameter settings
  orig_param = sc@misc$parameters
  
  # Keep some parameter settings from object and project defined
  basic_param_keep = param[c("data", "path_out")]
  
  # Integrate parameter
  param = modifyList(x = param, val = orig_param)
  param = modifyList(x = param, val = basic_param_keep)
} else {
  message("No predefined parameters")
}


