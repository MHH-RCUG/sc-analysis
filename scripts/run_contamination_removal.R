### Configuration
################################################################################
# Set configs
source("config/start_settings.R")

# Set parameter
source(file.path(param$path_to_git,"config/set_parameter.R"))



### Run script
################################################################################
for (i in param$cellranger_samples) {
  param$cellranger_dir = samples[i]
  source(file.path(param$path_to_git,"scripts/contamination_removal/soupx.R"))
}



