### Configuration
################################################################################
#param=list()

# Set paths
#param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
#param$path_to_basic_settings=file.path(param$path_to_git,"basic_settings.R")
#param$path_to_advanced_settings=file.path(param$path_to_git,"advanced_settings.R")
#setwd(param$path_to_git)

source("config/start_settings.R")

# Set environment
#renv::load(file.path(param$path_to_git,"env/basic"))



### Run inspect_rds.R
################################################################################
source(file.path(param$path_to_git,"scripts/read_data/inspect_rds.R"))



### Investigate plots
################################################################################
# Print plot
# All plots are saved in p_list[[]]. You can retrieve the respective plot by inserting the name, e.g. p_list[["clusters"]]
# Putting the cursor in the middle the brackets and pressing TAB will display all possibilities

#p_list[[]]

# Save plot
# As above, you have to select a individual plot from the list, e.g. p_list[["clusters"]]
# If a file with the same name is existing in the folder already, it will be overwritten
# Change the file name as fitting, but do not change path or ".tiff"
# You can change the width and height as it will fit to your plot

#tiff(filename = paste0(param$path_out, "plotname.tiff"), width = 2400, height = 1200, res = 300)
#p_list[[]]
#dev.off()