### Configuration
################################################################################
param=list()

# Location of rds object
param$data = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/output/Testdata/cluster_analysis/data/sc.rds"
# Output path (only needed for saving plots)
param$path_out = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/output/Testdata/additional_visualisation/"
# Location of inspect_rds.R script
# Needed packages: Seurat, ggplot2, magrittr
param$scriptpath = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/scripts/read_data/"


### Run inspect_rds.R
################################################################################
source(file.path(param$scriptpath,"inspect_rds.R"))


### Investigate plots
################################################################################
# Print plot
# All plots are saved in p_list[[]]. You can retrieve the respective plot by inserting the name, e.g. p_list[["clusters"]]
# Putting the cursor in the middle the brackets and pressing TAB will display all possibilities
p_list[[]]

# Save plot
# As above, you have to select a individual plot from the list, e.g. p_list[["clusters"]]
# If a file with the same name is existing in the folder already, it will be overwritten
# Change the file name as fitting, but do not change path or ".tiff"
# You can change the width and height as it will fit to your plot
tiff(filename = paste0(param$path_out, "plotname.tiff"), width = 2400, height = 1200, res = 300)
p_list[[]]
dev.off()
