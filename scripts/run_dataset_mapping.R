### Requirements
################################################################################
# Both datasets (reference and query) must be Seurat objects (rds files) with a @meta.data$annotation annotation column with a cell type annotation or at least a cluster number!



### Configuration
################################################################################
param=list()

# set paths
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
setwd(param$path_to_git)
param$scriptname = "modules/dataset_mapping/dataset_mapping_seurat.Rmd"

# set environment
renv::load(file.path(param$path_to_git,"env/basic"))

# set standard parameter
source(file.path(param$path_to_git,"scripts/set_standard_parameter.R"))



### Parameter
################################################################################

########## Basic settings ########## 

# set project name
param$project_name = "BG-HFO_vs_HFO"

# set reference object
param$ref = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/data/HFO/sc.rds"
# set query object
param$data = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/data/BGHFO/sc.rds"

# set output folder
param$path_out = file.path(param$path_to_git,"output", param$project_name, "dataset_mapping")

# Pre-annotated cell types; column in reference dataset
param$celltype = "annotation"   # at the moment required
param$reduction = 'umap'
param$predicted_score_threshold = 0.9
param$percent_predicted_cells_threshold = 0.1



########## Advanced settings ########## 

# Set colors
param$col_annotation_ref = c('Cluster 1-CM' = "#99e600FF", 'Cluster 3-Mesenchyme' = "#99e6e6FF", 'Cluster 4-Liver' = "#e6b800FF",
                             'Cluster 6-Apoptotic' = "#000000FF", 'Cluster 2-AFE' = "#99BBFFFF", 'Cluster 7-ECs-Blood' = "#ffff1aFF", 
                             'Cluster 5-hPSCs' = "#CFCFCFFF")
param$col_annotation = c('EPP' = "#4682B4FF", 'VE' = "#E6CCFFFF", 'PE/ST' = "#EEC900FF", 'MES' = "#CFCFCFFF", 
                             'AE' = "#FFB3D9FF", 'MKs' = "#80FFE5FF", 'PFE' = "#80FFE5FF", 'CM' = "#EE2C2CFF",
                             'Ery' = "#B33C00FF", 'Mo/Mo'	= "#99BBFFFF", 'EC' = "#CD00CDFF", 'HPC' = "#FFCC99FF",
                             'ELC' = "#EE2C2CFF")
param$col = "#0086b3"
param$col_bg = "#D3D3D3"

# compute_combined_umap_projection parameters
param$pc_n = 13
param$umap_k = 30



### Run markdown
################################################################################
# Create output directories
if (!file.exists(param$path_out)) dir.create(param$path_out, recursive=TRUE, showWarnings=FALSE)
# Run markdown
rmarkdown::render(file.path(param$path_to_git,param$scriptname), param = param, output_file = file.path(param$path_out, paste0(param$project_name,"_", gsub("Rmd", "html", basename(param$scriptname)))))
