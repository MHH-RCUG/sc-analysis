### Configuration
################################################################################
param=list()
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
setwd(param$path_to_git)

renv::load(file.path(param$path_to_git,"env/basic"))



### Parameter
################################################################################

########## Basic settings ########## 

# set project name
param$project_name = "BG-HFO_vs_HFO"

# set reference object
param$ref = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/input_data/HFO/sc.rds"
# set query object
param$query = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/input_data/BGHFO/sc.rds"

# Pre-annotated cell types; column in reference dataset
param$celltype = "annotation"
param$reduction = 'umap'
param$predicted_score_threshold = 0.9
param$percent_predicted_cells_threshold = 0.1

# Set standard colors
param$colors = FALSE


########## Advanced settings ########## 

# Set dot size for umaps/tsne
param$pt_size = 0.5

# Set colors
param$col_palette_annotation = "ggsci::pal_igv"
param$col_annotation_ref = c('Cluster 1-CM' = "#99e600FF", 'Cluster 3-Mesenchyme' = "#99e6e6FF", 'Cluster 4-Liver' = "#e6b800FF",
                             'Cluster 6-Apoptotic' = "#000000FF", 'Cluster 2-AFE' = "#99BBFFFF", 'Cluster 7-ECs-Blood' = "#ffff1aFF", 
                             'Cluster 5-hPSCs' = "#CFCFCFFF")
param$col_annotation = c('EPP' = "#4682B4FF", 'VE' = "#E6CCFFFF", 'PE/ST' = "#EEC900FF", 'MES' = "#CFCFCFFF", 
                             'AE' = "#FFB3D9FF", 'MKs' = "#80FFE5FF", 'PFE' = "#80FFE5FF", 'CM' = "#EE2C2CFF",
                             'Ery' = "#B33C00FF", 'Mo/Mo'	= "#99BBFFFF", 'EC' = "#CD00CDFF", 'HPC' = "#FFCC99FF",
                             'ELC' = "#EE2C2CFF")
param$col = "navy"

# compute_combined_umap_projection parameters
param$pc_n = 13
param$umap_k = 30



### Standard parameter 
################################################################################

# set output folder
param$path_out = file.path(param$path_to_git,"output_data", param$project_name, "dataset_mapping")

# set rmarkdown parameter
param$author = Sys.info()[["user"]]
param$scriptname = "modules/dataset_mapping/dataset_mapping_seurat.Rmd"

# set standard colors
if (param$colors == TRUE) {
  print(TRUE)
  param$col_annotation = NULL
  param$col_annotation_ref = NULL
  param$col_palette_annotation = "ggsci::pal_igv"
}



### Run markdown
################################################################################
# Create output directories
if (!file.exists(param$path_out)) dir.create(param$path_out, recursive=TRUE, showWarnings=FALSE)
# Run markdown
rmarkdown::render(file.path(param$path_to_git,param$scriptname), param = param, output_file = file.path(param$path_out, paste0(param$project_name,"_", gsub("Rmd", "html", basename(param$scriptname)))))
