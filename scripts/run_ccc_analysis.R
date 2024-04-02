### Requirements
################################################################################
# Must be a Seurat object (rds file) with a @meta.data$annotation annotation column with a cell type annotation or at least a cluster number!
# Annotations must be known and some specified in param$sender and param$receiver



### Configuration
################################################################################
param=list()

# set paths
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
setwd(param$path_to_git)
param$scriptname = "modules/ccc_analysis/ccc_analysis.Rmd"

# set environment
renv::load(file.path(param$path_to_git,"env/ccc"))

# set standard parameter
source(file.path(param$path_to_git,"scripts/set_standard_parameter.R"))



### Parameter
################################################################################

########## Basic settings ########## 

# set project name
param$project_name = "BG-HFO"

# set query object
# at the moment, only one object
param$data = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/data/BGHFO/sc.rds"

# set output folder
param$path_out = file.path(param$path_to_git,"output", param$project_name, "ccc_analysis")

# Set sender and receiver cell types
param$sender = c("MES","CM", "PE/ST", "HPCs", "AE")
param$receiver = c("HPC","AE", "CM", "PE/ST")



########## Advanced settings ########## 

# Available methods: "connectome", "logfc", "natmi", "sca", "cellphonedb", "cytotalk", "call_squidpy", "call_cellchat", "call_connectome", "call_sca", "call_italk", "call_natmi"  
param$liana_methods = c("cellphonedb", "logfc", "sca", "natmi", "connectome")

param$liana_agg_rank_threshold = 0.01



### Run markdown
################################################################################
# Create output directories
if (!file.exists(param$path_out)) dir.create(param$path_out, recursive=TRUE, showWarnings=FALSE)
# Run markdown
rmarkdown::render(file.path(param$path_to_git,param$scriptname), param = param, output_file = file.path(param$path_out, paste0(param$project_name,"_", gsub("Rmd", "html", basename(param$scriptname)))))
