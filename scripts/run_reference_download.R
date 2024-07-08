### Configuration
################################################################################
param=list()

# Set paths
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
setwd(param$path_to_git)

# Set environment
renv::load(file.path(param$path_to_git,"env/basic"))
#source(file.path(param$path_to_git,'/config/configuration.R'))



### Parameter
################################################################################

########## Basic settings ########## 

# Reference
param$species="mouse"



########## Advanced settings ########## 

# Set advanced parameter
source(file.path(param$path_to_git,"/advanced_settings.R"))



### Read gene annotation 
################################################################################
# Record advanced parameter settings
param = modifyList(x = param, val = param_advset)

# Does also download a reference if not existing
source(file.path(param$path_to_git,"scripts/read_data/read_gene_annotation.R"))

