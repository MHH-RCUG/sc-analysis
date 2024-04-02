### Requirements
################################################################################
#


### Configuration
################################################################################
param=list()

# set paths
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
setwd(param$path_to_git)

# set environment
renv::load(file.path(param$path_to_git,"env/basic"))
source(file.path(param$path_to_git,'/scripts/configuration.R'))




### Parameter
################################################################################

########## Basic settings ########## 

# Reference
param$species="human"


########## Advanced settings ########## 

# Reference
# Default is Ensembl release 98 which is the same as used for cellranger references
param$annot_version=NULL
param$annot_main=c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession")
param$mart_attributes=c(c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession"), 
                        c("chromosome_name", "start_position", "end_position", "percentage_gene_gc_content", "gene_biotype", "strand", "description"))
param$biomart_mirror=NULL





### Read gene annotation 
################################################################################
# Does also download a reference if not existing
source(file.path(param$path_to_git,"modules/read_data/read_gene_annotation.R"))

