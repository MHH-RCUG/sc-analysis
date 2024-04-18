### Configuration
################################################################################
param=list()

# Set paths
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
setwd(param$path_to_git)

# Set environment
renv::load(file.path(param$path_to_git,"env/basic"))
source(file.path(param$path_to_git,'/scripts/configuration.R'))



### Parameter
################################################################################

########## Basic settings ########## 

# Reference
param$species="mouse"



########## Advanced settings ########## 

# Reference
# Default is Ensembl release 98 which corresponds to 2020-A reference package of 10x Genomics Cell Ranger
# Ensembl release 110 which corresponds to 2024-A reference package of 10x Genomics Cell Ranger
### ATTENTION: Translation cc genes between human and mouse does not work; Error in getLDS() since version 105 (https://github.com/grimbough/biomaRt/issues/66)
### Means versions > 105 do not work for mouse
param$annot_version=98
param$annot_main=c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession")
param$mart_attributes=c(c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession"), 
                        c("chromosome_name", "start_position", "end_position", "percentage_gene_gc_content", "gene_biotype", "strand", "description"))
param$biomart_mirror=NULL



### Read gene annotation 
################################################################################
# Does also download a reference if not existing
source(file.path(param$path_to_git,"modules/read_data/read_gene_annotation.R"))

