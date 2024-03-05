### Parameter
################################################################################
param=list()

# Standard parameter
param$path_to_git='/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis'
renv::activate(file.path(param$path_to_git,"env/basic"))
setwd(param$path_to_git)

# Git directory and files to source must be done first, then all helper functions can be sourced
git_files_to_source = c("functions_biomart.R")
git_files_to_source = file.path(param$path_to_git, "R", git_files_to_source)
file_exists = purrr::map_lgl(git_files_to_source, file.exists)
if (any(!file_exists)) stop(paste("The following files could not be found:", paste(git_files_to_source[!file_exists], collapse=", "), ". Please check the git directory at '", param$path_to_git, "'.!"))
invisible(purrr::map(git_files_to_source, source))


##########

# Reference
param$species="human"


########## Advanced settings ########## 

# Reference
# Default is Ensembl release 98 which is the same as used for cellranger references
param$annot_version=99
param$annot_main=c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession")
param$mart_attributes=c(c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession"), 
                        c("chromosome_name", "start_position", "end_position", "percentage_gene_gc_content", "gene_biotype", "strand", "description"))
param$biomart_mirror=NULL






################################################################################

# Read gene annotation 
# Does also download a reference if not existing
source(file.path(param$path_to_git,"modules/read_gene_annotation/read_gene_annotation.R"))

