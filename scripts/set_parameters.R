### Parameter
################################################################################
param=list()

# Reference
param$species="mouse"


### Advanced settings
#param$mart_dataset=NULL
#param$annot_version=NULL
#param$mart_dataset="hsapiens_gene_ensembl"
#param$annot_version=98
param$mart_dataset="mmusculus_gene_ensembl"
param$annot_version=100
param$annot_main=c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession")
param$mart_attributes=c(c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession"), 
                        c("chromosome_name", "start_position", "end_position", "percentage_gene_gc_content", "gene_biotype", "strand", "description"))
param$biomart_mirror=NULL


# Standard parameter
#param$path_to_git="/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis"
param$path_to_git="."
# Git directory and files to source must be done first, then all helper functions can be sourced
git_files_to_source = c("functions_biomart.R")
git_files_to_source = file.path(param$path_to_git, "R", git_files_to_source)
file_exists = purrr::map_lgl(git_files_to_source, file.exists)
if (any(!file_exists)) stop(paste("The following files could not be found:", paste(git_files_to_source[!file_exists], collapse=", "), ". Please check the git directory at '", param$path_to_git, "'.!"))
invisible(purrr::map(git_files_to_source, source))

setwd('/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis')

################################################################################

# Download reference if not existing
if (param$species=="human" && is.null(param$annot_version)) {
  param$mart_dataset=="hsapiens_gene_ensembl"
  param$annot_version=98
}
if (param$species=="mouse" && is.null(param$annot_version)) {
  param$mart_dataset=="mmusculus_gene_ensembl"
  param$annot_version=98
}

param$path_reference=file.path(param$path_to_git, "references", param$mart_dataset, param$annot_version)
param$reference=paste0(param$mart_dataset, ".v", param$annot_version, ".annot.txt")
param$file_annot = file.path(param$path_reference, param$reference)
param$file_cc_genes = file.path(param$path_reference, "cell_cycle_markers.xlsx")

if (!file.exists(param$file_annot)) {
  source('./modules/download_references/download_references.R')
}

if (file.exists(param$file_cc_genes)) {
  source('./modules/download_references/download_references.R')
}
  

# Read gene annotation 
source('./modules/read_gene_annotation/read_gene_annotation.R')

