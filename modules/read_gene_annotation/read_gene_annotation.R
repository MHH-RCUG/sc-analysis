Sys.setenv(RENV_ACTIVATE_PROJECT = FALSE)
#source("/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/modules/download_references/renv/activate.R")
source("./modules/read_gene_annotation/renv/activate.R")

## Read gene annotation
# We read gene annotation from file. 
# We generate several dictionaries to translate between Ensembl IDs, gene symbols, Entrez Ids, and Seurat gene names. 

### Parameter
################################################################################
param=list()

# Reference
#param$mart_dataset="hsapiens_gene_ensembl"
#param$annot_version=98
param$mart_dataset="mmusculus_gene_ensembl"
param$annot_version=103

# Standard parameter
#param$path_to_git="/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis"
param$path_to_git="."
param$path_reference=file.path(param$path_to_git, "references", param$mart_dataset, param$annot_version)
param$reference=paste0(param$mart_dataset, ".v", param$annot_version, ".annot.txt")
param$file_annot = file.path(param$path_reference, param$reference)
param$file_cc_genes = file.path(param$path_reference, "cell_cycle_markers.xlsx")

param$annot_main=c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession")


# Git directory and files to source must be done first, then all helper functions can be sourced
git_files_to_source = c("functions_biomart.R")
git_files_to_source = file.path(param$path_to_git, "R", git_files_to_source)
file_exists = purrr::map_lgl(git_files_to_source, file.exists)
if (any(!file_exists)) stop(paste("The following files could not be found:", paste(git_files_to_source[!file_exists], collapse=", "), ". Please check the git directory at '", param$path_to_git, "'.!"))
invisible(purrr::map(git_files_to_source, source))


### read_ensembl_annotation
################################################################################

# Load from file
annot_ensembl = read.delim(param$file_annot)

# Double-check if we got all required annotation, in case annotation file was read
check_annot_main = all(param$annot_main %in% colnames(annot_ensembl))
if (!check_annot_main) {
  stop("The annotation table misses at least one of the following columns: ", paste(param$annot_main, collapse=", "))
}

# Translate IDs
IDs_out = suppressWarnings(TranslateIDs(annot_ensembl, param$annot_main)) 
ensembl_to_seurat_rowname = IDs_out[[1]]
seurat_rowname_to_ensembl = IDs_out[[2]]
seurat_rowname_to_entrez = IDs_out[[3]]
annot_ensembl = IDs_out[[4]]


### read_cc_genes
################################################################################
# Use biomart to translate human cell cycle genes to the species of interest

# Load from file
genes_s = openxlsx::read.xlsx(param$file_cc_genes, sheet=1)
genes_g2m = openxlsx::read.xlsx(param$file_cc_genes, sheet=2)

# Convert Ensembl ID to Seurat-compatible unique rowname
genes_s = data.frame(Human_gene_name=genes_s$Human_gene_name, Species_gene_name=unname(ensembl_to_seurat_rowname[genes_s$Species_ensembl_id]))
genes_g2m = data.frame(Human_gene_name=genes_g2m$Human_gene_name, Species_gene_name=unname(ensembl_to_seurat_rowname[genes_g2m$Species_ensembl_id]))

