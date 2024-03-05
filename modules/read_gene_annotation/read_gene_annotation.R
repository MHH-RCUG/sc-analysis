## Read gene annotation
# We read gene annotation from file. 
# We generate several dictionaries to translate between Ensembl IDs, gene symbols, Entrez Ids, and Seurat gene names. 

### Set reference
################################################################################
if (param$species=="human") {
  param$mart_dataset="hsapiens_gene_ensembl"
} else {
  if (param$species=="mouse") {
    param$mart_dataset="mmusculus_gene_ensembl"
  } else {
    param$mart_dataset=param$mart_dataset
  }
}

if (is.null(param$annot_version)) {
  param$annot_version=98
} else {
  param$annot_version=param$annot_version
}

param$path_reference=file.path(param$path_to_git, "references", param$mart_dataset, param$annot_version)
param$reference=paste0(param$mart_dataset, ".v", param$annot_version, ".annot.txt")
param$file_annot = file.path(param$path_reference, param$reference)
param$file_cc_genes = file.path(param$path_reference, "cell_cycle_markers.xlsx")


### Download reference if not existing
################################################################################
if (!file.exists(param$file_annot) | !file.exists(param$file_cc_genes)) {
  source(file.path(param$path_to_git, "modules/download_references/download_references.R"))
}


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

