Sys.setenv(RENV_ACTIVATE_PROJECT = FALSE)
#source("/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/modules/download_references/renv/activate.R")
source("./modules/download_references/renv/activate.R")

## Download references
# We download references from Ensembl and write the resulting table to file.
library(magrittr)

# Create output directories
if (!file.exists(param$path_reference)) dir.create(param$path_reference, recursive=TRUE, showWarnings=FALSE)


### download_ensembl_annotation
################################################################################

# Download from Ensembl and a tab separated txt will be created
if (!file.exists(param$file_annot)) {
  annot_mart = suppressWarnings(GetBiomaRt(biomart="ensembl", 
                                           dataset=param$mart_dataset, 
                                           mirror=param$biomart_mirror, 
                                           version=param$annot_version))
  annot_ensembl = biomaRt::getBM(mart=annot_mart, attributes=param$mart_attributes, useCache=FALSE)
  write.table(annot_ensembl, file=param$file_annot, sep="\t", col.names=TRUE, row.names=FALSE, append=FALSE)
  message("Gene annotation file was created at: ", param$file_annot)
  # Note: depending on the attributes, there might be more than one row per gene
}


### download_cc_genes
################################################################################
# Use biomart to download and translate human cell cycle genes to the species of interest and save them in a file

# Obtain from Ensembl
# Note: both mart objects must point to the same mirror for biomarT::getLDS to work
if (!file.exists(param$file_cc_genes)) {
  mart_human = suppressWarnings(GetBiomaRt(biomart="ensembl", 
                                           dataset="hsapiens_gene_ensembl", 
                                           mirror=param$biomart_mirror, 
                                           version=param$annot_version))
  mart_myspecies = suppressWarnings(GetBiomaRt(biomart="ensembl", 
                                           dataset=param$mart_dataset, 
                                           mirror=GetBiomaRtMirror(mart_human), 
                                           version=param$annot_version)) 
    
  # S phase marker
  genes_s = biomaRt::getLDS(attributes=c("ensembl_gene_id", "external_gene_name"), 
                            filters="external_gene_name", 
                            values=Seurat::cc.genes.updated.2019$s.genes, 
                            mart=mart_human, 
                            attributesL=c("ensembl_gene_id", "external_gene_name"), 
                            martL=mart_myspecies, 
                            uniqueRows=TRUE)
  colnames(genes_s) = c("Human_ensembl_id", "Human_gene_name", "Species_ensembl_id", "Species_gene_name")
  genes_s = genes_s %>% dplyr::arrange(Human_gene_name)
  
  # G2/M marker
  genes_g2m = biomaRt::getLDS(attributes=c("ensembl_gene_id", "external_gene_name"), 
                              filters="external_gene_name", 
                              values=Seurat::cc.genes.updated.2019$g2m.genes, 
                              mart=mart_human, 
                              attributesL=c("ensembl_gene_id", "external_gene_name"), 
                              martL=mart_myspecies, 
                              uniqueRows=TRUE)
  colnames(genes_g2m) = c("Human_ensembl_id", "Human_gene_name", "Species_ensembl_id", "Species_gene_name")
  genes_g2m = genes_g2m %>% dplyr::arrange(Human_gene_name)
  
  # Write to file
  openxlsx::write.xlsx(list(S_phase=genes_s,G2M_phase=genes_g2m), file=param$file_cc_genes)
}

