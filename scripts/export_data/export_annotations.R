### Export annotation
################################################################################

# Export annotation as excel file
openxlsx::write.xlsx(x=data.frame(seurat_id=rownames(sc), ensembl_gene_id=seurat_rowname_to_ensembl[rownames(sc)], row.names=rownames(sc)) %>%
                       dplyr::inner_join(annot_ensembl, by="ensembl_gene_id"),
                     file=file.path(param$path_out, "data", "gene_annotation.xlsx"), 
                     rowNames=FALSE, colNames=TRUE)
