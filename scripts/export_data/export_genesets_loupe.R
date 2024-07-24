### Export gene sets for loupe
################################################################################

# Export gene sets (CC genes, known markers, per-cluster markers up- and downregulated, ...)
gene_lists = Misc(sc, "gene_lists")

# Remove empty gene sets
gene_lists = gene_lists[purrr::map_int(gene_lists, length) > 0]
loupe_genesets = purrr::map_dfr(names(gene_lists), function(n) {return(data.frame(List=n, Name=gene_lists[[n]]))})
loupe_genesets$Ensembl = seurat_rowname_to_ensembl[loupe_genesets$Name]
write.table(loupe_genesets, file=file.path(param$path_out, "data", "Loupe_genesets.csv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep=",")