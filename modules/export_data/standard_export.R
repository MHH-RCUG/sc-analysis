### Export data
################################################################################

### Add metadata to misc
# Add colour lists for orig.dataset
col = GenerateColours(num_colours=length(levels(sc$orig.dataset)), names=levels(sc$orig.dataset), palette=param$col_palette_samples, alphas=1)
sc = ScAddLists(sc, lists=list(orig.dataset=col), lists_slot="colour_lists")

# Add experiment details
Seurat::Misc(sc, "experiment") = list(project_id=param$project_id, date=Sys.Date(), species=gsub("_gene_ensembl", "", param$mart_dataset))

# Add parameter
Seurat::Misc(sc, "parameters") = param

# Add technical output (note: cannot use Misc function here)
sc@misc$technical = data.frame(ScrnaseqSessionInfo(param$path_to_git))




### Export to Loupe Cell Browser
if (all(param$path_data$type == "10x")) { 
  
  # Export reductions (umap, pca, others)
  for(r in Seurat::Reductions(sc)) {
    write.table(Seurat::Embeddings(sc, reduction=r)[,1:2] %>% as.data.frame() %>% tibble::rownames_to_column(var="Barcode"),
                file=file.path(param$path_out, "data", paste0("Loupe_projection_", r, ".csv")), col.names=TRUE, row.names=FALSE, quote=FALSE, sep=",")
  }
  
  # Export categorical metadata
  loupe_meta = sc@meta.data
  idx_keep = sapply(1:ncol(loupe_meta), function(x) !is.numeric(loupe_meta[,x]))
  write.table(x=loupe_meta[, idx_keep] %>% tibble::rownames_to_column(var="barcode"), 
              file=file.path(param$path_out, "data", "Loupe_metadata.csv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep=",")
  
  # Export gene sets (CC genes, known markers, per-cluster markers up- and downregulated, ...)
  gene_lists = Misc(sc, "gene_lists")
  
  # Remove empty gene sets
  gene_lists = gene_lists[purrr::map_int(gene_lists, length) > 0]
  loupe_genesets = purrr::map_dfr(names(gene_lists), function(n) {return(data.frame(List=n, Name=gene_lists[[n]]))})
  loupe_genesets$Ensembl = seurat_rowname_to_ensembl[loupe_genesets$Name]
  write.table(loupe_genesets, file=file.path(param$path_out, "data", "Loupe_genesets.csv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep=",")
}



### Export as andata object (h5ad file)
# Convert Seurat v5 single cell object to anndata object
scCustomize::as.anndata(x = sc, file_path = file.path(param$path_out, "data"), file_name = "sc_anndata.h5ad", )

# scesy only worked for v3 and v4 objects
# convert to v3
#sc_v3 = sc
#sc_v3[["RNA"]] = as(sc_v3[["RNA"]], "Assay")
# Convert Seurat v3 single cell object to anndata object
#adata = sceasy::convertFormat(sc_v3, from="seurat", to="anndata", outFile=NULL, assay=DefaultAssay(sc))
# Write to h5ad file
#adata$write(file.path(param$path_out, "data", "sc.h5ad"), compression="gzip")



### Export seurat object
saveRDS(sc, file=file.path(param$path_out, "data", "sc.rds"))



### Export count matrix
invisible(ExportSeuratAssayData(sc, 
                                dir=file.path(param$path_out, "data", "counts"), 
                                assays=param$assay_raw, 
                                slot="counts",
                                include_cell_metadata_cols=colnames(sc[[]]),
                                metadata_prefix=paste0(param$project_id, ".")))


### Export metadata
openxlsx::write.xlsx(x=sc[[]] %>% tibble::rownames_to_column(var="Barcode"), file=file.path(param$path_out, "data", "cell_metadata.xlsx"), rowNames=FALSE, colNames=TRUE)


### Export annotation as excel file
openxlsx::write.xlsx(x=data.frame(seurat_id=rownames(sc), ensembl_gene_id=seurat_rowname_to_ensembl[rownames(sc)], row.names=rownames(sc)) %>%
                       dplyr::inner_join(annot_ensembl, by="ensembl_gene_id"),
                     file=file.path(param$path_out, "data", "gene_annotation.xlsx"), 
                     rowNames=FALSE, colNames=TRUE)
