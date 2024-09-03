### Convert data
################################################################################

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
}



### Export as andata object (h5ad file)
# Convert Seurat v5 single cell object to anndata object
# Does not work for SCT at the moment
# However, it would contain counts and data layer
#scCustomize::as.anndata(x = sc, file_path = file.path(param$path_out, "data"), file_name = "sc_anndata.h5ad", assay="RNA", 
#                        main_layer = "data", other_layers = "counts", transer_dimreduc = TRUE, verbose = FALSE)

# scesy only worked for v3 and v4 objects
# Convert to V3/4/Assay structure
sc_v3 = scCustomize::Convert_Assay(seurat_object = sc, convert_to = "V3", assay = "RNA")

# Convert Seurat v3 single cell object to anndata object
adata = sceasy::convertFormat(sc_v3, from="seurat", to="anndata", outFile=NULL, assay=DefaultAssay(sc))
# Write to h5ad file
adata$write(file.path(param$path_out, "data", "sc_anndata.h5ad"), compression="gzip")



### Export count matrix
invisible(ExportSeuratAssayData(sc, 
                                dir=file.path(param$path_out, "data", "counts"), 
                                assays=param$assay_raw, 
                                slot="counts",
                                include_cell_metadata_cols=colnames(sc[[]]),
                                metadata_prefix=paste0(param$project_id, ".")))


### Export metadata
openxlsx::write.xlsx(x=sc[[]] %>% tibble::rownames_to_column(var="Barcode"), file=file.path(param$path_out, "data", "cell_metadata.xlsx"), rowNames=FALSE, colNames=TRUE)

