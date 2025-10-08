#' Annotate clusters using the AddMetaData function to add annotations to your Seurat object based on marker genes.
#' Checks for existence of all marker genes. Does not work well with long list of marker genes!
#' 
#'@param sc Seurat object. Output from sc_analysis workflow.
#'@param marker_genes A list of marker genes: e.g. marker_genes = list("B cells" = c("CD19", "CD79A"), "Macrophages" = c("CD68", "CD14"))
AnnotateClustersMarkerBased = function(sc, marker_genes) {
  cell_type_annotations = vector("character", length = ncol(sc))
  for (i in 1:ncol(sc)) {
    cell_id = colnames(sc)[i]
    cluster_id = as.integer(sc$seurat_clusters[[i]])
    cluster_markers = markers %>% dplyr::filter(cluster == cluster_id) %>% dplyr::select(gene) %>% as.vector()
    for (cell_type in names(marker_genes)) {
      if (all(marker_genes[[cell_type]] %in% cluster_markers$gene)) {
        cell_type_annotations[i] = cell_type
        break
      }
    }
  }
  sc = AddMetaData(sc, metadata = cell_type_annotations, col.name = "CellType")
  return(sc)
}



#' Annotate clusters using the AddMetaData function to add annotations to your Seurat object based on module scores.
#'
#'@param sc Seurat object. Output from sc_analysis workflow.
#'@param marker_genes A list of marker genes: e.g. marker_genes = list("B cells" = c("CD19", "CD79A"), "Macrophages" = c("CD68", "CD14"))
AnnotateClustersScoreBased = function(sc, marker_genes) {
  cell_type_annotations = vector("character", length = ncol(sc))
  for (i in 1:ncol(sc)) {
    cell_id = colnames(sc)[i]
    scores = sapply(paste0("module_", names(marker_genes),"1"), function(cell_type) {
      sc@meta.data[cell_id, cell_type]
    })
    cell_type_annotations[i] = names(scores)[which.max(scores)]
  }
  sc = AddMetaData(sc, metadata = cell_type_annotations, col.name = "CellType_scored")
  return(sc)
}




