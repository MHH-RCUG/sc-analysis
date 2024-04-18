### Read rds object
################################################################################

### Load Seurat S4 objects 
# Test if file is defined
if (is.null(param$data)) {
  message("Dataset is not specified")
} else {
  # Test if file exists
  if (file.exists(param$data)) {
    # Read object
    message(paste0("Load dataset:", param$data))
    sc = base::readRDS(param$data)
    
    # Transfer original params to loaded object
    if ("parameters" %in% list_names(sc@misc[])) {
      orig_param = sc@misc$parameters
      
      param_keep = param[c("path_to_git", "scriptname", "author", "project_id", "data", "path_out", "file_annot", "file_cc_genes")]
      param = modifyList(x = param, val = orig_param)
      param = modifyList(x = param, val = param_keep)
      param = modifyList(x = param, val = param_advset)
    }
    
    ### Set colors
    # Set sample colors based on orig.ident
    if (!is.null(sc@meta.data[["orig.ident"]])) {
      if (length(unique(sc@meta.data[["orig.ident"]]))!=length(param$col_samples)) {
        message(paste0("No or wrong number of distinct colors for samples provieded. Generate colors using color palette", param$col_palette_samples))
        sample_colours_out = suppressWarnings(SetSampleColours(sc, param$col_palette_samples)) 
        sc = sample_colours_out[[1]]
        param$col_samples = sample_colours_out[[2]]
      }
    }
 
    # Set colors for clusters based on seurat_clusters 
    if (!is.null(sc@meta.data[["seurat_clusters"]])) {
      if (length(unique(sc@meta.data[["seurat_clusters"]]))!=length(param$col_clusters)) {
        message(paste0("No or wrong number of distinct colors for clusters provieded. Generate colors using color palette", param$col_palette_clusters))
        cluster_colours_out = suppressWarnings(SetClusterColours(sc, param$col_palette_clusters)) 
        sc = cluster_colours_out[[1]]
        param$col_clusters = cluster_colours_out[[2]]
      }
    }

    # Set colors for cell types based on annotation 
    if (!is.null(sc@meta.data[["annotation"]])) {
      if (length(unique(sc@meta.data[["annotation"]]))!=length(param$col_annotation)) {
        message(paste0("No or wrong number of distinct colors for cell tpye annotation provieded. Generate colors using color palette", param$col_palette_annotation))
        annotation_colours_out = suppressWarnings(SetAnnotationColours(sc, param$col_palette_annotation)) 
        sc = annotation_colours_out[[1]]
        param$col_annotation = annotation_colours_out[[2]]
      }
    }

  sc
  } else {
  message("Dataset does not exist")
  }
}


### Load reference Seurat S4 objects if specified
# Test if file is defined
if (!is.null(param$refdata)) {
  # Test if file exists
  if (file.exists(param$refdata)) {
   # Read object
    message(paste0("Load dataset:", param$refdata)) 
    scR = base::readRDS(param$refdata)
    
    # Transfer original params to loaded object
    if ("parameters" %in% list_names(scR@misc[])) {
      orig_paramR = scR@misc$parameters
      
      if (!is.null(orig_paramR$col_samples)) {
        param$col_samples_ref = orig_paramR$col_samples
      }
      if (!is.null(orig_paramR$col_clusters)) {
        param$col_clusters_ref = orig_paramR$col_clusters
      }
      if (!is.null(orig_paramR$col_annotation)) {
        param$col_annotation_ref = orig_paramR$col_annotation
      }
      param = modifyList(x = param, val = param_advset)
    }
    
    ### Set colors
    # Set sample colors based on orig.ident
    if (!is.null(scR@meta.data[["orig.ident"]])) {
      if (length(unique(scR@meta.data[["orig.ident"]]))!=length(param$col_samples_ref)) {
        message(paste0("No or wrong number of distinct colors for samples provieded. Generate colors using color palette", param$col_palette_samples))
        sample_colours_out = suppressWarnings(SetSampleColours(scR, param$col_palette_samples)) 
        scR = sample_colours_out[[1]]
        param$col_samples_ref = sample_colours_out[[2]]
      }
    }
    
    # Set colors for clusters based on seurat_clusters 
    if (!is.null(scR@meta.data[["seurat_clusters"]])) {
      if (length(unique(scR@meta.data[["seurat_clusters"]]))!=length(param$col_clusters_ref)) {
        message(paste0("No or wrong number of distinct colors for clusters provieded. Generate colors using color palette", param$col_palette_clusters))
        cluster_colours_out = suppressWarnings(SetClusterColours(scR, param$col_palette_clusters)) 
        scR = cluster_colours_out[[1]]
        param$col_clusters_ref = cluster_colours_out[[2]]
      }
    }
    
    # Set colors for cell types based on annotation 
    if (!is.null(scR@meta.data[["annotation"]])) {
      if (length(unique(scR@meta.data[["annotation"]]))!=length(param$col_annotation_ref)) {
        message(paste0("No or wrong number of distinct colors for cell tpye annotation provieded. Generate colors using color palette", param$col_palette_annotation))
        annotation_colours_out = suppressWarnings(SetAnnotationColours(scR, param$col_palette_annotation)) 
        scR = annotation_colours_out[[1]]
        param$col_annotation_ref = annotation_colours_out[[2]]
      }
    }
    
  scR
  } else {
  message("Reference dataset does not exist")
  }
}  
