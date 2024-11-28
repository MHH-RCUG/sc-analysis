### Inspect rds object
################################################################################

### Create output directories
if (!file.exists(param$path_out)) dir.create(param$path_out, recursive=TRUE, showWarnings=FALSE)


### Required libraries
library(Seurat) # main
library(ggplot2) # plots
library(magrittr) # %>% operator


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

  } else {
    message("Dataset does not exist")
  }
}


# Transfer original params to loaded object
if ("parameters" %in% names(sc@misc)) {
  # Retrieve previous parameter settings
  orig_param = sc@misc$parameters
  
  # Keep some parameter settings from object and project defined
  basic_param_keep = param[c("data", "path_out")]
  
  # Integrate parameter
  param = modifyList(x = param, val = orig_param)
  param = modifyList(x = param, val = basic_param_keep)
} else {
  message("No predefined parameters")
}
  

### Generate multiple plots 
p_list = list()
p_list[["clusters"]] = Seurat::DimPlot(sc, reduction="umap", group.by="seurat_clusters", pt.size = param$pt_size) + 
  scale_color_manual(values=param$col_clusters) 
p_list[["clusters_separately"]] = Seurat::DimPlot(sc, reduction="umap", group.by="seurat_clusters", split.by = "orig.ident", pt.size = param$pt_size, ncol = 3) + 
  scale_color_manual(values=param$col_clusters) 
p_list[["sample"]] = Seurat::DimPlot(sc, reduction="umap", group.by="orig.ident", pt.size = param$pt_size) + 
  scale_color_manual(values=param$col_samples) 
p_list[["sample_separately"]] = Seurat::DimPlot(sc, reduction="umap", group.by="orig.ident", split.by = "orig.ident", pt.size = param$pt_size, ncol = 3) + 
  scale_color_manual(values=param$col_samples) 
p_list[["cluster_annotation"]] = Seurat::DimPlot(sc, reduction="umap", group.by="annotation", pt.size = param$pt_size) + 
  scale_color_manual(values=param$col_annotation) 
p_list[["cluster_annotation_separately"]] = Seurat::DimPlot(sc, reduction="umap", group.by="annotation", split.by = "orig.ident", pt.size = param$pt_size, ncol = 3) + 
  scale_color_manual(values=param$col_annotation) 
p_list[["cell_annotation"]] = Seurat::DimPlot(sc, reduction="umap", group.by="SingleR.labels", pt.size=param$pt_size) 
p_list[["cell_annotation_separately"]] = Seurat::DimPlot(sc, reduction="umap", group.by="SingleR.labels", split.by = "orig.ident", pt.size=param$pt_size, ncol = 3) 

qc_feature = c(paste0("nCount_", param$assay_raw), paste0("nFeature_", param$assay_raw), "percent_mt", "percent_ribo")
for (n in seq(qc_feature)) {
  name = qc_feature[n]
  p_list[[name]] = suppressMessages(Seurat::FeaturePlot(sc, features=qc_feature[n]) + 
                          scale_colour_gradient(low="lightgrey", high=param$col))
  if (qc_feature[n]==paste0("nCount_", param$assay_raw) | qc_feature[n]==paste0("nFeature_", param$assay_raw)) {
    p_list[[name]] = suppressMessages(p_list[[n]] + scale_colour_gradient(low="lightgrey", high=param$col, trans="log10"))
  }
}

### Save plot p_list as rds
# Can be loaded with p_list <- readRDS("p_list.rds")
saveRDS(p_list, file.path(param$path_out,"p_list.rds"))


