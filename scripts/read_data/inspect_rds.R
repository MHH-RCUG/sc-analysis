### Inspect rds object
################################################################################

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
  if ("colour_lists" %in% names(sc@misc)) {
    if ("orig.ident" %in% names(sc@misc$colour_lists)) {
      orig_param$col_samples = sc@misc$colour_lists$orig.ident
    }
    if ("seurat_clusters" %in% names(sc@misc$colour_lists)) {
      orig_param$col_clusters = sc@misc$colour_lists$seurat_clusters
    }
    if ("annotation" %in% names(sc@misc$colour_lists)) {
      orig_param$col_annotation = sc@misc$colour_lists$annotation
    }
  }
  
  # Keep some parameter settings from object and project defined
  basic_param_keep = param[c("data", "path_out")]
  
  # Integrate parameter
  param = modifyList(x = param, val = orig_param)
  param = modifyList(x = param, val = basic_param_keep)
} else {
  message("No predefined parameters")
}

# Set layout
if (length(levels(sc$orig.ident))<=4) {
  samples_ncol = ceiling(length(levels(sc$orig.ident)))
} else if (length(levels(sc$orig.ident))<=8) {
  samples_ncol = ceiling(length(levels(sc$orig.ident))/2)
} else if (length(levels(sc$orig.ident))<=12) {
  samples_ncol = ceiling(length(levels(sc$orig.ident))/3)
} else {
  samples_ncol = ceiling(length(levels(sc$orig.ident))/4)
}

  

### Generate multiple plots 
p_list = list()
if ("seurat_clusters" %in% colnames(sc@meta.data)) {
  p_list[["clusters"]] = Seurat::DimPlot(sc, reduction="umap", group.by="seurat_clusters", pt.size = param$pt_size) + 
    AddStyle(title="Clusters", xlab = "UMAP 1", ylab = "UMAP 2") +
    scale_color_manual(values=param$col_clusters) + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8))
  p_list[["clusters_separately"]] = Seurat::DimPlot(sc, reduction="umap", group.by="seurat_clusters", split.by = "orig.ident", pt.size = param$pt_size, ncol = samples_ncol) + 
    AddStyle(title="Clusters per sample", xlab = "UMAP 1", ylab = "UMAP 2") +
    scale_color_manual(values=param$col_clusters) + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8))
}
if ("orig.ident" %in% colnames(sc@meta.data)) {
  p_list[["sample"]] = Seurat::DimPlot(sc, reduction="umap", group.by="orig.ident", pt.size = param$pt_size) + 
    AddStyle(title="Samples", xlab = "UMAP 1", ylab = "UMAP 2") +
    scale_color_manual(values=param$col_samples) + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8))
  p_list[["sample_separately"]] = Seurat::DimPlot(sc, reduction="umap", group.by="orig.ident", split.by = "orig.ident", pt.size = param$pt_size, ncol = samples_ncol) + 
    AddStyle(title="Samples separately", xlab = "UMAP 1", ylab = "UMAP 2") +
    scale_color_manual(values=param$col_samples) + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8))
}
if ("annotation" %in% colnames(sc@meta.data)) {
  p_list[["cluster_annotation"]] = Seurat::DimPlot(sc, reduction="umap", group.by="annotation", pt.size = param$pt_size) + 
    AddStyle(title="Cluster annotation", xlab = "UMAP 1", ylab = "UMAP 2") +
    scale_color_manual(values=param$col_annotation) + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8)) 
  p_list[["cluster_annotation_separately"]] = Seurat::DimPlot(sc, reduction="umap", group.by="annotation", split.by = "orig.ident", pt.size = param$pt_size, ncol = samples_ncol) + 
    AddStyle(title="Cluster annotation per sample", xlab = "UMAP 1", ylab = "UMAP 2") +
    scale_color_manual(values=param$col_annotation) + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8)) 
}
if ("SingleR.labels" %in% colnames(sc@meta.data)) {
  p_list[["cell_annotation"]] = Seurat::DimPlot(sc, reduction="umap", group.by="SingleR.labels", pt.size=param$pt_size) +
    AddStyle(title="Cell annotation", xlab = "UMAP 1", ylab = "UMAP 2") + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8)) 
  p_list[["cell_annotation_separately"]] = Seurat::DimPlot(sc, reduction="umap", group.by="SingleR.labels", split.by = "orig.ident", pt.size=param$pt_size, ncol = samples_ncol) +
    AddStyle(title="Cell annotation per sample", xlab = "UMAP 1", ylab = "UMAP 2") + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8)) 
}  
 

qc_feature = c(paste0("nCount_", param$assay_raw), paste0("nFeature_", param$assay_raw))
if ("percent_mt" %in% colnames(sc@meta.data)) {
  qc_feature = c(qc_feature, "percent_mt")
}  
if ("percent_ribo" %in% colnames(sc@meta.data)) {
  qc_feature = c(qc_feature, "percent_ribo")
}  

for (n in seq(qc_feature)) {
  name = qc_feature[n]
  p_list[[name]] = suppressMessages(Seurat::FeaturePlot(sc, features=qc_feature[n]) + 
                                      AddStyle(title=name, xlab = "UMAP 1", ylab = "UMAP 2") + 
                                      scale_colour_gradient(low=param$col_bg, high=param$col)) 
  if (qc_feature[n]==paste0("nCount_", param$assay_raw) | qc_feature[n]==paste0("nFeature_", param$assay_raw)) {
    p_list[[name]] = suppressMessages(p_list[[name]] + scale_colour_gradient(low=param$col_bg, high=param$col, trans="log10"))
  }
}


