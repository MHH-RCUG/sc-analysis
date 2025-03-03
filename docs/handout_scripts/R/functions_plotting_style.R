#' Out plotting style.
#'
#'@param title The plot title.
#'@param col A vector of colours to use.
#'@param fill A vector of fill colours to use.
#'@param legend_title The legend title.
#'@param legend_position The legend position.
#'@param xlab The title of the x-axis.
#'@param ylab The title of the y-axis.
#'@return None, add as theme.
AddStyle = function(title=NULL, col=NULL, fill=NULL, legend_title=NULL, legend_position=NULL, xlab=NULL, ylab=NULL) {
  list(
    theme_light() + theme(panel.border = element_blank()), 
    if (!is.null(title)) ggtitle(title), 
    if (length(col) > 0) scale_colour_manual(values=col),
    if (length(fill) > 0)  scale_fill_manual(values=fill),
    if (!is.null(legend_title)) {
      labs(color=legend_title, fill=legend_title)
    } else {
      theme(legend.title = element_blank()) 
    },
    if (!is.null(legend_position)) theme(legend.position=legend_position),
    if (!is.null(xlab)) xlab(xlab),
    if (!is.null(ylab)) ylab(ylab)
  )
}

#' Transform list of feature plots into combined plot with predefined plotting style.
#'
#'@param p A list of feature plots.
#'@param title The plot title.
#'@param ncol The number of plots per row.
#'@return Patchwork plot with theme.
AddStyleMultiFeaturePlot = function(p, title=NULL, ncol=length(p)) {
  for (i in list_indices(p):(length(p)-1) ) {
    p[[i]] = p[[i]] + AddStyle(xlab = " ", ylab = " ") + theme_classic() + NoLegend()
  }  
  p[[length(p)]] = p[[length(p)]] + AddStyle(xlab = " ", ylab = " ", legend_position = "right") + theme_classic()
  p = patchwork::wrap_plots(p, ncol = ncol) + plot_annotation(title = title)
  return(p)
}

#' Transform list of violin plots into combined plot with predefined plotting style.
#'
#'@param p A list of violin plots.
#'@param title The plot title.
#'@param ncol The number of plots per row.
#'@return Patchwork plot with theme.
AddStyleMultiVln = function(p, title=NULL, ncol=length(p)) {
  for (i in list_indices(p):(length(p)-1) ) {
    p[[i]] = p[[i]] + AddStyle(xlab = "", ylab = "", legend_position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 0.8))
  }  
  p[[length(p)]] = p[[length(p)]] + AddStyle(xlab = "", ylab = "", legend_position = "right") + theme(axis.text.x = element_text(angle = 45, hjust = 0.8))
  p = patchwork::wrap_plots(p, ncol = ncol) + plot_annotation(title = title)
  return(p)
}



#' Generate a list of standard plots with predefined and enhanced plotting style.
#' Part of the handout scripts for additional visualization. 
#'
#'@param sc Seurat object. Output from sc_analysis workflow.
#'@param param Parameter list. Output from sc_analysis workflow.
#'@param ncol The number of plots per row.
#'@param plottheme ggplot theme, e.g. theme_bw(), theme_classic(), theme_light(), theme_void(), theme_minimal(). Default: theme_light() 
#'@param clustercolors Color of clusters.
#'@param samplecolors Color of samples.
#'@param celltypecolors Color of cell type annotation. 
#'@return List with plots.
ADStandardPlots = function(sc, param, ncol=NULL, plottheme=theme_light(), clustercolors=NULL, samplecolors=NULL, celltypecolors=NULL) {
 
  if (is.null(clustercolors)) clustercolors=param$col_clusters
  if (is.null(samplecolors)) samplecolors=param$col_samples
  if (is.null(celltypecolors)) celltypecolors=param$col_annotation
  
  if (length(levels(sc$orig.ident))<=4) {
    samples_ncol = ceiling(length(levels(sc$orig.ident)))
  } else if (length(levels(sc$orig.ident))<=8) {
    samples_ncol = ceiling(length(levels(sc$orig.ident))/2)
  } else if (length(levels(sc$orig.ident))<=12) {
    samples_ncol = ceiling(length(levels(sc$orig.ident))/3)
  } else {
    samples_ncol = ceiling(length(levels(sc$orig.ident))/4)
  }
  if (is.null(ncol)) ncol=samples_ncol
  
  
  p_list = list()
  p_list[["clusters"]] = Seurat::DimPlot(sc, reduction="umap", group.by="seurat_clusters", pt.size = param$pt_size) + 
    AddStyle(title="Clusters", xlab = "UMAP 1", ylab = "UMAP 2") +
    scale_color_manual(values=clustercolors) + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8))  +
    plottheme
  p_list[["clusters_separately"]] = Seurat::DimPlot(sc, reduction="umap", group.by="seurat_clusters", split.by = "orig.ident", pt.size = param$pt_size, ncol = ncol) + 
    AddStyle(title="Clusters per sample", xlab = "UMAP 1", ylab = "UMAP 2") +
    scale_color_manual(values=clustercolors) + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8)) +
    plottheme
  p_list[["sample"]] = Seurat::DimPlot(sc, reduction="umap", group.by="orig.ident", pt.size = param$pt_size) + 
    AddStyle(title="Samples", xlab = "UMAP 1", ylab = "UMAP 2") +
    scale_color_manual(values=samplecolors) + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8)) +
    plottheme
  p_list[["sample_separately"]] = Seurat::DimPlot(sc, reduction="umap", group.by="orig.ident", split.by = "orig.ident", pt.size = param$pt_size, ncol = ncol) + 
    AddStyle(title="Samples separately", xlab = "UMAP 1", ylab = "UMAP 2") +
    scale_color_manual(values=samplecolors) + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8)) +
    plottheme
  p_list[["cluster_annotation"]] = Seurat::DimPlot(sc, reduction="umap", group.by="annotation", pt.size = param$pt_size) + 
    AddStyle(title="Cluster annotation", xlab = "UMAP 1", ylab = "UMAP 2") +
    scale_color_manual(values=celltypecolors) +
    plottheme
  p_list[["cluster_annotation_separately"]] = Seurat::DimPlot(sc, reduction="umap", group.by="annotation", split.by = "orig.ident", pt.size = param$pt_size, ncol = ncol) + 
    AddStyle(title="Cluster annotation per sample", xlab = "UMAP 1", ylab = "UMAP 2") +
    scale_color_manual(values=celltypecolors) + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8)) +
    plottheme
  p_list[["cell_annotation"]] = Seurat::DimPlot(sc, reduction="umap", group.by="SingleR.labels", pt.size=param$pt_size) +
    AddStyle(title="Cell annotation", xlab = "UMAP 1", ylab = "UMAP 2") + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8)) +
    plottheme
  p_list[["cell_annotation_separately"]] = Seurat::DimPlot(sc, reduction="umap", group.by="SingleR.labels", split.by = "orig.ident", pt.size=param$pt_size, ncol = ncol) +
    AddStyle(title="Cell annotation per sample", xlab = "UMAP 1", ylab = "UMAP 2") + 
    guides(colour = guide_legend(override.aes = list(size=3), nrow = 8)) +
    plottheme
  
  qc_feature = c(paste0("nCount_", param$assay_raw), paste0("nFeature_", param$assay_raw), "percent_mt", "percent_ribo")
  for (n in seq(qc_feature)) {
    name = qc_feature[n]
    p_list[[name]] = suppressMessages(Seurat::FeaturePlot(sc, features=qc_feature[n]) + 
                                        AddStyle(title=name, xlab = "UMAP 1", ylab = "UMAP 2") + 
                                        scale_colour_gradient(low=param$col_bg, high=param$col)) +
                                        plottheme
    if (qc_feature[n]==paste0("nCount_", param$assay_raw) | qc_feature[n]==paste0("nFeature_", param$assay_raw)) {
      p_list[[name]] = suppressMessages(p_list[[name]] + scale_colour_gradient(low=param$col_bg, high=param$col, trans="log10"))
    }
  }
  
    return(p_list)
}



#' Generate Feature Plots with predefined and enhanced plotting style.
#' Part of the handout scripts for additional visualization. 
#'
#'@param sc Seurat object. Output from sc_analysis workflow.
#'@param param Parameter list. Output from sc_analysis workflow.
#'@param markers Vector with gene names.
#'@param ncol The number of plots per row.
#'@param plottheme ggplot theme, e.g. theme_bw(), theme_classic(), theme_light(), theme_void(), theme_minimal(). Default: theme_light() 
#'@param featurecolor Color (scaled) for features.
#'@param bgcolor Background color for cells without feature expression (lower end of scale)
#'@param labelcluster Define whether to show cluster labels. Default: TRUE.
#'@return Patchwork plot.
ADFeaturePlots = function(sc, param, markers, ncol=4, plottheme=theme_light(), featurecolor=NULL, bgcolor=NULL, labelcluster=TRUE) {
  
  if (is.null(featurecolor)) featurecolor=param$col
  if (is.null(bgcolor)) bgcolor=param$col_bg
  
  p_list = list()
  for (i in 1:length(markers)) { 
    p_list[[i]] = Seurat::FeaturePlot(sc, features=markers[i], pt.size = param$pt_size,
                                      cols=c(bgcolor, featurecolor),  
                                      combine=TRUE, label=labelcluster) + 
      AddStyle(title=markers[i], 
               xlab="", ylab="") + 
      plottheme + theme(legend.position = "bottom")
      
  }
  # Combine all plots
  p = patchwork::wrap_plots(p_list, ncol=ncol) + 
    patchwork::plot_annotation(title="UMAP, cells coloured by normalised gene expression data")
  return(p)
}



#' Generate Feature Plots with predefined and enhanced plotting style.
#' Part of the handout scripts for additional visualization. 
#'
#'@param sc Seurat object. Output from sc_analysis workflow.
#'@param param Parameter list. Output from sc_analysis workflow.
#'@param markers Vector with gene names.
#'@param group Define how to group the plot, e.g. by "clusters", "samples", or "celltypes". Default: "clusters".
#'@param ncol The number of plots per row.
#'@param plottheme ggplot theme, e.g. theme_bw(), theme_classic(), theme_light(), theme_void(), theme_minimal(). Default: theme_light() 
#'@param clustercolors Color of clusters.

#'@return Patchwork plot.
ADViolinPlots = function(sc, param, markers, group="clusters", ncol=4, plottheme=theme_light(), clustercolors=NULL, samplecolors=NULL, celltypecolors=NULL) {
  
  if (is.null(clustercolors)) clustercolors=param$col_clusters
  if (is.null(samplecolors)) samplecolors=param$col_samples
  if (is.null(celltypecolors)) celltypecolors=param$col_annotation
  
  if (group=="clusters") {
    group_ident = "seurat_clusters"
    group_color = clustercolors
  } else if (group=="samples") {
    group_ident = "orig.ident"
    group_color = samplecolors
  } else if (group=="celltypes") {
    group_ident = "annotation"
    group_color = celltypecolors
  }
  
  
  # Plot violin plots per marker gene, and combine it all at the end
  p_list = list()
  
  for(i in 1:length(markers)) { 
    p_list[[i]] = Seurat::VlnPlot(sc, features=markers[i], assay=param$norm, pt.size=0, cols=group_color, group.by = group_ident) + 
      AddStyle(title=markers[i], xlab="") + 
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
      plottheme
  }
  # Combine all plots
  p = patchwork::wrap_plots(p_list, ncol=ncol) + 
    patchwork::plot_annotation(title="Violin plot of for normalised gene expression data") & theme(legend.position="none")
  return(p)
}







