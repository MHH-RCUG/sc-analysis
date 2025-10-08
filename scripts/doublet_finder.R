### DoubletFinder
################################################################################

# Load object including metadata and environment
source("scripts/run_inspect_rds.R")
param$path_out = file.path(param$path_to_git,"output", param$project_id, "doubletfinder_new")
# Create output directories
if (!file.exists(param$path_out)) dir.create(param$path_out, recursive=TRUE, showWarnings=FALSE)

# load library
library(DoubletFinder)




### For all using function
source(file.path(param$path_to_git, "scripts/doublet_finder_function.R"))
DefaultAssay(sc) = "RNA"

sc_subsets <- SplitObject(sc, split.by = "orig.ident")
sc_subsets <- lapply(sc_subsets, run_doubletfinder_custom)

doubletfinder_res <- data.frame(dplyr::bind_rows(sc_subsets)) # merge to a single dataframe
rownames(doubletfinder_res) <- doubletfinder_res$row_names # assign cell IDs to row names to ensure match
doubletfinder_res$row_names <- NULL
sc <- AddMetaData(sc, doubletfinder_res, col.name = "doublet_finder")

p = DimPlot(sc, reduction = "umap", group.by = "doublet_finder", pt.size = param$pt_size, cols = c("darkred", "cadetblue3")) + 
  AddStyle(title="Cells coloured by cluster identity", legend_position="right", legend_title="", xlab = "UMAP 1", ylab = "UMAP 2")
tiff(filename = file.path(param$path_out,"doublets.tiff"), width = 1600, height = 1200, res = 300)
p
dev.off()

p = Seurat::DimPlot(sc, reduction="umap", group.by="seurat_clusters", pt.size=param$pt_size, split.by = "doublet_finder", ncol = 2) +
  #scale_color_manual(values=param$col_clusters) +
  AddStyle(title="Cells coloured by cluster identity", legend_position="right", legend_title="", xlab = "UMAP 1", ylab = "UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size=3), nrow = 12))
tiff(filename = file.path(param$path_out,"doublets_per_cluster.tiff"), width = 4000, height = 1400, res = 300)
p
dev.off()


# Check how doublets singlets differ in QC measures per sample.
p1 = VlnPlot(sc, group.by = 'seurat_clusters', split.by = "doublet_finder", cols = c("darkred", "cadetblue3"),
        features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0, log = TRUE) + 
  theme(legend.position = '')
p2 = VlnPlot(sc, group.by = 'seurat_clusters', split.by = "doublet_finder", cols = c("darkred", "cadetblue3"),
             features = "percent_mt", pt.size = 0, log = FALSE) + 
  theme(legend.position = 'right')
p = p1 | p2 + patchwork::plot_layout(2,1)
tiff(filename = file.path(param$path_out,"Vlnplot_doublets.tiff"), width = 6000, height = 1200, res = 300)
p
dev.off()


### Save object
saveRDS(sc, file = file.path(param$path_out,"sc_doublets.rds"))




### Save object
# Subset
sc_subset <- subset(sc, doublet_finder == "Singlet")

# First save the idents you want to keep 
idents.sc = data.frame("orig.ident" = sc_subset$orig.ident)

# Then create a new Seurat object with the counts layer and idents of your object :
sc_subset = CreateSeuratObject(counts = LayerData(sc_subset, assay = "RNA", layer = "counts"), project = param$project_id, meta.data = idents.sc)

### Save subsetted object
saveRDS(sc_subset, file = file.path(param$path_out,"sc_subset_singlets.rds"))



Idents(sc_subset) = "seurat_clusters"
sc_subset = subset(x = sc_subset, idents = c(10, 16, 19), invert = TRUE)

# First save the idents you want to keep
idents.sc = data.frame("orig.ident" = sc_subset$orig.ident)

# Then create a new Seurat object with the counts layer and idents of your object :
sc_subset = CreateSeuratObject(counts = LayerData(sc_subset, assay = "RNA", layer = "counts"), project = param$project_id, meta.data = idents.sc)

### Save subsetted object
saveRDS(sc_subset, file = file.path(param$path_out,"sc_subset_clusters.rds"))

