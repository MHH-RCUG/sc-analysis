### Annotate cells with modules
################################################################################

# Load object including metadata and environment
source("scripts/run_inspect_rds.R")
param$path_out = file.path(param$path_to_git,"output", param$project_id, "geneset_modules")
# Create output directories
if (!file.exists(param$path_out)) dir.create(param$path_out, recursive=TRUE, showWarnings=FALSE)


# Create a list of known marker genes for each cell type:
#marker_genes = list(
#  "B cells" = c("CD19", "CD79A", "MS4A1"),
#  "Macrophages" = c("CD68", "CD14", "CD163"),
  # Add more cell types and their marker genes as needed
#)

# OR 

# Read list
all_markers = readxl::read_xlsx("/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/output/BGHFO/celltype_marker_ETH.xlsx")
marker_genes = as.list(all_markers)


# Find markers of sc object
#markers <- FindAllMarkers(sc, 
#                          test.use = sc@misc$markers$test,  # You can change the test based on your preference
#                          only.pos = TRUE,      # If you want only positive markers
#                          min.pct = sc@misc$markers$min_pct,       # Minimum percentage of cells expressing the gene
#                          logfc.threshold = sc@misc$markers$log2FC)  # Minimum log fold change

# OR 

# Load previously generated results
markers = sc@misc$markers$results



### Annotate with module score

# Calculate module scores for feature expression programs in single cells
# https://satijalab.org/seurat/reference/addmodulescore

# Add module scores for each cell type: 
# Use the AddModuleScore function to calculate the module scores for each cell type.
# Adapted from sc <- AddModuleScore(sc, features = marker_genes, name = "module_")
for (cell_type in names(marker_genes)) {
  sc = AddModuleScore(sc, features = list(marker_genes[[cell_type]]), name = paste0("module_", cell_type))
}

# Annotate Seurat object with cell types based on module scores
sc = AnnotateClustersScoreBased(sc, marker_genes)

# Rename column name
namevector = colnames(sc[[]])
namevector = gsub("module_|1", "", namevector)
colnames(sc@meta.data) = namevector
# Rename CellType_scored column content
sc$CellType_scored = gsub("module_|1", "", sc$CellType_scored)

# For each cluster, calculate mean of module score of all each cells in cluster
#comb_score = data.frame(sc[[c("seurat_clusters", paste0("module_", names(marker_genes),"1"))]])
comb_score = data.frame(sc[[c("seurat_clusters", names(marker_genes))]])
comb_score_tbl = comb_score %>% dplyr::group_by(seurat_clusters) %>% dplyr::summarise(across(everything(), mean))
# Save table with module score 
write.csv(comb_score_tbl, file.path(param$path_out,"comb_score_tbl.csv"))

# Extract highest module score and name for each cluster
cluster_annotations = list()
for (i in levels(sc$seurat_clusters)) {
  cluster_score = comb_score_tbl %>% dplyr::filter(seurat_clusters == i) %>% as.data.frame() %>% dplyr::select(2:length(colnames(comb_score_tbl)))
  cluster_annotations[i] = names(cluster_score)[which.max(cluster_score)]
}


### Manually correct annotations
cluster_annotations$`1` = "late ECs"
cluster_annotations$`3` = "HSPCs"
cluster_annotations$`6` = "AE"
cluster_annotations$`8` = "AE"
cluster_annotations$`9` = "RUNX1+ EC"
#cluster_annotations$`14` = "undefined"



### Transfer annotations
# Retrieve seurat clusters classification for each cell
cluster_idents = sc$seurat_clusters
# Transform cluster_annotations into a vector
cluster_annotations_vector = as.vector(unlist(cluster_annotations[]))
# Exchange the cluster numbers with cluster_annotations in the cluster_idents factor
levels(cluster_idents) = cluster_annotations_vector
# store cluster_idents factor in annotation meta data slot
sc[["annotation"]] = cluster_idents
# Set Idents to seurat_clusters
#Idents(sc) = 'annotation'


### Visualize the annotated cell types
# Each cell in a DimPlot
p = DimPlot(sc, reduction = "umap", group.by = "CellType_scored")
tiff(filename = file.path(param$path_out,"CellType_scored.tiff"), width = 2400, height = 1200, res = 300)
p
dev.off()

# Feature plot for each module
p_list = list()
for (cell_type in names(marker_genes)) {
  #p_list[[]] = VlnPlot(sc, features = paste0("module_", cell_type,"1"), pt.size = 0)
  #p_list[[cell_type]] = FeaturePlot(sc, features = paste0("module_", cell_type,"1"), reduction = "umap", cols = c("lightgrey", "navy"), pt.size = 1, min.cutoff = 0, max.cutoff = 1)
  p_list[[cell_type]] = FeaturePlot(sc, features = cell_type, reduction = "umap", cols = c("lightgrey", "navy"), pt.size = 1, min.cutoff = 0, max.cutoff = 1)
}
p = patchwork::wrap_plots(p_list, ncol=4)
tiff(filename = file.path(param$path_out,"Modules.tiff"), width = 4400, height = 1800, res = 300)
p
dev.off()

# Annotation per cluster in a DimPlot
p = DimPlot(sc, reduction = "umap", group.by = "annotation")
tiff(filename = file.path(param$path_out,"CellType_scored_cluster.tiff"), width = 2400, height = 1200, res = 300)
p
dev.off()


### Save object with annotations
saveRDS(sc, file = file.path(param$path_out,"sc_annotated.rds"))

source(file.path(param$path_to_git,'scripts/export_data/convert.R'))


