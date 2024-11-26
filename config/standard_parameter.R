### Standard parameter 
################################################################################

### Set rmarkdown parameter
param$author = Sys.info()[["user"]]



### Dataset settings
# Use test dataset
param$download_test_datasets = NULL

# Path to count matrices; data frame with name, type, and path
param$path_data = NULL

# Path to existing RDS file 
param$data = NULL
# Path to existing reference dataset for dataset mapping
param$refdata = NULL

# Set assay
# Data type ("RNA", "Spatial")
param$assay_raw = "RNA"

# Downsample data to at most n cells per sample AFTER filtering (mainly for tests)
# NULL to deactivate
param$downsample_cells_n = NULL

# Downsample all samples equally according to the smallest sample
# TRUE/FALSE
# Overwritten by downsample_cells_n
param$downsample_cells_equally = FALSE



### Filter
# Filter for cells
param$cell_filter = list(nFeature_RNA=c(20, NA), nCount_RNA=c(200, NA), percent_mt=c(0, 20))
# Filter for features
param$feature_filter = list(min_counts=1, min_cells=5)
# Samples to drop
# Cells from these samples will be dropped after initial QC
# Example: samples_to_drop = c("<name of dataset>_<names of subsample>")
param$samples_to_drop = c()
# Drop samples with too few cells
param$samples_min_cells = 10



### Normalization
# Which normalization should be used for analysis? ("RNA", "SCT")
param$norm="RNA"

# Whether or not to remove cell cycle effects
param$cc_remove = FALSE
# Should all cell cycle effects be removed, or only the difference between proliferating cells (G2M and S phase)?
# Read https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html, for an explanation
param$cc_remove_all = FALSE
# Whether or not to re-score cell cycle effects after data from different samples have been merged/integrated
param$cc_rescore_after_merge = TRUE

# Additional (unwanted) variables that will be regressed out for visualisation and clustering ("nCount_RNA", "percent_mt")
param$vars_to_regress = c()

# How to combine multiple datasets; all parameter saved in one list 
# - method = "merge" or "integrate"
# "merge": Concatenate data e.g. when samples were multiplexed on the same chip.
# "integrate": Anchors are computed for all pairs of datasets. This will give all datasets the same weight during dataset integration but can be computationally intensive
# Additional options for the "integrate" method:
#   - integration_function: "CCAIntegration" or "RPCAIntegration"
#   - dimensions: Number of dimensions to consider for integration
#   - reference: Use one or more datasets (separate by comma) as reference and compute anchors for all other datasets. Computationally faster but less accurate.
#   - use_reciprocal_pca: Compute anchors in PCA space. Even faster but less accurate. Recommended for big datasets.
#   - k_filter: How many neighbors to use when filtering anchors (default: min(200, minimum number of cells in a sample))
#   - k_weight: Number of neighbors to consider when weighting anchors (default: min(100, minimum number of cells in a sample))
#   - k_anchor: How many neighbors to use when picking anchors (default: min(10, minimum number of cells in a sample))
#   - k_score: How many neighbors to use when scoring anchors (default: min(30, minimum number of cells in a sample))
param$integrate_samples = list(method="merge", dimensions=30, k_anchor=10, k_weight=100, reference=NULL, integration_function="CCAIntegration")

# Similarity between samples ("homogene" or "heterogene")
# "heterogene" (default): If samples are biologically heterogeneous or under different treatments.
# "homogene": If samples (with roughly the same celltype composition) are technically noisy (i.e. have batch effect) with only simple shifts in mean expression.
param$experimental_groups = "homogene"


### Dimensional reduction
param$pc_n = 20
# k nearest neighbors to find clusters
# k nearest neighbors to construct the UMAP
# Scanpy uses 15 for both by default
# Seurat uses 20 for cluster_k, and 30 for umap_k by default
param$cluster_k = 20
param$umap_k = 30

# Cluster resolutions to compute; multiple values possible (comma separated); Empty vector if not needed
param$cluster_resolution_test = c(0.5, 0.7, 0.8)
# Cluster resolution to use for analysis
param$cluster_resolution = 0.6



### Set reference
# All default references are set in read_gene_annotation.R
# Here only parameter introduction that can be changed in advanced_settings
param$file_annot = NULL
param$file_cc_genes = NULL 
param$mart_dataset = NULL
param$mt = NULL
param$annot_version=NULL
param$annot_main=NULL
param$mart_attributes=NULL
param$biomart_mirror=NULL



# List of marker genes
param$file_known_markers = NULL



### Marker genes and differential expression testing
# Thresholds to define marker genes
param$marker_padj = 0.05
param$marker_log2FC = log2(2)
param$marker_pct = 0.25

# Additional (unwanted) variables to account for in statistical tests
param$latent_vars = NULL

# Contrasts to find differentially expressed genes (R data.frame or Excel file)
# Required columns:
# condition_column: Categorial column in the cell metadata; specify "orig.ident" for sample and "seurat_clusters" for cluster
# condition_group1: Condition levels in group 1, multiple levels concatenated by the plus character
#                     Empty string = all levels not in group2 (cannot be used if group2 is empty)
# condition_group2: Condition levels in group 2, multiple levels concatenated by the plus character
#                     Empty string = all levels not in group1 (cannot be used if group1 is empty)
#
# Optional columns:
# subset_column: Categorial column in the cell metadata to subset before testing (default: NA)
#                  Specify "orig.ident" for sample and "seurat_clusters" for cluster 
# subset_group: Further subset levels (default: NA)
#                 For the individual analysis of multiple levels separate by semicolons
#                 For the joint analysis of multiple levels concatenate by the plus character 
#                 For the individual analysis of all levels empty string ""
# assay: Seurat assay to test on; can also be a Seurat dimensionality reduction (default: "RNA")
# slot: In case assay is a Seurat assay object, which slot to use (default: "data")
# padj: Maximum adjusted p-value (default: 0.05)
# log2FC: Minimum absolute log2 fold change (default: 0)
# min_pct: Minimum percentage of cells expressing a gene to test (default: 0.1)
# test: Type of test; "wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"; (default: "wilcox")
# downsample_cells_n: Downsample each group to at most n cells to speed up tests (default: NA)
# latent_vars: Additional variables to account for; multiple variables need to be concatenated by semicolons; will overwrite the default by param$latent_vars (default: none).
#param$deg_contrasts = data.frame(condition_column=c("orig.ident", "orig.ident", "Phase"),
#                             condition_group1=c("pbmc_10x", "pbmc_10x", "G1"),
#                             condition_group2=c("pbmc_smartseq2_sample1", "pbmc_smartseq2_sample1", "G2M"),
#                             subset_column=c(NA, "seurat_clusters", "seurat_clusters"),
#                             subset_group=c(NA, "", "1;2"),
#                             downsample_cells_n=c(NA, 50, 30))
#param$deg_contrasts = data.frame(condition_column=c("orig.ident"),
#                             condition_group1=c("sample1"),
#                             condition_group2=c("sample2"),
#                             subset_column=c("seurat_clusters"),
#                             subset_group=c(""),
#                             downsample_cells_n=c(50))
param$deg_contrasts = NULL

# Enrichr site ("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr")
param$enrichr_site = "Enrichr"
# Enrichr libraries
param$enrichr_dbs = NULL

# Default set in read_gene_annotation.R
# Here only parameter introduction that can be changed in advanced_settings
param$enrichr_dbs=NULL

# P-value threshold for functional enrichment tests
param$enrichr_padj = 0.05

# Cell type annotation database
param$annotation_dbs = NULL


### Set standard colors
# See https://r-charts.com/color-palettes/ and https://nanx.me/ggsci/articles/ggsci.html for palette characteristics 
# Colour palette used for samples
param$col_palette_samples = "ggsci::pal_igv"
# Colour palette used for cluster
param$col_palette_clusters = "ggsci::pal_igv"
# Colour palette used for annotated cell types
param$col_palette_annotation = "ggsci::pal_ucscgb"
# Defined colours for samples
param$col_samples = NULL
param$col_samples_ref = NULL
# Defined colours for seurat_clusters
param$col_clusters = NULL
param$col_clusters_ref = NULL
# Defined colours for annotated cell types
param$col_annotation = NULL
param$col_annotation_ref = NULL
# Feature Plot colors - Highlights
param$col = "#0086b3"
# Feature Plot colors - Background
param$col_bg = "#D3D3D3"

# Set dot size for umaps/tsne
param$pt_size = 0.5


### ccc analysis
# Set sender and receiver cell types
param$sender = NULL
param$receiver = NULL

# Available methods: "connectome", "logfc", "natmi", "sca", "cellphonedb", "cytotalk", "call_squidpy", "call_cellchat", "call_connectome", "call_sca", "call_italk", "call_natmi"  
param$liana_methods = c("connectome", "logfc", "natmi", "sca", "cellphonedb")

# Threshold for liana agg rank (default: 0.01)
param$liana_agg_rank_threshold = 0.01


### For dataset mapping
# Set reference object
param$refdata = NULL

# Pre-annotated cell types; column in reference dataset
param$celltype = "annotation"   # at the moment required

# Reduction to use 'umap' or 'tsne'; must exist in ref dataset
param$reduction = 'umap'

# Predicted score threshold
param$predicted_score_threshold = 0.9
# Minimum fraction of cell with respective cell identity 
param$percent_predicted_cells_threshold = 0.1


### For cell annotation clustifyr
# Set clustifyr reference
param$clustifyr_ref = NULL
# Name of annotation column in reference dataset
param$cluster_col_clustifyr = NULL





### Other parameters
# Advanced parameter list
param_advset = list()



