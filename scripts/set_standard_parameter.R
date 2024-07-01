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
param$cell_filter = list(nFeature_RNA=c(20, NA), nCount_RNA=c(200, NA), percent_mt=c(0, 50))
# Filter for features
param$feature_filter = list(min_counts=1, min_cells=5)
# Samples to drop
# Cells from these samples will be dropped after initial QC
# Example: samples_to_drop = c("<name of dataset>_<names of subsample>")
param$samples_to_drop = c()
# Drop samples with too few cells
param$samples_min_cells = 10


### Normalisation
# Which normalisation should be used for analysis? ("RNA", "SCT")
param$norm="RNA"

# Whether or not to remove cell cycle effects
param$cc_remove = FALSE
# Should all cell cycle effects be removed, or only the difference between profilerating cells (G2M and S phase)?
# Read https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html, for an explanation
param$cc_remove_all = FALSE
# Whether or not to re-score cell cycle effects after data from different samples have been merged/integrated
param$cc_rescore_after_merge = TRUE

# Additional (unwanted) variables that will be regressed out for visualisation and clustering ("nCount_RNA", "percent_mt")
param$vars_to_regress = c()

# Similarity between samples ("homogene" or "heterogene")
param$experimental_groups = "homogene"

# How to combine multiple datasets (method = "merge" or "integrate")
# "merge": Concatenate data e.g. when samples were multiplexed on the same chip.
# "integrate": Anchors are computed for all pairs of datasets. This will give all datasets the same weight during dataset integration but can be computationally intensive
#
# Additional options for the "integrate" method:
#   - integration_function: "CCAIntegration" or "RPCAIntegration"
#   - dimensions: Number of dimensions to consider for integration
#   - reference: Use one or more datasets (separate by comma) as reference and compute anchors for all other datasets. Computationally faster but less accurate.
#   - use_reciprocal_pca: Compute anchors in PCA space. Even faster but less accurate. Recommended for big datasets.
#   - k_filter: How many neighbors to use when filtering anchors (default: min(200, minimum number of cells in a sample))
#   - k_weight: Number of neighbors to consider when weighting anchors (default: min(100, minimum number of cells in a sample))
#   - k_anchor: How many neighbors to use when picking anchors (default: min(5, minimum number of cells in a sample))
#   - k_score: How many neighbors to use when scoring anchors (default: min(30, minimum number of cells in a sample))
#
# e.g. param$integrate_samples=list(method="integrate", dimensions=30, k_anchor=20, reference=NULL, integration_function="CCAIntegration")
param$integrate_samples = list(method="merge")



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
param$file_annot = NULL
param$file_cc_genes = NULL 
# Default is Ensembl release 98 which corresponds to 2020-A reference package of 10x Genomics Cell Ranger
param$annot_version=98
param$annot_main=c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession")
param$mart_attributes=c(c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession"), 
                        c("chromosome_name", "start_position", "end_position", "percentage_gene_gc_content", "gene_biotype", "strand", "description"))
param$biomart_mirror=NULL

# List of marker genes
param$file_known_markers = NULL


param$vars_to_regress = c()
param$latent_vars = c()



### Set standard colors
# Colour palette used for samples
param$col_palette_samples = "ggsci::pal_lancet"
# Colour palette used for cluster
param$col_palette_clusters = "ggsci::pal_igv"
# Colour palette used for annotated cell types
param$col_palette_annotation = "ggsci::pal_igv"
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


### Other parameters
# Advanced parameter list
param_advset = list()



