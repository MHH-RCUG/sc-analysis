### Advanced settings
################################################################################

### Dataset settings
# Set assay
# Data type ("RNA", "Spatial")
param_advset$assay_raw = NULL

# Downsample data to at most n cells per sample AFTER filtering (mainly for tests)
# NULL to deactivate
param_advset$downsample_cells_n = NULL

# Downsample all samples equally according to the smallest sample
# TRUE/FALSE
# Overwritten by downsample_cells_n
param_advset$downsample_cells_equally = NULL



### Filter
# Filter for cells
#param_advset$cell_filter = list(nFeature_RNA=c(20, NA), nCount_RNA=c(200, NA), percent_mt=c(0, 50))
param_advset$cell_filter = NULL
# Filter for features
#param_advset$feature_filter = list(min_counts=1, min_cells=5)
param_advset$feature_filter = NULL
# Samples to drop
# Cells from these samples will be dropped after initial QC
# Example: samples_to_drop = c("<name of dataset>_<names of subsample>")
param_advset$samples_to_drop = NULL
# Drop samples with too few cells
param_advset$samples_min_cells = NULL



### Normalization
# Which normalization should be used for analysis? ("RNA", "SCT")
param_advset$norm=NULL

# Whether or not to remove cell cycle effects
param_advset$cc_remove = NULL
# Should all cell cycle effects be removed, or only the difference between proliferating cells (G2M and S phase)?
# Read https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html, for an explanation
param_advset$cc_remove_all = NULL
# Whether or not to re-score cell cycle effects after data from different samples have been merged/integrated
param_advset$cc_rescore_after_merge = NULL

# Additional (unwanted) variables that will be regressed out for visualisation and clustering ("nCount_RNA", "percent_mt")
param_advset$vars_to_regress = NULL

# How to combine multiple datasets (method = "merge" or "integrate")
# "merge": Concatenate data e.g. when samples were multiplexed on the same chip.
# "integrate": Anchors are computed for all pairs of datasets. This will give all datasets the same weight during dataset integration but can be computationally intensive
param_advset$integrate_samples = NULL

if (!is.null(param_advset$integrate_samples)) {
  if (param_advset$integrate_samples[["method"]]=="integrate") {
  # Additional options for the "integrate" method:
  #   - integration_function: "CCAIntegration" or "RPCAIntegration"
  #   - dimensions: Number of dimensions to consider for integration
  #   - reference: Use one or more datasets (separate by comma) as reference and compute anchors for all other datasets. Computationally faster but less accurate.
  #   - use_reciprocal_pca: Compute anchors in PCA space. Even faster but less accurate. Recommended for big datasets.
  #   - k_filter: How many neighbors to use when filtering anchors (default: min(200, minimum number of cells in a sample))
  #   - k_weight: Number of neighbors to consider when weighting anchors (default: min(100, minimum number of cells in a sample))
  #   - k_anchor: How many neighbors to use when picking anchors (default: min(5, minimum number of cells in a sample))
  #   - k_score: How many neighbors to use when scoring anchors (default: min(30, minimum number of cells in a sample))
  param_advset$integrate_samples = list(dimensions=30, k_anchor=20, reference=NULL, integration_function="CCAIntegration")
  # Similarity between samples ("homogene" or "heterogene")
  param_advset$experimental_groups = "homogene"
  }
}



### Dimensional reduction
param_advset$pc_n = NULL
# k nearest neighbors to find clusters
# k nearest neighbors to construct the UMAP
# Scanpy uses 15 for both by default
# Seurat uses 20 for cluster_k, and 30 for umap_k by default
param_advset$cluster_k = NULL
param_advset$umap_k = NULL

# Cluster resolutions to compute; multiple values possible (comma separated); Empty vector if not needed
param_advset$cluster_resolution_test = NULL
# Cluster resolution to use for analysis
param_advset$cluster_resolution = NULL



### Set reference
param_advset$file_annot = NULL
param_advset$file_cc_genes = NULL 
# Default is Ensembl release 98 which corresponds to 2020-A reference package of 10x Genomics Cell Ranger
# Ensembl release 110 which corresponds to 2024-A reference package of 10x Genomics Cell Ranger
# ATTENTION: Translation cc genes between human and mouse does not work; Error in getLDS() since version 105 (https://github.com/grimbough/biomaRt/issues/66)
# Means versions > 105 do not work for mouse
param_advset$annot_version=NULL
#param$annot_main=c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession")
param_advset$annot_main=NULL
#param_advset$mart_attributes=c(c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession"), 
#                        c("chromosome_name", "start_position", "end_position", "percentage_gene_gc_content", "gene_biotype", "strand", "description"))
param_advset$mart_attributes=NULL
param_advset$biomart_mirror=NULL

# List of marker genes
param_advset$file_known_markers = NULL



### Marker genes and differential expression testing
# Thresholds to define marker genes
param_advset$marker_padj = NULL
param_advset$marker_log2FC = NULL
param_advset$marker_pct = NULL

# Additional (unwanted) variables to account for in statistical tests
param_advset$latent_vars = NULL

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
#param_advset$deg_contrasts = data.frame(condition_column=c("orig.ident", "orig.ident", "Phase"),
#                             condition_group1=c("pbmc_10x", "pbmc_10x", "G1"),
#                             condition_group2=c("pbmc_smartseq2_sample1", "pbmc_smartseq2_sample1", "G2M"),
#                             subset_column=c(NA, "seurat_clusters", "seurat_clusters"),
#                             subset_group=c(NA, "", "1;2"),
#                             downsample_cells_n=c(NA, 50, 30))
#param_advset$deg_contrasts = data.frame(condition_column=c("orig.ident"),
#                                 condition_group1=c("sample1"),
#                                 condition_group2=c("sample2"),
#                                 subset_column=c("seurat_clusters"),
#                                 subset_group=c(""),
#                                 downsample_cells_n=c(50))
param_advset$deg_contrasts = NULL

# Enrichr site ("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr")
param_advset$enrichr_site = NULL

# P-value threshold for functional enrichment tests
param_advset$enrichr_padj = NULL

# Enrichr libraries
param_advset$enrichr_dbs = NULL



### Set colors
# Colour palette used for samples
param_advset$col_palette_samples = NULL
# Colour palette used for cluster
param_advset$col_palette_clusters = NULL
# Colour palette used for annotated cell types
param_advset$col_palette_annotation = NULL
# Defined colours for samples
param_advset$col_samples = NULL
param_advset$col_samples_ref = NULL
# Defined colours for seurat_clusters
param_advset$col_clusters = NULL
param_advset$col_clusters_ref = NULL
# Defined colours for annotated cell types
param_advset$col_annotation = NULL
param_advset$col_annotation_ref = NULL
# Feature Plot colors - Highlights
param_advset$col = NULL
# Feature Plot colors - Background
param_advset$col_bg = NULL

# Set dot size for umaps/tsne
param_advset$pt_size = NULL




