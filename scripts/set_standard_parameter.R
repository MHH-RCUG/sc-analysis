### Standard parameter 
################################################################################

### set rmarkdown parameter
param$author = Sys.info()[["user"]]

### set datasets
param$ref = NULL

### Dimensional reduction
param$pc_n = 30
param$umap_k = 30

### set standard colors
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

# set dot size for umaps/tsne
param$pt_size = 0.5

### set reference
# Default is Ensembl release 98 which is the same as used for cellranger references
param$annot_version=NULL
param$annot_main=c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession")
param$mart_attributes=c(c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession"), 
                        c("chromosome_name", "start_position", "end_position", "percentage_gene_gc_content", "gene_biotype", "strand", "description"))
param$biomart_mirror=NULL





