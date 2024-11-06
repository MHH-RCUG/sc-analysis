### Generate clustifyr reference from data downloaded from https://cells.ucsc.edu,
library(clustifyr)
library(magrittr)

### Load refererence data
# Download directly into R; does not work!
#expr_matrix = data.table::fread("http://cells.ucsc.edu/<Cell Browser dataset ID>/exprMatrix.tsv.gz")
#meta_data = data.frame(data.table::fread("http://cells.ucsc.edu/<Cell Browser dataset ID>/meta.tsv"), row.names=1)

# Load downloaded data into R
# Load only a subset of cells
expr_matrix = data.table::fread(file.path(param$ref_data_path,"exprMatrix.tsv.gz"), select = c(1:300000))
meta_data = read.table(file.path(param$ref_data_path,"meta.tsv"), header=T, sep="\t", as.is=T, row.names=1)

# Format gene names
genes = expr_matrix[,1][[1]]
genes = gsub(".+[|]", "", genes)
expr_matrix = data.frame(expr_matrix[,-1], row.names=genes)

# Format cellnames
rownames(meta_data) = gsub("-", ".", rownames(meta_data))

# Subset metadata
subset_cells = colnames(expr_matrix)
meta_data = meta_data %>% dplyr::filter(row.names(meta_data) %in% subset_cells)

gc()

expr_matrix = as(as.matrix(expr_matrix), "dgCMatrix")

# Create SingleCellExperiment object
sce_ref = SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(expr_matrix)), colData = meta_data)
# Remove matrix and metadata files to free storage
rm(list= c("expr_matrix", "meta_data"))
# Normalize data
sce_ref = scuttle::logNormCounts(sce_ref)


### Generate reference
ref = object_ref(input = sce_ref, cluster_col = "cell_type")
# Remove object
rm("sce_ref")
gc()
# Save reference
saveRDS(ref, file = file.path(param$path_to_git,"references", paste0(param$ref_data_name,"_clustifyr_reference.rds")))


# Create Seurat object
#sc_ref = Seurat::CreateSeuratObject(counts = expr_matrix, project = "reference", meta.data=meta_data)
#... Optionally, perform some standard preprocessing steps in Seurat

# Generation of clustifyr reference does not work with Seurat v5 object
#ref <- seurat_ref(seurat_object = sc_ref, cluster_col = "cell_type")

