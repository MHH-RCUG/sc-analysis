### UCSC reference from data downloaded from https://cells.ucsc.edu,
################################################################################

### Load libraries
library(magrittr)



### Load reference data
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

suppressWarnings(gc())

expr_matrix = as(as.matrix(expr_matrix), "dgCMatrix")

