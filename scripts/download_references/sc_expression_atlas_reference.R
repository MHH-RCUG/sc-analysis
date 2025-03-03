### Single Cell Expression Atlas reference from data downloaded from https://www.ebi.ac.uk/gxa/sc/experiments/
################################################################################

### Load libraries
library(Matrix)



### Load reference data
# Load Single Cell Expression Atlas data

# Download data from browser
# https://www.ebi.ac.uk/gxa/sc/experiments/
# Raw counts files (MatrixMarket archive) e.g. E-MTAB-8077-quantification-raw-files folder with expression matrix daata
# Experiment design file (TSV format) e.g. ExpDesign-E-MTAB-8077.tsv withh metadata

# Load the original matrix (Matrix Market format)
expr_matrix = readMM(file.path(param$ref_data_path, paste0(param$ref_data_name, "-quantification-raw-files"), paste0(param$ref_data_name, ".aggregated_filtered_counts.mtx")))
# Load the column names
expr_matrix_cols = readLines(file.path(param$ref_data_path, paste0(param$ref_data_name, "-quantification-raw-files"), paste0(param$ref_data_name, ".aggregated_filtered_counts.mtx_cols")))
# Load the row names
expr_matrix_rows = readLines(file.path(param$ref_data_path, paste0(param$ref_data_name, "-quantification-raw-files"), paste0(param$ref_data_name, ".aggregated_filtered_counts.mtx_rows")))
expr_matrix_rows = gsub("\\t.*", "", expr_matrix_rows)
# Add column and row names
colnames(expr_matrix) <- expr_matrix_cols
rownames(expr_matrix) <- expr_matrix_rows

expr_matrix = as(as.matrix(expr_matrix), "dgCMatrix")

# Read metadata
meta_data = read.table(file.path(param$ref_data_path, paste0("ExpDesign-", param$ref_data_name, ".tsv")), header=T, sep="\t", as.is=T, row.names=1)


