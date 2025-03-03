### Generate clustifyr reference from data downloaded from https://cells.ucsc.edu,
################################################################################

### Load libraries
library(clustifyr)
library(Seurat)


### Load data 
if (param$ref_data_orig == "ucsc") {
  if (file.exists(file.path(param$ref_data_path,"exprMatrix.tsv.gz"))
      && file.exists(file.path(param$ref_data_path,"meta.tsv"))) {
    source(file.path(param$path_to_git,"scripts/download_references/ucsc_reference.R"))
  } else {
    message(paste0("One or more of the following files are missing: \n",
                   file.path(param$ref_data_path,"exprMatrix.tsv.gz"), "\n", 
                   file.path(param$ref_data_path,"exprMatrix.tsv.gz")
    ))
  }
  
} else if (param$ref_data_orig == "sc_expr_atlas") {
  if (file.exists(file.path(param$ref_data_path, paste0(param$ref_data_name, "-quantification-raw-files"), paste0(param$ref_data_name, ".aggregated_filtered_counts.mtx")))
      && file.exists(file.path(param$ref_data_path, paste0(param$ref_data_name, "-quantification-raw-files"), paste0(param$ref_data_name, ".aggregated_filtered_counts.mtx_cols")))
      && file.exists(file.path(param$ref_data_path, paste0(param$ref_data_name, "-quantification-raw-files"), paste0(param$ref_data_name, ".aggregated_filtered_counts.mtx_rows")))
      && file.exists(file.path(param$ref_data_path, paste0("ExpDesign-", param$ref_data_name, ".tsv")))) {
    source(file.path(param$path_to_git,"scripts/download_references/sc_expression_atlas_reference.R"))
  } else {
    message(paste0("One or more of the following files are missing: \n",
                   file.path(param$ref_data_path, paste0(param$ref_data_name, "-quantification-raw-files"), paste0(param$ref_data_name, ".aggregated_filtered_counts.mtx")), "\n", 
                   file.path(param$ref_data_path, paste0(param$ref_data_name, "-quantification-raw-files"), paste0(param$ref_data_name, ".aggregated_filtered_counts.mtx_cols")), "\n", 
                   file.path(param$ref_data_path, paste0(param$ref_data_name, "-quantification-raw-files"), paste0(param$ref_data_name, ".aggregated_filtered_counts.mtx_rows")), "\n", 
                   file.path(param$ref_data_path, paste0("ExpDesign-", param$ref_data_name, ".tsv"))
    ))
  }
  
} else if (param$ref_data_orig == "seurat_object") {
  if (file.exists(file.path(param$ref_data_path, paste0(param$ref_data_name, ".rds")))) {
    
  } else {
    message(paste0("The following files is missing: \n",
                   file.path(param$ref_data_path, paste0(param$ref_data_name, ".rds"))
    ))
  }
  
} else {
  message("Reference dataset source does not exist.")
}



### Create output directories
param$path_out = file.path(param$path_to_git,"references",param$ref_data_name)
if (!file.exists(param$path_out)) dir.create(param$path_out, recursive=TRUE, showWarnings=FALSE)



### Generate object
if (param$ref_data_orig == "seurat_object") {
  
  # Read Seurat object
  sc_ref = readRDS(file.path(param$ref_data_path, paste0(param$ref_data_name, ".rds")))
  
  # Convert Seurat object to SingleCellExperiment class
  sce_ref = as.SingleCellExperiment(sc_ref)
  
} else {
  # Create Seurat object
  sc_ref = Seurat::CreateSeuratObject(counts = expr_matrix, project = "reference", meta.data=meta_data, min.cells = 5, min.features = 200)
  #... Optionally, perform some standard preprocessing steps in Seurat
  saveRDS(sc_ref, file = file.path(param$path_out, paste0(param$ref_data_name,"_sc_ref.rds")))
  
  # Clean up
  rm(list= c("sc_ref"))
  suppressWarnings(gc())
  
  # Create SingleCellExperiment object
  sce_ref = SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(expr_matrix)), colData = meta_data)
  # Clean up
  rm(list= c("expr_matrix", "meta_data"))
  suppressWarnings(gc())
  # Normalize data
  sce_ref = scuttle::logNormCounts(sce_ref)
}



### Generate reference
ref = object_ref(input = sce_ref, cluster_col = param$cluster_col_clustifyr)
# Clean up
rm("sce_ref")
suppressWarnings(gc())
# Save reference
saveRDS(ref, file = file.path(param$path_out, paste0(param$ref_data_name,"_clustifyr_reference.rds")))






