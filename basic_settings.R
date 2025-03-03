### Basic settings
################################################################################
# Needs to be set


### Project specification
# Set project name
param$project_id = "Testdata"


### Dataset
# A) Either load existing object
    # Set path to rds object
    param$data = NULL
# B) Or use test dataset included in ./modules/download_test_datasets/
    # Set name
    # Does only work for 10x data automatically at the moment
    param$download_test_datasets=NULL
# C) Or load count data
    # Data frame of sample names, type, and path to folder with count matrix
    # E.g. param$path_data = data.frame(name=c("sample1","sample2"), 
    #                                   type=c("10x","10x"), 
    #                                   path=c("/filtered_feature_bc_matrix", "/filtered_feature_bc_matrix")
    param$path_data = NULL
                                 
  
# Reference
param$species="human"


### For ccc
# Set sender and receiver cell types
# Set as vector e.g. param$sender = c("Tcells", "Monocytes")
param$sender = NULL
param$receiver = NULL


### For dataset mapping
# Set path to reference rds object
param$refdata = NULL


### For cell annotation clustifyr
# Set clustifyr reference e.g. "http://cells.ucsc.edu/?ds=covid19-influenza-response", "ref_hema_microarray()" 
# from https://bioconductor.org/packages/release/data/experiment/vignettes/clustifyrdatahub/inst/doc/clustifyrdatahub.html,
# or file.path(param$path_to_git,"/references/fetal-immune-pfi_clustifyr_reference.rds")
param$clustifyr_ref = NULL
# Name of annotation column in reference dataset e.g. "cell_type" for ucsc pulled datasets
param$cluster_col_clustifyr = NULL


### For generation of reference
# Dataset origin: "ucsc", "sc_expr_atlas", or "seurat_object"
param$ref_data_orig = NULL
# Name of reference 
# Downloaded from UCSC https://cells.ucsc.edu or Single Cell Expression Atlas https://www.ebi.ac.uk/gxa/sc/experiments/
# For Single Cell Expression Atlas, name should correspond to the dataset names/identifier e.g. E-MTAB-8077
# Name will be given to reference output file
param$ref_data_name = NULL
# Folder containing downloaded datasets files
# For Single Cell Expression Atlas: 
#   - '<param$ref_data_name>-quantification-raw-files' with 
#       - '<param$ref_data_name>.aggregated_filtered_counts.mtx'
#       - '<param$ref_data_name>.aggregated_filtered_counts.mtx_cols'
#       - '<param$ref_data_name>.aggregated_filtered_counts.mtx_rows'
#   - 'ExpDesign-<param$ref_data_name>.tsv'
# For UCSC:
#   - 'exprMatrix.tsv.gz'
#   - 'meta.tsv' file
param$ref_data_path = NULL 
# Name of annotation column in reference dataset e.g. "cell_type" for ucsc pulled datasets
param$cluster_col_clustifyr = NULL


### For contamination_removal
# path to samples. Output from cellranger counts. Top level cellranger output directory containing the "filtered_gene_bc_matrices" and "raw_gene_bc_matrices" folders.
# Set as vector for multiple samples.
param$cellranger_samples = NULL
