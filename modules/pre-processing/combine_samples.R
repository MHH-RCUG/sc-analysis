##### Combining samples

### Preparation
# When merging, feature meta-data is removed by Seurat entirely; save separately for each assay except for SCT and add again afterwards
# Note: not sure whether this is still needed - discuss
assay_names = setdiff(unique(purrr::flatten_chr(purrr::map(list_names(sc), function(n) { Seurat::Assays(sc[[n]]) } ))), "SCT")

# Loop through all assays and accumulate meta data
sc_feature_metadata = purrr::map(values_to_names(assay_names), function(a) {
  # "feature_id", "feature_name", "feature_type" are accumulated for all assays and stored just once
  # This step is skipped for assays that do not contain all three types of feature information
  contains_neccessary_columns = purrr::map_lgl(list_names(sc), function(n) { 
    all(c("feature_id", "feature_name", "feature_type") %in% colnames(sc[[n]][[a]][[]])) 
  })
  
  if (all(contains_neccessary_columns)) {
    feature_id_name_type = purrr::map(sc, function(s) return(s[[a]][[c("feature_id", "feature_name", "feature_type")]]) )
    feature_id_name_type = purrr::reduce(feature_id_name_type, function(df_x, df_y) {
      new_rows = which(!rownames(df_y) %in% rownames(df_x))
      if (length(new_rows) > 0) return(rbind(df_x, df_y[new_rows, ]))
      else return(df_x)
    })
    feature_id_name_type$row_names = rownames(feature_id_name_type)
  } else {
    feature_id_name_type = NULL
  }
  
  # For all other meta-data, we prefix column names with the dataset
  other_feature_data = purrr::map(list_names(sc), function(n) {
    df = sc[[n]][[a]][[]]
    if (contains_neccessary_columns[[n]]) df = df %>% dplyr::select(-dplyr::one_of(c("feature_id", "feature_name", "feature_type"), c()))
    if (ncol(df) > 0) colnames(df) = paste(n, colnames(df), sep=".")
    df$row_names = rownames(df)
    return(df)
  })
  
  # Now join everything by row_names by full outer join
  if (!is.null(feature_id_name_type)) {
    feature_data = purrr::reduce(c(list(feature_id_name_type=feature_id_name_type), other_feature_data), dplyr::full_join, by="row_names")
  } else {
    feature_data = purrr::reduce(other_feature_data, dplyr::full_join, by="row_names")
  }
  rownames(feature_data) = feature_data$row_names
  feature_data$row_names = NULL
  
  return(feature_data)
})

# When merging, cell meta-data are merged but factors are not kept
sc_cell_metadata = suppressWarnings(purrr::map_dfr(sc, function(s){ s[[]] }) %>% as.data.frame())
sc_cell_metadata_factor_levels = purrr::map(which(sapply(sc_cell_metadata, is.factor)), function(n) {
  return(levels(sc_cell_metadata[, n, drop=TRUE]))
})

# Save variable features
# necessary for runPCA of merged, SCTransformed data
# According to https://github.com/satijalab/seurat/issues/5205 
# and https://github.com/satijalab/seurat/issues/6185
# This integrates the variable features from a list of objects (separately SCTransform)
# and ranks features by the number of data sets they are deemed variable in.
SCTVariableFeatures = SelectIntegrationFeatures(sc, nfeatures = 3000)





### Merge data
# Data for different samples can be merged if no integration is needed, 
#   for example, when samples were multiplexed on the same chip
if (param$integrate_samples[["method"]]=="merge") {
  sc = merge(x=sc[[1]], y=sc[2:length(sc)], project=param$project_id, merge.data = TRUE)
  

    # Join layers (counts, data, scaled.data) of RNA assay 
    sc[["RNA"]] = JoinLayers(sc[["RNA"]])
  
  
  message("Data values for all samples have been merged. This means that data values have been concatenated, not integrated.")
}



### Integrate data
# CCA computationally intensive.
if (param$integrate_samples[["method"]]=="integrate") {
  
  if (param$integrate_samples[["integration_function"]]=="CCAIntegration") {integration_function='cca'} 
  else if (param$integrate_samples[["integration_function"]]=="RPCAIntegration") {integration_function='rpca'}
  
  # As in https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html
  sc = RunIntegration(sc, 
                      assay=param$norm,
                      reduction=integration_function,
                      ndims=param$integrate_samples[["dimensions"]], 
                      verbose=FALSE, 
                      reference=param$integrate_samples[["reference"]], 
                      k_filter=param$integrate_samples[["k.filter"]], 
                      k_weight=param$integrate_samples[["k.weight"]], 
                      k_anchor=param$integrate_samples[["k.anchor"]],
                      k_score=param$integrate_samples[["k.score"]])
  
  if (param$norm == "RNA") { 
    # Join layers (counts, data, scaled.data) of RNA assay 
    sc[["RNA"]] = JoinLayers(sc[["RNA"]])
  }
  
  # Is this necessary???
  # Add sample as latent_vars for marker detection
  param$latent_vars = c(param$latent_vars, "orig.ident")
  
  message("Data values for all samples have been integrated.")
}



### Add metadata  
# Add feature metadata
for (a in Seurat::Assays(sc)) {
  if (a %in% names(sc_feature_metadata)) {
    sc[[a]] = Seurat::AddMetaData(sc[[a]], sc_feature_metadata[[a]][rownames(sc[[a]]),, drop=FALSE])
  }
}

# Fix cell metadata factors
for (f in names(sc_cell_metadata_factor_levels)) {
  sc[[f]] = factor(sc[[f, drop=TRUE]], levels=sc_cell_metadata_factor_levels[[f]])
}

# Add sample/dataset colours again to misc slot
sc = ScAddLists(sc, lists=list(orig.ident=param$col_samples), lists_slot="colour_lists")



### Scaling and dimensional reduction of combined data
# Re-score cell cycle effects after integration
if (param$cc_rescore_after_merge) {
  sc = CCScoring(sc=sc, genes_s=genes_s[,2], genes_g2m=genes_g2m[,2])
  if (any(is.na(sc[["S.Score"]])) | any(is.na(sc[["G2M.Score"]])))  {
    param$cc_remove=FALSE
    param$cc_remove_all=FALSE
    param$cc_rescore_after_merge=FALSE
    param$vars_to_regress = setdiff(param$vars_to_regress, c("S.Score", "G2M.Score", "CC.Difference"))
    param$latent_vars = setdiff(param$latent_vars, c("S.Score", "G2M.Score", "CC.Difference"))
  }
}

# Set default assay (will be the integrated version)
if (param$integrate_samples[["method"]]=="merge") {DefaultAssay(sc) = ifelse(param$norm=="RNA", "RNA", "SCT") } 
if (param$integrate_samples[["method"]]=="integrate") {DefaultAssay(sc) = paste0(param$norm, "integrated") } 

if (param$norm == "RNA") { 
  # Find variable features in RNA assay
  sc = Seurat::FindVariableFeatures(sc, selection.method="vst", nfeatures=3000, verbose=FALSE)
  # Scale RNA assay
  sc = Seurat::ScaleData(sc, features=rownames(sc), vars.to.regress=param$vars_to_regress, verbose=FALSE)
  
} else if (param$norm == "SCT") {
  if (param$experimental_groups == "homogene") { 
    # Re-running of SCT not recommended https://github.com/satijalab/seurat/issues/6054 
    # and https://github.com/satijalab/seurat/issues/7407
    # Set variable features
    VariableFeatures(sc) = SCTVariableFeatures
  
  } else if (param$experimental_groups == "heterogene") {
    # (Re)-run SCTransform
    min_cells_overall = max(purrr::map_int(param$feature_filter, function(f) as.integer(f[["min_cells"]])))
    sc = suppressWarnings(SCTransform(sc, 
                                      assay=param$assay_raw,
                                      vars.to.regress=param$vars_to_regress, 
                                      variable.features.n = 3000,
                                      min_cells=min_cells_overall, 
                                      verbose=FALSE, 
                                      return.only.var.genes=FALSE,
                                      method=ifelse(packages_installed("glmGamPoi"), "glmGamPoi", "poisson")))
  }
}

# Run PCA
sc = Seurat::RunPCA(sc, features=Seurat::VariableFeatures(object=sc), verbose=FALSE, npcs=min(50, ncol(sc)))

print(sc)
