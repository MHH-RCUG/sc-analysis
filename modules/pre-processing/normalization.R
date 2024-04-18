##### Normalization

### (part1_normalization) and (part2_cc_scores)
# Calculate cell cycle score if it should be regressed out
if (param$cc_remove) {
  # Normalize data the original way
  # This is required to score cell cycle (https://github.com/satijalab/seurat/issues/1679)
  sc = purrr::map(sc, Seurat::NormalizeData, normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE)
  
  # Determine cell cycle effect per sample 
  sc = purrr::map(list_names(sc), function(n) {
    sc[[n]] = CCScoring(sc=sc[[n]], genes_s=genes_s[,2], genes_g2m=genes_g2m[,2], name=n)
    if (any(is.na(sc[[n]][["S.Score"]])) | any(is.na(sc[[n]][["G2M.Score"]]))) {
      param$cc_remove=FALSE
      param$cc_remove_all=FALSE
      param$cc_rescore_after_merge=FALSE
   }
    return(sc[[n]])
  })

  # If cell cycle effects should be removed, we first score cells 
  # The effect is then removed in the following chunk 
  # Add to vars that need to regressed out during normalization
  if (param$cc_remove_all) {
    # Remove all signal associated to cell cycle
    param$vars_to_regress = unique(c(param$vars_to_regress, "S.Score", "G2M.Score"))
    param$latent_vars = unique(c(param$latent_vars, "S.Score", "G2M.Score"))
  } else {
    # Don't remove the difference between cycling and non-cycling cells 
    param$vars_to_regress = unique(c(param$vars_to_regress, "CC.Difference"))
    param$latent_vars = unique(c(param$latent_vars, "CC.Difference"))
  }  
}

# Assign cell cycle score also already here if there is only one sample    
if (length(sc) == 1) {
  sc = purrr::map(list_names(sc), function(n) {
    sc[[n]] = CCScoring(sc=sc[[n]], genes_s=genes_s[,2], genes_g2m=genes_g2m[,2], name=n)
    if (any(is.na(sc[[n]][["S.Score"]])) | any(is.na(sc[[n]][["G2M.Score"]]))) {
      param$cc_remove=FALSE
      param$cc_remove_all=FALSE
      param$cc_rescore_after_merge=FALSE
    }
    return(sc[[n]])
  })
}



### (part3_normalization)
if (param$norm == "RNA") { 
  # For Log Normalize data
  sc = purrr::map(sc, Seurat::NormalizeData, normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE)
  
  # For integration or if only one sample
  # Multiple Log Normalized samples that should be merged, will be scaled after merging
  # vignette: https://satijalab.org/seurat/articles/integration_introduction.html
  if (param$integrate_samples[["method"]]=="integrate" | (length(sc) == 1)) {
    # Find variable features from normalised data (unaffected by scaling)
    sc = purrr::map(sc, Seurat::FindVariableFeatures, selection.method="vst", nfeatures=3000, verbose=FALSE)
    # Scale RNA assay
    sc = purrr::map(sc, Seurat::ScaleData, features=rownames(sc[["RNA"]]),
                    vars.to.regress=param$vars_to_regress, 
                    verbose=FALSE)
    # Run PCA
    sc = purrr::map(sc, Seurat::RunPCA, verbose=FALSE, npcs=min(50, ncol(sc)))
  }
  
} else if (param$norm == "SCT") {
  # For SCTransform
  # This is a new normalisation method that replaces previous Seurat functions "NormalizeData", "FindVariableFeatures", and "ScaleData". 
  # vignette: https://satijalab.org/seurat/v3.0/sctransform_vignette.html
  # vignette: https://satijalab.org/seurat/articles/integration_introduction.html
  # Normalised data end up here: sc@assays$SCT@data
  # Note: It is recommended to run SCTransform on each experiment individually before Merging or integrating. 
  #   See: https://github.com/satijalab/seurat/issues/6054
  # Note: It is not guaranteed that all genes are successfully normalised with SCTransform. 
  #   Consequently, some genes might be missing from the SCT assay. 
  #   See: https://github.com/ChristophH/sctransform/issues/27
  # Note: The performance of SCTransform can be improved by using "glmGamPoi" instead of "poisson" as method for initial parameter estimation.
  sc = purrr::map(list_names(sc), function(n) { 
    suppressWarnings(SCTransform(sc[[n]], 
                assay=param$assay_raw,
                vars.to.regress=param$vars_to_regress, 
                min_cells=param$feature_filter[[n]][["min_cells"]], 
                verbose=FALSE, 
                return.only.var.genes=FALSE,
                method=ifelse(packages_installed("glmGamPoi"), "glmGamPoi", "poisson"))) 
  })
  
    
    # For integration 
    # vignette: https://satijalab.org/seurat/articles/integration_introduction.html
    if (param$integrate_samples[["method"]]=="integrate") {
      sc = purrr::map(sc, Seurat::RunPCA, assay="SCT", verbose=FALSE, npcs=min(50, ncol(sc[[n]])))
    }
}

