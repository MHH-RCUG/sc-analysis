### Integrate 
# On merged data sets 
# sc = merge(x=sc[[1]], y=sc[2:length(sc)], ...)
# sc = JoinLayers(sc)
# sc = RunPCA(sc, ...)

# Split RNA layer(s)
#DefaultAssay(sc) = "RNA"

sc[["RNA"]] <- split(sc[["RNA"]], f = sc$orig.ident)

# As in https://satijalab.org/seurat/articles/integration_introduction.html
# and https://satijalab.org/seurat/articles/seurat5_integration
# Integration functions:
integration_function = c("CCAIntegration", "RPCAIntegration")
reduction_name = c("integrated.cca", "integrated.rpca")

if (param$norm == "RNA") { 
  for (i in seq(integration_function)) {
    sc = IntegrateLayers(object = sc, 
                         method = integration_function[i],
                         orig.reduction = "pca", 
                         normalization.method = "LogNormalize",
                         features = VariableFeatures(sc),
                         new.reduction = reduction_name[i],
                         verbose = FALSE)
    
    sc = RunUMAP(sc, reduction = reduction_name[i], dims=1:param$pc_n, reduction.name = paste0(gsub("integrated", "umap", reduction_name[i])))
  }
  
} else if (param$norm == "SCT") {
  #DefaultAssay(sc) = "SCT"
  # For streamlined (one-line) integrative analysis the SCT model needs to be calculated for each sample separately!
  # Recalculate for each layer again after merging and splitting
  #sc = suppressWarnings(SCTransform(sc, 
  #                                  assay=param$assay_raw,
  #                                  vars.to.regress=param$vars_to_regress, 
  #                                  variable.features.n = 3000,
  #                                  min_cells=param$feature_filter[[i]][["min_cells"]], 
  #                                  verbose=FALSE, 
  #                                  return.only.var.genes=FALSE,
  #                                  method=ifelse(packages_installed("glmGamPoi"), "glmGamPoi", "poisson")))
  # Run PCA
  sc = Seurat::RunPCA(sc, features=Seurat::VariableFeatures(object=sc), verbose=FALSE, npcs=min(50, ncol(sc)))
  
  for (i in seq(integration_function)) {
    sc = IntegrateLayers(object = sc, 
                         method = integration_function[i],
                         orig.reduction = "pca", 
                         normalization.method = "SCT",
                         features = VariableFeatures(sc),
                         new.reduction = reduction_name[i],
                         verbose = FALSE)
    
    sc = RunUMAP(sc, reduction = reduction_name[i], dims=1:param$pc_n, reduction.name = paste0(gsub("integrated", "umap", reduction_name[i])))
  }
}

# Re-join RNA layer(s)
sc[["RNA"]] = JoinLayers(sc[["RNA"]])
