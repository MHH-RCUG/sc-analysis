### Integrate 
# On merged data sets 
# sc = merge(x=sc[[1]], y=sc[2:length(sc)], ...)
# sc = JoinLayers(sc)
# sc = RunPCA(sc, ...)

# Split RNA layer(s)
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
    sc = RunUMAP(sc, reduction = reduction_name[i], dims = 1:30, reduction.name = paste0(gsub("integrated", "umap", reduction_name[i])))
  }
  
} else if (param$norm == "SCT") {
  for (i in seq(integration_function)) {
    sc = IntegrateLayers(object = sc, 
                         method = integration_function[i], 
                         normalization.method = "SCT",  
                         features = VariableFeatures(sc),
                         new.reduction = reduction_name[i], 
                         verbose = FALSE)
    
    sc = RunUMAP(sc, reduction = reduction_name[i], dims = 1:30, reduction.name = paste0(gsub("integrated", "umap", reduction_name[i])))
  }
}

# Re-join RNA layer(s)
sc[["RNA"]] = JoinLayers(sc[["RNA"]])
