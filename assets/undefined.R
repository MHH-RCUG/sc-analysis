

sc = readRDS("/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/input_data/HFO/sc.rds")
scR = readRDS("/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/input_data/BGHFO/sc.rds")




# Set groups
Idents(sc) = "groups"
sc = RenameIdents(sc, 'SeuratProject' = "HFO")
sc$groups = Idents(sc)

sc$annotation = gsub("_", "-", sc$annotation)


Seurat::Idents(scR) = "orig.ident"
scR = Seurat::RenameIdents(scR, 
                   'N3774' = "BG_HFO", 
                   'N3919m' = "BG_HFO")
scR$groups = Seurat::Idents(scR)
Seurat::Idents(scR) = "seurat_clusters"

Seurat::Idents(scR) = "seurat_clusters"
scR = Seurat::RenameIdents(scR, 
                  '1' = "AE", 
                  '2' = "MES", 
                  '3' = "CM", 
                  '15' = "CM", 
                  '4' = "MKs", 
                  '12' = "MKs",
                  '16' = "MKs",
                  '5' = "PE/ST", 
                  '6' = "PE/ST", 
                  '13' = "PE/ST", 
                  '14' = "PE/ST", 
                  '7' = "Ery",
                  '9' = "VE",
                  '10' = "Mo/Mo",
                  '11' = "HPC",
                  '17' = "EPP",
                  '18' = "PFE",
                  '19' = "ELC",
                  '8' = "EC")
scR$annotation = Seurat::Idents(scR)
Seurat::Idents(scR) = "seurat_clusters"

saveRDS(sc, "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/input_data/HFO/sc.rds")
saveRDS(scR, "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/input_data/BGHFO/sc.rds")



param$col_samples = c('reference' = "#5050FFFF", 'query' = "#CE3D32FF")
param$col_samples[[1]]



#renv::use_python(type = "virtualenv", name = file.path(param$path_to_git,"env/basic/virtualenvs/r-reticulate"))
#reticulate::use_virtualenv(file.path(param$path_to_git,"env/basic/virtualenvs/r-reticulate"))


