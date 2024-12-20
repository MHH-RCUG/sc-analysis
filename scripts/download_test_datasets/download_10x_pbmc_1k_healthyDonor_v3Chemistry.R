# Dataset provided by 10x: 1k Peripheral blood mononuclear cells (PBMCs) from a healthy donor
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3

# Data output in a input_data subfolder of the directory where it is run 

unlink("filtered_feature_bc_matrix", recursive=T)
unlink("counts", recursive=T)
dir.create("counts", showWarnings=FALSE)
dir.create("counts/sample1", showWarnings=FALSE)

# download and untar
url = 'https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz'
curl::curl_download(url=url, destfile=basename(path=url))
untar(tarfile = basename(path=url))
unlink(basename(path=url))
file.rename("filtered_feature_bc_matrix","counts/sample1")

# get metrics summary
url = 'https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_metrics_summary.csv'
curl::curl_download(url=url, destfile='counts/metrics_summary.csv')

# write marker file
openxlsx::write.xlsx(data.frame(bcell=c("ENSG00000105369", "ENSG00000156738"),
           tcell=c("ENSG00000167286", NA),
           tcell_cd8p=c("ENSG00000153563", "ENSG00000172116"),
           nk=c("ENSG00000115523", "ENSG00000105374"),
           myeloid=c("ENSG00000101439", "ENSG00000090382"),
           monocytes=c("ENSG00000203747", NA),
           dendritic=c("ENSG00000179639", NA)),
           "counts/known_markers.xlsx")

setwd(param$path_to_git)