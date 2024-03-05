# Minimal dataset for testing.
#
# The 10x dataset "1k Peripheral blood mononuclear cells (PBMCs)" is used to create two artifical samples with 500 and 700 cells and 8000 genes each.
# Dataset: https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3

# Data output in a input_data subfolder of the directory where it is run 
# Create output directories
if (!file.exists("input_data")) dir.create("input_data", recursive=TRUE, showWarnings=FALSE)
setwd(file.path(param$path_to_git,"input_data"))

unlink("counts", recursive=T)
dir.create("counts", showWarnings=FALSE)
dir.create("counts/sample1", showWarnings=FALSE)
dir.create("counts/sample2", showWarnings=FALSE)
unlink("filtered_feature_bc_matrix", recursive=T)

# Download dataset
url = 'https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz'
curl::curl_download(url=url, destfile=basename(path=url))
untar(tarfile = basename(path=url))
unlink(basename(path=url))

# Read 10x data
tenx_data = Seurat::Read10X("filtered_feature_bc_matrix")
features = read.table(file.path("filtered_feature_bc_matrix", "features.tsv.gz"), sep="\t", header=FALSE)

# Discard genes that are not expressed at all
keep_genes = Matrix::rowSums(tenx_data) > 1
tenx_data = tenx_data[keep_genes, ]
features = features[keep_genes, ]

# Create two samples with 500 and 700 cells and 8000 genes
set.seed(11)
cells = sample(colnames(tenx_data), 1200)
cells_sample1 = cells[1:500]
cells_sample2 = cells[501:1200]

set.seed(11)
feature_rows_sampled = sample(1:nrow(features), 8000)
feature_rows_sampled = sort(feature_rows_sampled)

sample1_data = tenx_data[feature_rows_sampled, cells_sample1]
sample2_data = tenx_data[feature_rows_sampled, cells_sample2]

# Write to disk
mh = file.path("counts", "sample1", "matrix.mtx")
Matrix::writeMM(sample1_data, file=mh)
R.utils::gzip(mh, overwrite=TRUE)

bh = gzfile(file.path("counts", "sample1", "barcodes.tsv.gz"), open="wb")
write(colnames(sample1_data), file=bh)
close(bh)

fh = gzfile(file.path("counts", "sample1", "features.tsv.gz"), open="wb")
write.table(features[feature_rows_sampled, ], file=fh, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
close(fh)

mh = file.path("counts", "sample2", "matrix.mtx")
Matrix::writeMM(sample2_data, file=mh)
R.utils::gzip(mh, overwrite=TRUE)

bh = gzfile(file.path("counts", "sample2", "barcodes.tsv.gz"), open="wb")
write(colnames(sample2_data), file=bh)
close(bh)

fh = gzfile(file.path("counts", "sample2", "features.tsv.gz"), open="wb")
write.table(features[feature_rows_sampled, ], file=fh, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
close(fh)

unlink("filtered_feature_bc_matrix", recursive=T)

# write marker file
openxlsx::write.xlsx(data.frame(bcell=c("ENSG00000105369", "ENSG00000156738"),
                                tcell=c("ENSG00000167286", NA),
                                tcell_cd8p=c("ENSG00000153563", "ENSG00000172116"),
                                nk=c("ENSG00000115523", "ENSG00000105374"),
                                myeloid=c("ENSG00000101439", "ENSG00000090382"),
                                monocytes=c("ENSG00000203747", NA),
                                dendritic=c("ENSG00000179639", NA)),
                     "known_markers.xlsx")

setwd(param$path_to_git)