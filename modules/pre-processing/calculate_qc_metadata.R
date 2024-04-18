### Calculate percentages of mitochondrial and ribosomal genes as well as ERCC (qc_calculate_cells)
# Calculate percentage of counts in mitochondrial genes for each Seurat object
sc = purrr::map(sc, function(s) {
  if (!("percent_mt" %in% colnames(s@meta.data))) {
    mt_features = grep(pattern=param$mt, rownames(s), value=TRUE)
    return(Seurat::PercentageFeatureSet(s, features=mt_features, col.name="percent_mt", assay=param$assay_raw))
  } else {
    return(s)
  }
})

# Calculate percentage of counts in ribosomal genes for each Seurat object
sc = purrr::map(sc, function(s) {
  if (!("percent_ribo" %in% colnames(s@meta.data))) {
    ribo_features = grep(pattern="^RP[SL]", rownames(s), value=TRUE, ignore.case=TRUE)
    return(Seurat::PercentageFeatureSet(s, features=ribo_features, col.name="percent_ribo", assay=param$assay_raw))
  } else {
    return(s)
  }
})

# Calculate percentage of counts in ERCC for each Seurat object (if assay is available)
sc = purrr::map(sc, function(s) {
  if ("ERCC" %in% Seurat::Assays(s)) s$percent_ercc = s$nCount_ERCC/(s$nCount_ERCC + s$nCount_RNA)*100
  return(s)
})

# Combine cell metadata of the Seurat objects into one big metadata
sc_cell_metadata = suppressWarnings(purrr::map_dfr(sc, function(s) return(s[[]])) %>% as.data.frame())

# Names of all available QC metrics
cell_qc_features = c(paste0(c("nFeature_", "nCount_"), param$assay_raw), "percent_mt")
if ("percent_ercc" %in% colnames(sc_cell_metadata)) cell_qc_features = c(cell_qc_features, "percent_ercc")
cell_qc_features = values_to_names(cell_qc_features)




### Expand filter thresholds to all samples (qc_criteria_create)
# If filters were specified globally (i.e. not by sample), this chunk will copy them for each sample such that downstream filtering can work by sample
param$cell_filter = purrr::map(list_names(sc), function(s) {
  if (s %in% names(param$cell_filter)) {
    return(param$cell_filter[[s]])
  } else {
    return(param$cell_filter)
  }
})

param$feature_filter = purrr::map(list_names(sc), function(s) {
  if (s %in% names(param$feature_filter)) {
    return(param$feature_filter[[s]])
  } else {
    return(param$feature_filter)
  }
})

# List of all filter thresholds per QC metrics
cell_qc_thresholds = purrr::map(cell_qc_features, function(m) {
  tresh = purrr::map_dfr(names(param$cell_filter), function(n) {
    if (m %in% names(param$cell_filter[[n]])) {
      t = data.frame(orig.ident=n, min=param$cell_filter[[n]][[m]][1], max=param$cell_filter[[n]][[m]][2]) %>% 
        tidyr::pivot_longer(cols=2:3, names_to=c("threshold")) %>%
        dplyr::filter(!is.na(value))
      t$threshold = factor(t$threshold, levels=c("min", "max"))
      return(t)
    }
  })
})



### Calculate diverse feature caracteristics (qc_calculate_features)
# Only RNA assay at the moment
# counts_median uses sapply on the counts matrix, which converts the sparse matrix into a normal matrix
#   This might have to be adapted in future (Sparse Matrix Apply Function)
sc = purrr::map(list_names(sc), function(n) {
  # Calculate percentage of counts per gene in a cell
  counts_rna = Seurat::GetAssayData(sc[[n]], layer="counts", assay=param$assay_raw)
  total_counts = sc[[n]][[paste0("nCount_", param$assay_raw), drop=TRUE]]
  counts_rna_perc = Matrix::t(Matrix::t(counts_rna)/total_counts)*100
  
  # Calculate feature filters
  num_cells_expr = Matrix::rowSums(counts_rna >= 1)
  num_cells_expr_threshold = Matrix::rowSums(counts_rna >= param$feature_filter[[n]][["min_counts"]])
  
  # Calculate median of counts_rna_perc per gene 
  counts_median = apply(counts_rna_perc, 1, median)
  
  # Add all QC measures as metadata
  sc[[n]][[param$assay_raw]] = Seurat::AddMetaData(sc[[n]][[param$assay_raw]], data.frame(num_cells_expr, num_cells_expr_threshold, counts_median))
  return(sc[[n]])
})

