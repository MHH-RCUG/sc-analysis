## ---- read_data --------
The workflow can be run for pre-processed 10x Genomics, SmartSeq-2 or other data that are represented by a simple table with transcript counts per gene and cell. Similarly, a Seurat object can be loaded to inspect the stored scRNA seq data. 

<details>
  <summary>Pre-processing of 10x data with Cell Ranger</summary>
  
  We use the 10x Genomics Cell Ranger software to map 10x sequencing reads. The result is a count matrix that contains the number of unique transcripts per gene (rows) and cell (columns). To save storage space, the matrix is converted into a condensed format and described by the following 3 files:   
  <ul>
  <li>	“features.tsv.gz”: Tabular file with gene annotation (Ensembl gene ID and the gene symbol) to identify genes </li>
  <li>	“barcodes.tsv.gz”: File with cell barcodes to identify cells </li>
  <li> “matrix.mtx.gz”: A condensed version of the count matrix </li>
  </ul>
</details>

<details>
  <summary>Pre-processing of SmartSeq-2 data</summary>
  
  Sequencing reads from SmartSeq-2 (and other) experiments can be mapped with any suitable mapping software as long as a count matrix is generated that contains the number of mapped reads per gene (rows) and cell (columns). The first row of the count matrix should contain cell barcodes or names, the first column of the matrix should contain Ensembl gene IDs. 
</details>

<details>
  <summary>What is Seurat and which information is contained in a Seurat object?</summary>

  Seurat is an R-package that is used for the analysis of single-cell RNA-seq data. We read pre-processed data as described above, and convert it to a count matrix in R. In the case of 10x data, the count matrix contains the number of unique RNA molecules (UMIs) per gene and cell. In the case of SmartSeq-2 data, the count matrix contains the number of reads per gene and cell.

  In addition to the count matrix, the workflow stores additional information in the Seurat object, including but not limited to the normalized data, dimensionality reduction and cluster results.
</details>
<br>



## ---- pre-processing --------
<details>
  <summary>Why is pre-processing so important?</summary>
  
  Cells may have been damaged during sample preparation and might be only partially captured in the sequencing. In addition, free-floating mRNA from damaged cells can be encapsulated, adding to the background noise. The low-quality cells and free-floating mRNA interfere with downstream analyses. Therefore, cells and genes are filtered based on defined quality metrics. Data normalization eliminates cell-specific biases such as the absolute number of reads per cell, allowing us to systematically compare cells afterwards. Subsequent scaling corrects for the fact that genes have different lengths, allowing us to compare genes afterwards. Lastly, highly variable genes are determined, reducing computational overhead and noise from genes that are not interesting. 
</details>
<br>



## ---- qc --------
Three commonly used QC metrics include the number of unique genes detected in each cell ("`r paste0("nFeature_", param$assay_raw)`"), the total number of molecules detected in each cell ("`r paste0("nCount_", param$assay_raw)`"), and the percentage of counts that map to the mitochrondrial genome ("percent_mt"). In cell QC these covariates are filtered via thresholding as they are a good indicator of cell viability. However, it is crucial to consider the three QC covariates jointly as otherwise it might lead to misinterpretation of cellular signals.
  If ERCC spike-in controls were used, the percentage of counts mapping to them is also shown ("percent_ercc").

<details>
  <summary>Which cells can be considered to be of low quality?</summary>
    
  On the one hand, low-quality cells such as dying, degraded and damaged cells, and empty droplets will often have very few genes ("nFeature_RNA") and low total counts ("nCount_RNA"). On the other hand, cell doublets may show an aberrantly high number of genes. Since the total number of reads detected within a cell typically strongly correlates with the number of unique genes, we can focus on the number of unique genes for filtering. In addition, damaged cells often exhibit high mitochondrial ("percent_mt") or spike-in ("percent_ercc") content. As mitochondria have their own membranes, their RNA is often the last to degrade in damaged cells and can thus be found in high numbers. However, it is important to keep in mind that different cell types and cells isolated from various species may differ in their number of expressed genes and metabolism. For example, stem cells may express a higher number of unique genes, and metabolically active cells may express a higher number of mitochondrial genes. Hence, it is crucial to consider the three QC covariates jointly as otherwise it might lead to misinterpretation of cellular signals. E.g. a cell with a high fraction of mitochondrial counts might be a dying cell with a broken membrane in combination with a low number of counts, but in combination with a intermediate to high count depth it might be a cell involved in respiratory processes.  
</details>
  
<details>
  <summary>Impact of low-quality cells on downstream analyses</summary>
    
  First of all, low-quality cells of different cell types might cluster together due to similarities in their damage-induced gene expression profiles. This leads to artificial intermediate states or trajectories between subpopulations that would otherwise be distinct. Second, low-quality cells might mask relevant gene expression changes. The main differences captured by the first principal components will likely be based on cell quality rather than the underlying biology, reducing the power of dimensionality reduction. Likewise, highly variable genes might just be genes that differ between low- and high-quality cells. Lastly and equally important, the low number of total counts in low-quality cells might lead to artificial upregulation of genes due to wrong scaling effects during the normalisation process. 
</details>
<br>



## ---- normalization --------
We start by running a __standard log normalisation__, where counts for each cell are divided by the total counts for that cell and multiplied by 10,000. This is then natural-log transformed.   
  
<details>
  <summary>What do we need normalisation for?</summary>
    
  The number of raw sequencing reads per cell is influenced by technical differences in the capture, reverse transcription and sequencing of RNA molecules, particularly due to the difficulty of achieving consistent library preparation with minimal starting material. Thus, comparing gene expression between cells may reveal differences that are solely due to sampling effects. After low-quality cells were removed in the previous step, the primary goal of normalization is to remove technical sampling effects while preserving the true biological signal.
  
  Count depth scaling is the simplest and most commonly used normalization strategy. The underlying idea is that each cell initially contained an equal number of mRNA molecules, and differences arise due to sampling effects. For each cell, the number of reads per gene is divided by a cell-specific “size factor”, which is proportional to the total count depth of the cell. The resulting normalized data add up to 1 per cell, and is then typically multiplied by a factor of 10 (10,000 in this workflow).
  
  Finally, normalized data are log-transformed for three important reasons. First, distances between log-transformed expression data represent log fold changes. Log-transformation emphasizes contributions from genes with strong relative differences, for example a gene that is expressed at an average count of 50 in cell type A and 10 in cell type B rather than a gene that is expressed at an average count of 1100 in A and 1000 in B. Second, log-transformation mitigates the relation between mean and variance in the data. Lastly, log-transformation reduces that skewness of the data as many downstream analyses require the data to be normally distributed. 
</details>
<br>



## ---- cc-removal --------
How much do gene expression profiles in the dataset reflect the cell cycle phases the single cells were in? After initial normalisation, we determined the effects of cell cycle heterogeneity by calculating a score for each cell based on its expression of G2M and S phase markers. Scoring is based on the strategy described in `r Cite("10.1126/science.aad0501")`, and human gene symbols are translated to gene symbols of the species of interest using biomaRt.  
If specified, cell cycle effects can be removed during scaling. 
<details>
  <summary>How does removal of cell cycle effects affect the data?</summary>
    
  Note that removing all signal associated to cell cycle can negatively impact downstream analysis. For example, in differentiating processes, stem cells are quiescent and differentiated cells are proliferating (or vice versa), and removing all cell cycle effects can blur the distinction between these cells. An alternative approach is to remove the difference between G2M and S phase scores. This way, signals separating non-cycling and cycling cells will be maintained, while differences amongst proliferating cells will be removed. For a more detailed explanation, see the cell cycle vignette for Seurat `r Cite("https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html")`.
</details>   



## ---- normalization-method --------
Dependent on the normalisation of your choice, we either:   

  | a. Run standard functions to select __variable genes__, and __scale__ normalised gene counts. For downstream analysis it is beneficial to focus on genes that exhibit high cell-to-cell variation, that is they are highly expressed in some cells and lowly in others. To be able to compare normalised gene counts between genes, gene counts are further scaled to have zero mean and unit variance (z-score). 

<details>
  <summary>What do we need scaling for?</summary>
  
  After normalization, gene expression data can be compared between cells. However, expression of individual genes still cannot be compared. This is because genes have different lengths and, depending on the experimental set up, longer genes can be represented by a higher number of reads. To account for this effect, normalized data are further scaled using a z-transformation, resulting in the average expression of 0 and the variance of 1 for each gene across all cells. Note that additional unwanted sources of variations can be regressed out during the scaling process, such as cell cycle effects or the percentage of mitochondrial reads. 
</details>   
  
  | b. Run __SCTransform__, a new and more sophisticated normalisation method that replaces the previous functions (__normalisation, variable genes and scaling__). 

<details>
  <summary>What is __SCTransform__ special about?</summary>
  
  The standard log-transformation applied in step 1 assumes that count depth influences all genes equally. However, it has been shown that the use of a single size factor will introduce different effects on highly and lowly expressed genes (`r Cite("10.1186/s13059-019-1874-1")`). __SCTransform__ is a new statistical approach for the modelling, normalization and variance stabilization of single-cell RNA-seq data, and is an alternative to steps 1 and 3a described above. Note that __SCTransform__ has been developed for UMI count data and can therefore safely be applied to 10x but not SmartSeq-2 data. As for the scaling in step 3a, additional unwanted sources of variations can be regressed out during __SCTransform__. 
</details>   
  
Normalisation method used for this report: `r param$norm`; with additional sources of variance regressed out: `r param$vars_to_regress`.  

While raw data is typically used for statistical tests such as finding marker genes, normalised data is mainly used for visualising gene expression values. Scaled data include variable genes only, potentially without cell cycle effects, and are mainly used to determine the structure of the dataset(s) with Principal Component Analysis, and indirectly to cluster and visualise cells in 2D space.



## ---- variable_genes --------
<details>
  <summary>Why selecting variable genes?</summary>
  
  Experience shows that 1,000-2,000 genes with the highest cell-to-cell variation are often sufficient to describe the global structure of a single-cell dataset. For example, cell type-specific genes typically highly vary between cells. Housekeeping genes, on the other hand, are similarly expressed across cells and can be disregarded to differentiate between cells. Highly variable genes are typically the genes with a cell type specific expression profile, and are often the genes of interest in single-cell experiments. Housekeeping genes, with similar levels of expression across all cells, or genes with minor expression differences, might add random noise and mask relevant changes during downstream dimensionality reduction and clustering. We therefore aim to select a sensible set of variable genes that includes interesting biological signal and excludes noise. Here, the top 3,000 variable genes are selected and used for downstream analysis. 
</details>
  
<details>
  <summary>How are variable genes selected?</summary>
  
  To determine variable genes, we need to separate biological variability from technical variability. Technical variability arises especially for lowly expressed genes, where high variability corresponds to small absolute changes that we are not interested in. Here, we use the variance-stabilizing transformation (vst) method implemented in Seurat (`r Cite("10.1186/s13059-019-1874-1")`). This method first models the technical variability as a relationship between mean gene expression and variance using local polynomial regression. The model is then used to calculate the expected variance based on the observed mean gene expression. The difference between the observed and expected variance is called residual variance and likely reflects biological variability.
</details>



## ---- rel --------
To better understand the efficiency of the applied normalisation procedures, we plot the relative log expression of genes in at most `r n_cells_rle_plot` randomly selected cells per sample before and after normalisation. This type of plot reveals unwanted variation in your data. The concept is taken from `r Cite("10.1371/journal.pone.0191629")`. In brief, we remove variation between genes, leaving only variation between samples. If expression levels of most genes are similar in all cell types, sample heterogeneity is a sign of unwanted variation.

For each gene, we calculate its median expression across all cells, and then calculate the deviation from this median for each cell. For each cell, we plot the median expression (black), the interquartile range (<span style="color:lightgrey;font-weight:bold">lightgrey</span>), whiskers defined as 1.5 times the interquartile range (<span style="color:darkgrey;font-weight:bold">darkgrey</span>), and outliers (`r paste0('<span style="color:', param$col_samples, ';font-weight:bold">', param$col_samples, '</span>', collapse=', ')`)