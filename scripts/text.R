## ---- read_data --------
The workflow is appicable to single cell and nuclei RNAseq data pre-processed via 10x Genomics or SmartSeq-2 or for other data that are represented by a simple table with transcript counts per gene and cell. In the first step, a Seurat object of the data is generated and subsequently some metadata are added. Similarly, a Seurat object can be loaded to inspect the stored scRNA seq data. 

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
  Seurat v5 assays store data in layers. These layers can store raw, un-normalized counts (layer='counts'), normalized data (layer='data'), or z-scored/variance-stabilized (scaled) data (layer='scale.data').  
  In addition, the workflow stores additional information in the Seurat object, including but not limited to dimensionality reduction and cluster results.
</details>



## ---- pre-processing --------
<details>
  <summary>Why is pre-processing so important?</summary>
  
  Cells may have been damaged during sample preparation and might be only partially captured in the sequencing. In addition, free-floating mRNA from damaged cells can be encapsulated, adding to the background noise. The low-quality cells and free-floating mRNA interfere with downstream analyses. Therefore, cells and genes are filtered based on defined quality metrics. Data normalization eliminates cell-specific biases such as the absolute number of reads per cell, allowing us to systematically compare cells afterwards. Subsequent scaling corrects for the fact that genes have different lengths, allowing us to compare genes afterwards. Lastly, highly variable genes are determined, reducing computational overhead and noise from genes that are not interesting. 
</details>



## ---- qc --------
These commonly used QC metrics, namely the number of unique genes detected in each cell ("nFeature_"), the total number of molecules detected in each cell ("nCount_"), and the percentage of counts that map to the mitochrondrial genome ("percent_mt"), can be used to infer filter thresholds. If ERCC spike-in controls were used, the percentage of counts mapping to them is also shown ("percent_ercc").

<details>
  <summary>Which cells can be considered to be of low quality?</summary>
    
  In cell QC these covariates are filtered via thresholding as they are a good indicator of cell viability. However, it is crucial to consider the three QC covariates jointly as otherwise it might lead to misinterpretation of cellular signals.  
  
  - Cells with very high number of genes should possibly be excluded as those might be duplicates. Since the total number of reads detected within a cell typically strongly correlates with the number of unique genes, we can focus on the number of unique genes for filtering.  
  - Cells with very low count depth might be droplets with ambient RNA. If the dataset contains a high number of such events, it is advisable to perform additionally upstream correction for ambient RNA e.g. via SoupX `r knitcitations::citet("https://doi.org/10.1093/gigascience/giaa151")` as ambient RNA can confound the number of observed counts and adds background noise to the data.  
  - Low-quality cells such as dying, degraded and damaged cells will also often have very few genes ("nFeature_") and low total counts ("nCount_"). In addition, damaged cells often exhibit high mitochondrial ("percent_mt") or spike-in ("percent_ercc") content. As mitochondria have their own membranes, their RNA is often the last to degrade in damaged cells and can thus be found in high numbers. However, it is important to keep in mind that different cell types and cells isolated from various species may differ in their number of expressed genes and metabolism. For example, stem cells may express a higher number of unique genes, and metabolically active cells may express a higher number of mitochondrial genes. Hence, it is crucial to consider the three QC covariates jointly as otherwise it might lead to misinterpretation of cellular signals. E.g. a cell with a high fraction of mitochondrial counts might be a dying cell with a broken membrane in combination with a low number of counts, but in combination with a intermediate to high count depth it might be a cell involved in respiratory processes.  
</details>
  
<details>
  <summary>Impact of low-quality cells on downstream analyses</summary>
    
  First of all, low-quality cells of different cell types might cluster together due to similarities in their damage-induced gene expression profiles. This leads to artificial intermediate states or trajectories between subpopulations that would otherwise be distinct. Second, low-quality cells might mask relevant gene expression changes. The main differences captured by the first principal components will likely be based on cell quality rather than the underlying biology, reducing the power of dimensionality reduction. Likewise, highly variable genes might just be genes that differ between low- and high-quality cells. Lastly and equally important, the low number of total counts in low-quality cells might lead to artificial upregulation of genes due to wrong scaling effects during the normalisation process. 
</details>



## ---- normalization --------
Dependent on the normalisation of your choice, we either:   

  a. We perform __standard log normalisation__ by running standard functions to select __variable genes__, and __scale__ normalised gene counts.    
  b. Run __SCTransform__, a new and more sophisticated normalisation method that replaces the previous functions (__normalisation, variable genes and scaling__). 

<details>
  <summary>What do we need normalisation for?</summary>
    
  The number of raw sequencing reads per cell is influenced by technical differences in the capture, reverse transcription and sequencing of RNA molecules, particularly due to the difficulty of achieving consistent library preparation with minimal starting material. Thus, comparing gene expression between cells may reveal differences that are solely due to sampling effects. After low-quality cells were removed by filtering, the primary goal of normalization is to remove technical sampling effects while preserving the true biological signal. Hence, normalization of the raw counts accounts for differences in sequencing depth per cell.
  
  The __standard log normalisation__, a count depth scaling, is the simplest and most commonly used normalization strategy. The underlying idea is that each cell initially contained an equal number of mRNA molecules, and differences arise due to sampling effects. For each cell, the number of reads per gene is divided by a cell-specific “size factor”, which is proportional to the total count depth of the cell. The resulting normalized data add up to 1 per cell, and is then typically multiplied by a factor of 10 (10,000 in this workflow).
  Finally, normalized data are log-transformed for three important reasons. First, distances between log-transformed expression data represent log fold changes. Log-transformation emphasizes contributions from genes with strong relative differences, for example a gene that is expressed at an average count of 50 in cell type A and 10 in cell type B rather than a gene that is expressed at an average count of 1100 in A and 1000 in B. Second, log-transformation mitigates the relation between mean and variance in the data. Lastly, log-transformation reduces that skewness of the data as many downstream analyses require the data to be normally distributed. 
</details>

<details>
  <summary>What do we need scaling for?</summary>
  
  To be able to compare normalised gene counts between genes, gene counts are further scaled to have zero mean and unit variance (z-score).  
  After normalization, gene expression data can be compared between cells. However, expression of individual genes still cannot be compared. Highly expresed genes can overpower the signal of other less expresed genes with equal importance (overdispersed count values). Moreover, in SmartSeq-2 data, due to different gene length, longer genes can be represented by a higher number of reads. To account for this effect, normalized data are further scaled using a z-transformation, resulting in the average expression of 0 and the variance of 1 for each gene across all cells. Note that additional unwanted sources of variations can be regressed out during the scaling process, such as cell cycle effects or the percentage of mitochondrial reads. 
</details>   

<details>
  <summary>What is __SCTransform__ special about?</summary>
  
  The __standard log normalisation__ assumes that count depth influences all genes equally. However, it has been shown that the use of a single size factor will introduce different effects on highly and lowly expressed genes (`r Cite("10.1186/s13059-019-1874-1")`). __SCTransform__ is a new statistical approach for the modelling, normalization and variance stabilization of single-cell RNA-seq data. __SCTransform__ models the UMI counts (via a regularized negative binomial model) to remove the variation due to sequencing depth and adjusts the variance based on pooling information across genes with similar abundances and automatically regresses out sequencing depth (nCounts). Hence, __SCTransform__ can be applied to 10x (contains UMI counts) but not SmartSeq-2 data. Additional unwanted sources of variations can be regressed out during __SCTransform__. 
</details>   
  
While raw data is typically used for statistical tests such as finding marker genes, normalised data is mainly used for visualising gene expression values. Scaled data are mainly used to determine the structure of the dataset(s) with Principal Component Analysis, and indirectly to cluster and visualise cells in 2D space.



## ---- variable_genes --------
<details>
  <summary>Why selecting variable genes?</summary>
  
  For downstream analysis it is beneficial to focus on genes that exhibit high cell-to-cell variation, that is they are highly expressed in some cells and lowly in others.  
  Experience shows that 1,000-2,000 genes with the highest cell-to-cell variation are often sufficient to describe the global structure of a single-cell dataset. For example, cell type-specific genes typically highly vary between cells. Housekeeping genes, on the other hand, are similarly expressed across cells and can be disregarded to differentiate between cells. Highly variable genes are typically the genes with a cell type specific expression profile, and are often the genes of interest in single-cell experiments. Housekeeping genes, with similar levels of expression across all cells, or genes with minor expression differences, might add random noise and mask relevant changes during downstream dimensionality reduction and clustering. We therefore aim to select a sensible set of variable genes that includes interesting biological signal and excludes noise. Here, the top 3,000 variable genes are selected and used for downstream analysis. 
</details>
  
<details>
  <summary>How are variable genes selected?</summary>
  
  To determine variable genes, we need to separate biological variability from technical variability. Technical variability arises especially for lowly expressed genes, where high variability corresponds to small absolute changes that we are not interested in.  
  __standard log normalisation__:  
  Here, we use the variance-stabilizing transformation (vst) method implemented in Seurat (`r Cite("10.1186/s13059-019-1874-1")`). This method first models the technical variability as a relationship between mean gene expression and variance using local polynomial regression. The model is then used to calculate the expected variance based on the observed mean gene expression. The difference between the observed and expected variance is called residual variance and likely reflects biological variability.  
  __SCTransform__:  
  Also __SCTransform__ returns a list of top 3,000 variable genes.
</details>



## ---- regress_out --------
Next, we estimate the impact of covariates, such as the percentage of mitochondrial gene expression and cell cycle phase, to evaluate reliability of the applied normalization and scaling method, and determine possibly variables to regress out and integration method.
  
<details>
  <summary>What are unwanted sources of variations that could be regressed out?</summary>
  
  There can be sources of uninteresting variation in the data that is often specific to the dataset. Hence, checking and removing unwanted variation is another important step in pre-processing so that those artifacts do not drive clustering in the downstream analysis.  
    
  Variables provided in vars.to.regress are individually regressed against each feature, and the resulting residuals are then scaled and centered. This step allows controling for factors that may bias clustering.  
    
  __Count depth effect__  
  Although normalization scales count data to render gene counts comparable between cells, a count depth effect often remains in the data as no scaling method can infer the expression values of genes that were not detected. This count depth effect can be both a biological and a technical artefact (due to poor sampling).  
  __SCTransform__ regresses out variation due to the number of UMIs by default.  
    
  __Cell cycle score__  
  Cell cycle phase may be a source of significant variation in the data. It is essential to check, how much the gene expression profiles in the data reflect the cell cycle phases the single cells were in.
    
  __Mitochondrial gene expression__  
  mitochondrial gene expression is another factor which can greatly influence clustering. However, if the differences in mitochondrial gene expression represent a biological phenomenon that may help to distinguish cell clusters, then we advise not regressing the mitochondrial expression.  
</details>



## ---- dimensional_reduction --------
A single-cell dataset of 20,000 genes and 5,000 cells has 20,000 dimensions. At this point of the analysis, we have already reduced the dimensionality of the dataset to 3,000 variable genes. The biological manifold however can be described by far fewer dimensions than the number of (variable) genes, since expression profiles of different genes are correlated if they are involved in the same biological process. Dimension reduction methods aim to find these dimensions. There are two general purposes for dimension reduction methods: to __summarize__ a dataset, and to __visualize__ a dataset. 

We use Principal Component Analysis (PCA) to __summarize__ a dataset, overcoming noise and reducing the data to its essential components. Later, we use Uniform Manifold Approximation and Projection (UMAP) to __visualize__ the dataset, placing similar cells together in 2D space, see below. 

<details>
  <summary>PCA in a nutshell</summary>
  
  Principal Component Analysis is a way to summarize a dataset and to reduce noise by averaging similar gene expression profiles. The information for correlated genes is compressed into single dimensions called principal components (PCs) and the analysis identifies those dimensions that capture the largest amount of variation. This process gives a more precise representation of the patterns inherent to complex and large datasets.
    
  In a PCA, the first PC captures the greatest variance across the whole dataset. The next PC captures the greatest remaining amount of variance, and so on. This way, the top PCs are likely to represent the biological signal where multiple genes are affected by the same biological processes in a coordinated way. In contrast, random technical or biological noise that affects each gene independently are contained in later PCs. Downstream analyses can be restricted to the top PCs to focus on the most relevant biological signal and to reduce noise and unnecessary computational overhead. 
</details>



## ---- determinig_dimensionality --------
PCs include biological signal as well as noise, and we need to determine the number of PCs with which we include as much biological signal as possible and as little noise as possible. The following "Elbow plot" is designed to help us make an informed decision. 

<details>
  <summary>How do we determine the number of PCs for downstream analysis?</summary>
  
  The Elbow plot shows PCs ranked based on the percentage of variance they explain. The top PC captures the greatest variance across cells and each subsequent PC represents decreasing levels of variance. By visual inspection of the Elbow plot, we try to find the point at which we can explain most of the variance across cells. Commonly, the top 10 PCs are chosen. It may be helpful to repeat downstream analyses with a different number of PCs, although the results often do not differ dramatically. Note that it is recommended to rather choose too many than too few PCs.   
</details>



## ---- cc-removal --------
<details>
  <summary>How does removal of cell cycle effects affect the data?</summary>
    
  Note that removing all signal associated to cell cycle can negatively impact downstream analysis. Cell cycle signals can be informative of the biology. For example, in differentiating processes, stem cells are quiescent and differentiated cells are proliferating (or vice versa), and removing all cell cycle effects can blur the distinction between these cells. Moreover, biological signals must be understood in context. Dependencies exist between diffferent cellular processes within the organism. Thus, correcting for one process may unintentionally mask the signal of another. Vice versa, due to a correlation of cell size and some transcriptomic effect attributed to the cell cycle (McDavid et al, 2016), correcting for cell size via normalization, also partially corrects for cell cycle effects in scRNA‐seq data.  
  
  An alternative approach is to remove the difference between G2M and S phase scores. This way, signals separating non-cycling and cycling cells will be maintained, while differences in cell cycle phase amongst proliferating cells (which are often uninteresting), will be removed. For a more detailed explanation, see the cell cycle vignette for Seurat `r Cite("https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html")`. 
</details>   
  
<details>
  <summary>How are cycle effects regressed out?</summary>
  
  Prior to calling the function CellCycleScoring the data need to be standard log normalized, since counts need to be comparable between cells. 
  Once the data is normalized for sequencing depth, a score can be calculated for each cell based on its expression of G2M and S phase markers. Scoring is based on the strategy described in `r Cite("10.1126/science.aad0501")`, and human gene symbols are translated to gene symbols of the species of interest using biomaRt.  
  If specified, cell cycle effects can be removed during scaling. In the process, for each gene, Seurat models the relationship between gene expression and the S and G2M cell cycle scores. The scaled residuals of this model represent a ‘corrected’ expression matrix. 
</details>   



## ---- sample_combination --------
Whether samples should be combined via merging or integrating is very context dependent.
It is advisable first performing a dimensionaly reduction on the merged dataset and only proceed to integration if common cell types separate due to batch effect.
  
<details>
  <summary>Why integration?</summary>
  
  If the cells cluster by covariants, such as donors, batches, or condition, integrative analysis can help to remove biological differences, match shared cell types and states across datasets, and improve clustering and downstream analyses.   
  The goal of integration is to find corresponding cell states across samples. Vice versa, integration relys on the presence of at least a subset of corresponding cell states that are shared across datasets. However, given the increased degrees of freedom of non‐linear data integration approaches, this approach might lead to over‐correction of data, especially when a large proportion of cells are non-overlapping across datasets, and loss of biological resolution.
</details>  
  
<details>
  <summary>How does integration work?</summary>
  
  Integration uses the shared highly variable genes from each condition and then allignes the conditions to overlay cells that are similar between samples. The integrated values are a non-linear transformation of scale.data and saved in an additional layer. Moreover, the Seurat v5 integration procedure returns a dimensional reduction (i.e. integrated.cca) of the data that captures the shared sources of variance across multiple layers and can be used for visualization and unsupervised clustering. The harmonization of LogNormalized and SCT expression values across datasets is overall similar (see https://satijalab.org/seurat/articles/seurat5_integration and https://satijalab.org/seurat/articles/integration_introduction.html). 
</details>  
  
<details>
  <summary>The workflow in detail</summary>
  
  The worklows for the different normalization and sample combination methods are as followed:  
  <ol>
  <li>LogNormalize data  
  1a) LogNormalize - (CCsoring) - Merge - (Rerun CCsoring) - Scale - PCA  
  1b) LogNormalize - (CCsoring) - Scale - PCA - Integrate - (Rerun CCsoring) - Scale - PCA</li>  
    
  <li>SCTansform data  
  with homogene experimental groups  
  2a) (LogNormalize - CCsoring - ) SCTansform - Merge - (Rerun CCsoring) - PCA  
  2b) (LogNormalize - CCsoring - ) SCTansform - PCA - Integrate - (Rerun CCsoring) - PCA    
  with heterogene experimental groups  
  2c) (LogNormalize - CCsoring - ) SCTansform - Merge - (Rerun CCsoring) - Rerun SCTansform - PCA  
  2d) (LogNormalize - CCsoring - ) SCTansform - PCA - Integrate - (Rerun CCsoring) - Rerun SCTansform - PCA</li>  
  </ol>
    
  __Merging LogNormalized data (1a)__  
  LogNormalized data are merged based according to https://satijalab.org/seurat/archive/v4.3/merge. Merging will combine the Seurat objects based on the raw count matrices and also merge the normalized data matrices. Any previously scaled data matrices are erased. As the results of this layer dependent on the expression across all cells in an object, the previous values from one object before merging with another are no longer valied. Therefore, the layer is removed and the scaling needs to be re-run after merging objects.  
    
  __Merging SCTransform data (2a, 2c)__  
  The goal of SCTransform is to perform normalization within an experiment by learning the model of technical noise (within the experiment). Running SCTransform standardizes the expression of each gene in each cell relative to all the other cells in the input matrix.   
  __2a)__ If you have multiple samples that are technically noisy and you want to remove batch effects that are characterized by simple shifts in mean expression, it is recommend to run SCTransform on each sample individually. However, the samples should have roughly the same celltype compositions. The merged object has then multiple SCT models and the genes in the scale.data are the intersected genes.  
  __2c)__ However, scale.data calculated based on multiple different SCTransform models are not comparable. Performing SCTransform separately on very diverging samples results in loss of the baseline differences between them (similar to a gene-wise scaling before merging). Hence, if samples are biologically heterogeneous or under different treatments and, therefore, the SCTransform models of each separate sample very diverging, SCTranform should be (re-)run on merged samples. This way, SCTransform will learn one model using the median sequencing depth of the entire dataset to calculate the corrected counts.  
    
  __Integration (1b, 2b, 2d)__  
  As integration uses the shared highly variable genes, prior to integration LogNormalize amd identification of variable genes or SCTansform is required. Depending on the integration methode, dementional reduction of the data is needed additionally (e.g. for reciprocal PCA).  
    
  In most cases, a re-scoring of the cell cycle state is recommanded or at least not disadvantageous.  
</details>  



## ---- integration_methods --------
This workflow implements the following integration methods offered by Seurat offers:   
- Anchor-based CCA integration (method=CCAIntegration)  
- Anchor-based RPCA integration (method=RPCAIntegration)  
  
Anchor-based integration can be applied to either log-normalized or SCTransform-normalized datasets and can be also performed as reference-based integration where one or a subset of the datasets are listed as a ‘reference’. Reference-based integration, reduces computational time.  
  
<details>
  <summary>Canonical Correlation Analysis (CCA)</summary>
  
  CCA-based integration enables indentification of conserved cell types across datasets and integrative analysis when gene expression is shifted e.g. due to experimental conditions, disease states, or even species affiliation. If cell types are present in only one dataset, the cells will appear as a separate sample-specific cluster.  
    
  Integration comprises the following steps:  
  1. CCA, a form of PCA, identifies shared sources of variation (using the 3000 most variant genes from each sample) between the samples.  
  2. For each cell in one sample, nearest neighbors (MNNs) based on gene expression values are identified in the other sample(s). If the same match results from reciprical analysis the cell pair forms an anchor.  
  3. The similarity between anchor pairs, i.e. whether the adjacent cells also form corresponding anchor pairs, is assed and incorrect anchors are removed.  
  4. The datasets are integrated by transforming the cell expression values using anchors and corresponding scores.  
</details>  
  
<details>
  <summary>Reciprocal PCA (RPCA)</summary>
  
  For RPCA, instead of identifying shared sources of variation using variant genes, each dataset is projected into the other’s PCA space. Afterwards, similarly to CCa, mutual neighbors and anchors are identified.  
  RPCA-based integration is significantly faster and more conservative, i.e. cells in different biological states are less likely to align. Hence, RPCA is suitable for integration of multiple datasets, of datasets with little overlay of cell types, or datasets originate from the same platform. The strength of alignment can be elevated by increasing the k.anchor parameter (default 5).  
</details>  