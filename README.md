# UNDER DEVELOPMENT !!!

# Preamble  
  
**At moment, the workflow is under development and only running and tested on RStudio server at Medizinische Hochschule Hannover (MHH)** 
**The workflow is set up for the rcug user on hpc-rc11 on HPC at MHH!**  
**Before starting, see docs/Run_scrnaseq_workflow.docx for an instruction how to run the scrnaseq workflow on our system!**    
- includes renv, python env, and instruction for usage  
- includes installation documentation and instruction for usage on our system  
- includes documentation of package versions  
  
  
----------------------------------------------------------
  

# Table of contents
* [Introduction](#introduction) 
* [News](#news) 
* [Workflow summary](#workflow_summary)
* [Installation](#installation) 
* [Quick start](#quick_start)
* [Usage](#usage) 
* [Output](#output) 
* [Documentation](#documentation)
* [Credits](#credits)
* [Contributions and support](#contributions_and_support)
* [Citation](#citation)

# Introduction
<a name="introduction"/>

**sc-analysis** is a bioinformatics analysis workflow for single-cell RNA-seq analysis. The workflow is based on Seurat and the scrnaseq workflow (status May 4, 2022; Code on Zenodo: https://zenodo.org/record/7849063) created from Dresden-concept Genome Center URL "https://genomecenter.tu-dresden.de". In addition, the workflow utilizes diverse R packages for data processing, visualization, and downstream analysis.   

The workflow is composed of modules. While some modules execute core sc-RNA seq data processing steps, others modules are optional providing supporting functionalities or allowing further specific downstream analyses, such as dataset mapping or trajectory analysis. For more information regarding the scope of each module, refer to the module descriptions below. 
   
The workflow is appicable to single cell and nuclei RNAseq data pre-processed via 10x Genomics or SmartSeq-2 or for other data that are represented by a simple table with transcript counts per gene and cell. Similarly, a Seurat object can be loaded to inspect the stored scRNA seq data and perform downstream analysis.  

Example reports for the test dataset '10x_pbmc_small_split2samples' can be found in the GitHub repository under output/Testdata.  

  
# News
<a name="news"/>

  
# sc-analysis Workflow summary
<a name="workflow_summary"/>

The workflow comprises different modules that can be run sequentially or independently as long as the required data input and object structure is provided. The modules are categorized into 'Pre-processing core modules', 'Downstream analysis', and 'Supporting functionalities' modules.  

## Pre-processing core modules
Core modules perform substantial scRNA seq data pre-processing steps, allowing quality estimationn and guided desicion making for algorithm selection and parameter setting in an iterative process. Hence, although the modules can be run independent of each other, a subsequent conduction of the core modules is recommended to acertain appropriate quality assesment and pre-processing performance.

### Module: qc
Core module to estimate cell QC and filter parameter, investigate covariants, evaluate batch effects, and define normalisation, scaling, and sample combination strategy as well as number of principle components to use.
* Read data
   * Read gene annotation
   * Read scRNA-seq data
* Quality control
   * Determining filter thresholds
   * Genes with highest expression
* On downsampled data, normalization, combining multiple samples, and  dimensionality reduction
   * Evaluation of normalisation and intigration method 
   * Investigating covariants
   * Evaluation of batch effects
   * Determinig dimensionality of the dataset (number of PC to be used; Ellbowplot)

### Module: pre-processing
Core module to perform filtering, normalization, scaling, and sample combination as well as clustering and to determine suitable cluster resolution.
* Read data
   * Read gene annotation
   * Read scRNA-seq data
* Pre-processing
   * Filtering
   * Quality control post-filtering
   * Normalization, scaling, and combining multiple samples
   * Dimensionality reduction
* Clustering
   * Evaluation of resolutions
* Export data


### Module: cluster_analysis
Core module to evaluate and analyse cell clusters. 
* Visualisation with UMAP
* Distribution of cells in clusters
* Cell Cycle Effect
* Cluster QC
* Identification of marker genes



## Downstream analysis

### Module: dataset_mapping
Single cell transcriptomes can be difficult to annotate without extensive knowledge of the underlying biology. Hence, the biological knowledge (defined marker genes and cluster identities) can be propagated from a previously annotated dataset to the test dataset in an automated manner and aid in cluster identification. 
This module maps the cluster annotations from a reference dataset onto the query dataset. Reference and query dataset both need to be provided as Seurat objects.

### Module: ccc_analysis
Cell-cell communication (CCC) is a process by which cells react to stimuli during many biological processes. This module utilizes the LIANA tool to infer ligand-receptor interactions between cell types by running multiple CCC inference methods using a consensus resource and combines the results.

### Module: compositional_analysis


### Module: trajectory_analysis



## Supporting functionalities

### Module: download_references
Module to download reference genome from ENSMBL via BioMart data mining tool.

### Module: download_test_datasets
Module to download test datasets. Test datasets are automatically stored in the appropriate format within the data folder. 
* download_10x_pbmc_1k_healthyDonor_v3Chemistry
* download_10x_pbmc_5k_protein
* download_10x_pbmc_hto_GSE108313
* download_10x_pbmc_small_split2samples
* download_10x_SmartSeq2_pbmc_GSE132044

### Module: read_data
Module executing diverse steps to read data for subsequent processing and analysis.
* read_gene_annotation
* read_rds
* loupe_clustering





# Installation
<a name="installation"/>
  
# Quick start
<a name="quick_start"/>


The workflow is inialised for test dataset '10x_pbmc_small_split2samples'.  

The repository provides several other useful test dataset that you can use to get to know the functionality of the workflow. To run the workflow for another than the initial dataset, you need to select the respective data under the 'Parameter' 'Basic settings'. 

# Usage
<a name="usage"/>


# Output
<a name="output"/>
The core modules as well as the modules for downstream analysis generate RMarkdown reports in html format with comprehensive visualisations, tables and documentation describing the analysis steps.  
The output is created within subfolders carrying the names of the respective modules in a project folder under the general output folder

The output varies depending on the module. However, the most modules generate the following output data:
* RMarkdown report (html format)
* Output files   
   * Export to Loupe Cell Browser
   * Export to cellxgene browser
   * Other output files (e.g. respective tables or plots) 
* Parameter table
* Software versions
* Credits and References

  
# Documentation 
<a name="documentation"/>

Comprehensive documentation can be found in the `docs/` directory:
 
[Installation](docs/...)   
[Running the workflow](docs/...)   

  
# Credits
<a name="credits"/>

The workflow is based on the scrnaseq workflow (status May 4, 2022; Code on Zenodo: https://zenodo.org/record/7849063) developed by [Katrin Sameith](https://github.com/ktrns) and [Andreas Petzold](https://github.com/andpet0101) at the [Dresden-concept Genome Center (Dresden, Germany)](https://genomecenter.tu-dresden.de/about-us).
The workflow is based on the [Seurat](https://satijalab.org/seurat/) package and the vignettes were used as templates. 
Many thanks to all who have contributed.

  
# Contributions and Support
<a name="contributions_and_support"/>

  
# Citation
<a name="citation"/>

If you use this workflow to analyse your data, please cite it by mentioning the Research Core Unit Genomics (RCUG) of Hannover Medical School "https://www.mhh.de/en/genomics" and the Dresden-concept Genome Center "https://genomecenter.tu-dresden.de". 
