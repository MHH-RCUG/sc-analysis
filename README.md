# Preamble  
  
**The workflow is set up for the rcug user on hpc-rc11!**  
**Before starting, see docs/Run_scrnaseq_workflow.docx for an instruction how to run the scrnaseq workflow on our system!**    
- includes renv, python env, and instruction for usage  
- includes installation documentation and instruction for usage on our system   
- includes documentation of package versions  
  
  
----------------------------------------------------------
  

# Table of contents
* [Introduction](#introduction) 
* [News](#news) 
* [Workflow summary](#workflow_summary)
* [Quick start](#quick_start)
* [Documentation](#documentation)
* [Credits](#credits)
* [Contributions and support](#contributions_and_support)
* [Citation](#citation)

# Introduction
<a name="introduction"/>

**sc_analysis** is a bioinformatics analysis workflow for single-cell RNA-seq analysis. The workflow is based on Seurat and the scrnaseq workflow (status May 4, 2022; Code on Zenodo: https://zenodo.org/record/7849063) created from Dresden-concept Genome Center URL "https://genomecenter.tu-dresden.de".  
The workflow supports RNA sequencing data from one or more samples processed with 10X Genomics and SmartSeq-2. 
The workflow generates RMarkdown reports in html format with comprehensive visualisations, tables and documentation describing the analysis steps. 

Example reports for the test dataset '10x_pbmc_small_split2samples' can be find the GitHub repository.

  
# News
<a name="news"/>

  
# Workflow summary
<a name="workflow_summary"/>

## Workflow: sc_analysis 

## Module: qc
* Read data
   * Read gene annotation
   * Read scRNA-seq data
* Pre-processing
   * Quality control
   * Genes with highest expression

## Module: pre-processing
* Pre-processing
   * Filtering
   * Quality control post filtering
   * Normalisation
   * Combining multiple samples
   * Relative log expression
   * Dimensionality reduction
   * Dimensionality of the dataset

## Module: clustering
* Clustering

## Module: annotation
* Visualisation with UMAP
* Distribution of cells in clusters
* Cell Cycle Effect
* Cluster QC
* Known marker genes

## Module: download_references

## Module: download_test_datasets

## Module: loupe_clustering

## Module: compositional_analysis
* Differentially expressed genes

## Module: trajectory_analysis

## Output
* RMarkdown report (html format)
* Output files   
   * Export to Loupe Cell Browser
   * Export to cellxgene browser
   * Other output files (e.g. respective tables or plots) 
* Parameter table
* Software versions
* Credits and References

  
# Quick start
<a name="quick_start"/>

The workflow is inialised for test dataset '10x_pbmc_small_split2samples'.  

The repository provides several other useful test dataset that you can use to get to know the functionality of the workflow. To run the workflow for another than the initial dataset, you need to ... 

  
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
