### Configuration
################################################################################
# Location of rds object
DataPath = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/output/2024_040277_3//cluster_analysis/data/sc.rds"
# Output path (only needed for saving plots)
OutputPath = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/output/2024_040277_3/additional_visualisation/"
# Location of inspect_rds.R script
# Needed packages: Seurat, ggplot2, magrittr
ScriptPath = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/docs/handout_scripts/"



### Initialize additional visualization
################################################################################
source(file.path(ScriptPath,"R/initialize_additional_visualization.R"))



### Generate and save plots
################################################################################

### Generate Standard Plots
## Arguments
# sc: Seurat object. Output from sc_analysis workflow. DO NOT CHANGE!!!
# param: Parameter list. Output from sc_analysis workflow. DO NOT CHANGE!!!
# ncol: The number of plots per row.
# plottheme: ggplot theme, e.g. theme_bw(), theme_classic(), theme_light(), theme_void(), theme_minimal(). Default: theme_light() 
# clustercolors: Color of clusters.
# samplecolors: Color of samples.
# celltypecolors: Color of cell type annotation. 

## Output: 
# 1: "clusters",  2: "clusters_separately",
# 3: "sample",   4: "sample_separately",
# 5: "cluster_annotation",   6: "cluster_annotation_separately",
# 7: "cell_annotation",  8: "cell_annotation_separately",
# 9: "nCount_RNA",   10: "nFeature_RNA",   11: "percent_mt",  12: "percent_ribo"    


p_list = ADStandardPlots(sc = sc, param = param)


# Print plot
# All plots are saved in p_list[[]]. You can retrieve the respective plot by inserting the name, e.g. p_list[["clusters"]]
# Putting the cursor in the middle the brackets and pressing TAB will display all possibilities
p_list[[]]

# Save plot
# As above, you have to select a individual plot from the list, e.g. p_list[["clusters"]]
# If a file with the same name is existing in the folder already, it will be overwritten
# Change the file name as fitting, but do not change path or ".tiff"
# You can change the width and height as it will fit to your plot
tiff(filename = paste0(OutputPath, "plotname.tiff"), width = 2400, height = 1200, res = 300)
p_list[[]]
dev.off()




### Generate Feature Plots
## Arguments
# sc: Seurat object. Output from sc_analysis workflow. DO NOT CHANGE!!!
# param: Parameter list. Output from sc_analysis workflow. DO NOT CHANGE!!!
# markers: Vector with gene names using gene symbol, e.g. c("Kit", "Metrnl")
# ncol: The number of plots per row.
# featurecolor: Color (scaled) for features.
# bgcolor: Background color for cells without feature expression (lower end of scale)
# labelcluster: Define whether to show cluster labels. Default: TRUE.

genes = c("Kit", "Metrnl")
p = ADFeaturePlots(sc = sc, param = param, markers = genes)
p 


### Generate Violin Plots
## Arguments
# sc: Seurat object. Output from sc_analysis workflow. DO NOT CHANGE!!!
# param: Parameter list. Output from sc_analysis workflow. DO NOT CHANGE!!!
# markers: Vector with gene names using gene symbol, e.g. c("Kit", "Metrnl")
# group: Define how to group the plot, e.g. by "clusters", "samples", or "celltypes". Default: "clusters".
# ncol: The number of plots per row.
# plottheme: ggplot theme, e.g. theme_bw(), theme_classic(), theme_light(), theme_void(), theme_minimal(). Default: theme_light() 
# clustercolors: Color of clusters.

genes = c("Kit", "Metrnl")
p = ADViolinPlots(sc = sc, param = param, markers = genes)
p




### Save plots
# If a file with the same name is existing in the folder already, it will be overwritten
# Change the file name as fitting, but do not change path or ".tiff"
# You can change the width and height as it will fit to your plot
tiff(filename = paste0(OutputPath, "plotname.tiff"), width = 2400, height = 1200, res = 300)
p
dev.off()
