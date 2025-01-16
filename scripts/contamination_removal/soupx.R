### SoupX - mRNA contamination removal
################################################################################

### Create output directories
param$path_out = file.path(param$cellranger_dir, "SoupX")
if (!file.exists(param$path_out)) dir.create(param$path_out, recursive=TRUE, showWarnings=FALSE)


### Required libraries
library(SoupX)
library(Seurat)



### Load data
# Load data and estimate soup profile
# Top level cellranger output directory (the directory that contains the "filtered_gene_bc_matrices" and "raw_gene_bc_matrices" folder).
sc = load10X(param$cellranger_dir)
# Estimate rho
sc = autoEstCont(sc)
# Clean the data
out = adjustCounts(sc)



### Save output 
# Retrieve parameter
rhoProbes=seq(0,1,.001)
priorRho = sc$fit$priorRho
priorRhoStdDev = sc$fit$priorRhoStdDev
rhoEst = sc$fit$rhoEst
rhoFWHM = sc$fit$rhoFWHM
tmp = sc$fit$posterior

v2 = (priorRhoStdDev/priorRho)**2
k = 1 +v2**-2/2*(1+sqrt(1+4*v2))
theta = priorRho/(k-1)
xx=dgamma(rhoProbes,k,scale=theta)

# Save plot
tiff(filename = file.path(param$path_out,'rhoEst.tiff'), width = 2400, height = 1200, res = 300)
plot(rhoProbes,tmp,'l',
     xlim=c(0,1),
     ylim=c(0,max(c(xx,tmp))),
     frame.plot=FALSE,
     xlab='Contamination Fraction',
     ylab='Probability Density')
# Add prior
lines(rhoProbes,xx,lty=2)
abline(v=rhoProbes[which.max(tmp)],col='red')
legend(x='topright',
       legend=c(sprintf('prior rho %g(+/-%g)',priorRho,priorRhoStdDev),
                sprintf('post rho %g(%g,%g)',rhoEst,rhoFWHM[1],rhoFWHM[2]),
                'rho max'),
       lty=c(2,1,1),
       col=c('black','black','red'),
       bty='n')
dev.off()

# Print estimate
cat(sprintf("Estimated global rho of %.2f",rhoEst), file = file.path(param$path_out,'rhoEst.txt'), sep = "\n", append = TRUE )



### Estimation of Soup Correction 
# Get top markers
mrks = quickMarkers(sc$toc,sc$metaData$clusters,N=Inf)
# And only the most specific entry for each gene
mrks = mrks[order(mrks$gene,-mrks$tfidf),]
mrks = mrks[!duplicated(mrks$gene),]
# Order by tfidif maxness
mrks = mrks[order(-mrks$tfidf),]
# Apply tf-idf cut-off
tfidfMin=1
mrks = mrks[mrks$tfidf > tfidfMin,]
# Of the top markers, order by soup 
mrks = mrks[order(sc$soupProfile[mrks$gene,'est'],decreasing=TRUE),]
# And keep the top 20
nonExpressedGene = mrks$gene[seq(min(nrow(mrks),5))]
cat(sprintf("Found %d marker genes",nrow(mrks)), file = file.path(param$path_out,'rhoEst.txt'), sep = "\n", append = TRUE )

# Plot MarkerDistribution of 20 top Markers
tiff(filename = file.path(param$path_out,'MarkerDistribution.tiff'), width = 2400, height = 1200, res = 300)
plotMarkerDistribution(sc)
dev.off()

# Plot MarkerMap und ChangeMap for top 5 Markers
p_list = list()
for (genes in nonExpressedGene) {
  p1 = plotMarkerMap(sc, geneSet = genes) +
    ggplot2::ggtitle(genes)
  p2 = plotChangeMap(sc, out, geneSet = genes)
  p_list[[genes]] = patchwork::wrap_plots(p1+p2)
}
p = patchwork::wrap_plots(p_list, ncol=1)

tiff(filename = file.path(param$path_out,'MarkerMap.tiff'), width = 2400, height = 4500, res = 300)
print(p)
dev.off()



### Save corrected matrix
DropletUtils::write10xCounts(file.path(param$path_out,'strainedCounts'), out)


### Modify output for sc-analysis workflow
features_ids_types = read.delim(file.path(param$path_out,'strainedCounts', "genes.tsv"), header=FALSE, stringsAsFactors=FALSE)
features_ids_types$V3 = "Gene Expression"
write.table(features_ids_types, file=file.path(param$path_out,'strainedCounts', "features.tsv"), sep='\t', row.names = FALSE, col.names = FALSE)



