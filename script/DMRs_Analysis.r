library(methylKit)
library(GenomicRanges)
library(ggplot2)
library(ComplexHeatmap)
library(annotatr)
library(dplyr)

#import the data
file.list <- list(
  "Control_1_CpG_report.txt",
  "Control_2_CpG_report.txt",
  "Treated_1_CpG_report.txt",
  "Treated_2_CpG_report.txt"
)

sample.id <- c("WT_1","WT_2","Treated_1","Treated_2")
treatment <- c(0,0,1,1)    # 0 = control, 1 = treated

myobj <- methRead(
  file.list,
  sample.id = as.list(sample.id),  # Convert to explicit list
  treatment = treatment,
  assembly = "hg38",
  context = "CpG",
  pipeline = "bismarkCytosineReport",
  header = FALSE,
  mincov = 5
)

# quality filteration
filtered.obj <- filterByCoverage(
  myobj,
  lo.count = 10, lo.perc = NULL,
  hi.count = 5000, hi.perc = 99.9
)

# Normalization
norm.obj <- normalizeCoverage(filtered.obj)

# Unite CpG Sites Across Samples
meth <- unite(norm.obj, destrand = TRUE)
head(meth)

#PCA 
PCASamples(meth, screeplot = TRUE)
PCASamples(meth)

# Clustering
clusterSamples(meth)

#Methylation Percentage Distribution
# If you still have your original methylRawList object
getMethylationStats(myobj[[1]], plot = TRUE)  # For first sample
getMethylationStats(myobj[[2]], plot = TRUE)  # For second sample, etc.

# Or apply to all samples
for(i in 1:length(myobj)) {
  getMethylationStats(myobj[[i]], plot = TRUE)
}

# Get percent methylation matrix
perc.meth <- getData(meth)[, grep("numCs", names(getData(meth)))] / 
  getData(meth)[, grep("coverage", names(getData(meth)))] * 100

# Calculate basic stats
summary(perc.meth)

# Correlation between samples
getCorrelation(meth, plot = TRUE)

# PCA
pc <- prcomp(t(perc.meth))
plot(pc$x[,1], pc$x[,2], xlab="PC1", ylab="PC2")

# Differential methylation
myDiff <- calculateDiffMeth(meth)

# Differential Methylation Analysis (DMRs at CpG level)
myDiff <- calculateDiffMeth(meth)

#Extract significant CpGs:
sigCpG <- getMethylDiff(myDiff, qvalue = 0.01, difference = 25)  # 25% difference

# Identify Differentially Methylated Regions (DMRs)
DMRs <- tileMethylCounts(norm.obj, win.size = 1000, step.size = 1000)


# Unite the tiled data into methylBase
DMRs_united <- unite(DMRs, destrand = FALSE)

# Calculate differential methylation on the united object
DMRs.diff <- calculateDiffMeth(DMRs_united)

#  Get significant DMRs
sigDMRs <- getMethylDiff(DMRs.diff, qvalue = 0.01, difference = 10)

# check results
head(sigDMRs)
nrow(sigDMRs) 

# Annotation to Genes / Promoters
annotations <- build_annotations(genome = "hg38", 
                                 annotations = c("hg38_basicgenes",
                                                 "hg38_genes_promoters"))

annotated <- annotate_regions(
  regions = as(sigDMRs, "GRanges"),
  annotations = annotations,
  ignore.strand = TRUE
)

head(annotated)

#Visualizations

#Methylation Heatmap
meth.mat <- percMethylation(meth)
Heatmap(meth.mat, name = "Methylation (%)")

# Volcano plot of differentially methylated CpGs
df <- as.data.frame(myDiff)

ggplot(df, aes(x = meth.diff, y = -log10(qvalue))) +
  geom_point(alpha = 0.5) +
  theme_minimal()

#DMRs distribution across genome
plotChromosomes(DMRs.diff)

# Gene Set Enrichment Analysis (GSEA)
geneList <- unique(annotated$gene_id)

# Example: KEGG enrichment
library(clusterProfiler)
enrich <- enrichKEGG(geneList, organism = "hsa")
dotplot(enrich)

#Export Results
write.table(sigCpG, "significant_CpGs.txt", sep="\t", row.names=FALSE)
write.table(sigDMRs, "significant_DMRs.txt", sep="\t", row.names=FALSE)
write.table(annotated, "annotated_DMRs.txt", sep="\t", row.names=FALSE)
