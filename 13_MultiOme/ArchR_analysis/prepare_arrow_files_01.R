library(ArchR)
library(dplyr)
library(tidyr)
library(parallel)

#ArchR::installExtraPackages()
projdir = '/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis'

addArchRGenome("mm10")
addArchRThreads(threads = 45) 

inputFiles = c("/nas/homes/benyang/HiC/13_MultiOme/RawData/Aged/10x_analysis_5422-JL/Sample_5422-JL-1/atac_fragments.tsv.gz",
"/nas/homes/benyang/HiC/13_MultiOme/RawData/Aged_v2/atac_fragments.tsv.gz",
"/nas/homes/benyang/HiC/13_MultiOme/RawData/Young/10x_analysis_5644-JL/Sample_5644-JL-1/atac_fragments.tsv.gz",
"/nas/homes/benyang/HiC/13_MultiOme/RawData/6098-JL/10x_analysis_6098-JL/Sample_6098-JL-1/atac_fragments.tsv.gz")
names(inputFiles) = c("Aged","Aged_v2","Young","Young_v2")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1e3,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

projAging1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/output",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

saveArchRProject(ArchRProj = projAging1, outputDirectory = file.path(projdir, "Save-projAging1-01"), load = FALSE)

## What should the TSSEnrichment cutoff be?
df = getCellColData(projAging1, select = c("Sample", "log10(nFrags)", "TSSEnrichment"))

png(file.path("/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis", "QualityControl", "TSS_enrichment_score_density.png"), res=300, units='in', width=5, height=5)
plot(density(df[df$Sample=="Aged", "TSSEnrichment"]), col='red', ylim = c(0, 0.1), main="TSS Enrichment Score Density")
lines(density(df[df$Sample=="Young", "TSSEnrichment"]), col='blue')
lines(density(df[df$Sample=="Young_v2", "TSSEnrichment"]), col='green')
legend(x='topright', legend=c("Aged","Young","Young_v2"), col=c('red','blue','green'), lty=1)
abline(v=6)
dev.off()

projAging2 = projAging1[projAging1$TSSEnrichment>=6 & projAging1$BlacklistRatio<0.05, ]

projAging2 = filterDoublets(projAging2)
# Filtering 880 cells from ArchRProject!
#         Aged : 229 of 4787 (4.8%)
#         Young_v2 : 395 of 6292 (6.3%)
#         Young : 256 of 5064 (5.1%)

projAging1$Age = sapply(projAging1$Sample, function(x) ifelse(x=="Aged", "Aged", "Young"))
projAging2$Age = sapply(projAging2$Sample, function(x) ifelse(x=="Aged", "Aged", "Young"))

p1 = plotFragmentSizes(projAging2, groupBy="Age", pal=ggsci::pal_nejm()(2))
p2 = plotTSSEnrichment(projAging2, groupBy="Age", pal=ggsci::pal_nejm()(2))

df2 = getCellColData(projAging1, select = c("Age", "log10(nFrags)", "TSSEnrichment"))
plot_qc = function(df, age) {
  p = ggPoint(
    x = df[df$Age==age, 2], 
    y = df[df$Age==age, 3], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment"
  ) + geom_hline(yintercept = 6, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed") +
  guides(color = guide_colorbar(title="Density")) +
  theme(legend.position = "top", legend.text = element_text(size=12), legend.title = element_text(size=12),
        legend.key.width = unit(1, "cm"),
        axis.text = element_text(size=10), axis.title = element_text(size=12))

  return(p)
}

p3 = plot_qc(df2, "Aged")
p4 = plot_qc(df2, "Young")

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projAging2, addDOC = FALSE, width = 4, height = 4)

ggsave(p1, filename=file.path(getOutputDirectory(projAging2), "Plots", "Fragment_sizes_qc.png"), dpi=300, width=4, height=4)
ggsave(p2, filename=file.path(getOutputDirectory(projAging2), "Plots", "TSS_enrichment_qc.png"), dpi=300, width=4, height=4)
ggsave(p3, filename=file.path(getOutputDirectory(projAging2), "Plots", "young_sample_qc.png"), dpi=300, width=5, height=5)
ggsave(p4, filename=file.path(getOutputDirectory(projAging2), "Plots", "aged_sample_qc.png"), dpi=300, width=5, height=5)

saveArchRProject(ArchRProj = projAging2, outputDirectory = "Save-projAging2-01", load = FALSE)