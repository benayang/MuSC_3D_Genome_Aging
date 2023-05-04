library(ArchR)
library(dplyr)
library(tidyr)
library(parallel)

#ArchR::installExtraPackages()
projdir = '/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis'

addArchRGenome("mm10")
addArchRThreads(threads = 45) 

projAging2 = readRDS(file.path(projdir, "Save-projAging2-01", "Save-ArchR-Project.rds"))

# Add LSI and first round of UMAP embedding ------------------------------------------------------

projAging3 <- addIterativeLSI(
    ArchRProj = projAging2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 3, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.05, 0.1),
        sampleCells = 10000, 
        n.start = 10,
        algorithm = 2
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    scaleDims = F
)

projAging3 <- addUMAP(
    ArchRProj = projAging3, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

# Harmony for stronger batch correction ------------------------------------------------------

projAging3 <- addHarmony(
    ArchRProj = projAging3,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    plot_convergence = T,
    max.iter.cluster = 400,
    force = T
)

projAging3 <- addUMAP(
    ArchRProj = projAging3, 
    reducedDims = "Harmony", 
    name = "UMAP_Harmony", 
    nNeighbors = 30, 
    minDist = 0.3, 
    metric = "cosine",
    force = T
)

p1 = plotEmbedding(projAging3, "UMAP")
p2 = plotEmbedding(projAging3, "UMAP_Harmony")
ggAlignPlots(p1, p2, type='h')

# MAGIC imputation for gene score plotting ------------------------------------------------------

projAging3 <- addImputeWeights(projAging3)

markerGenes  <- c("Fabp4","Cdh5","Pecam1",
                   "Gsn","Col3a1","Pdgfra",
                   "Thbs4","Fmod","Tnmd",
                   "Pax7","Sdc4","Myod1",
                   "Myh11","Synpo2","Acta2",
                   "Rgs5","Pdgfrb","Kcnj8",
                   "Cd74","H2-Aa","Bach2",
                   "Ptpn22","Il7r","Ptprc",
                   "Lyz2","Ctsb","Ctsz",
                   "S100a9","S100a8","Mmp9",
                   "Cdh19","Ptn","Postn",
                   "Bnc2","Apod","Itgb4",
                   "Mpz","Plp1","Mbp")

p = plotEmbedding(
    ArchRProj = projAging3, 
    colorBy = "GeneScoreMatrix", 
    name = "Pax7", 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(projAging3)
)

p = plotEmbedding(
    ArchRProj = projAging3, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(projAging3)
)

pc <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})

plotPDF(plotList = pc, 
    name = "Plot-UMAP-Marker-Genes-W-Imputation-GeneScore.pdf", 
    ArchRProj = projAging3, 
    addDOC = FALSE, width = 4, height = 4)

# Add Seurat clusters ------------------------------------------------------

projAging3 <- addClusters(
    input = projAging3,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Harmony_Clusters",
    resolution = 0.4
)

p3 = plotEmbedding(projAging3, embedding="UMAP_Harmony", groupBy="cellColData", name="Harmony_Clusters")

cM <- confusionMatrix(paste0(projAging3$Harmony_Clusters), paste0(projAging3$Sample))

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
p

# Look at marker genes per cluster ------------------------------------------------------

markersGS <- getMarkerFeatures(
    ArchRProj = projAging3, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Harmony_Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

marker_hmp = ComplexHeatmap::draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(marker_hmp, name = "marker_peak_heatmap.pdf", ArchRProj = projAging3, addDOC = FALSE, width = 8, height = 4) 

pma <- plotMarkers(seMarker = markersGS, name = "C12", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")

plotPDF(pma, name = "marker_peak_diff_MuSC_Cluster12.pdf", ArchRProj = projAging3, addDOC = FALSE, width = 4, height = 4) 

plotPDF(p1, p2, p3, p$Pax7, p$Tnmd, p$Pdgfra, p$Pdgfrb, p$Acta2, p$Mmp9, p$Cdh5, p$Pecam1, p$Ptn, p$Lyz2, p$Mpz,
         name = "plots.pdf", ArchRProj = projAging3, addDOC = FALSE, width = 4, height = 4) 

# p2c <- lapply(p2, function(x){
#     x + guides(color = FALSE, fill = FALSE) + 
#     theme_ArchR(baseSize = 6.5) +
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
#     theme(
#         axis.text.x=element_blank(), 
#         axis.ticks.x=element_blank(), 
#         axis.text.y=element_blank(), 
#         axis.ticks.y=element_blank()
#     )
# })

saveArchRProject(ArchRProj = projAging3, outputDirectory = "Save-projAging3-01", load = FALSE)