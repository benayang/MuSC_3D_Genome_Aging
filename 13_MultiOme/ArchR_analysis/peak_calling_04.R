library(ArchR)
library(dplyr)
library(tidyr)
library(parallel)
library(RColorBrewer)
library(scales)

#ArchR::installExtraPackages()
projdir = '/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis'

addArchRGenome("mm10")
addArchRThreads(threads = 45) 

projAging4 = readRDS(file.path(projdir, "Save-projAging4-01", "Save-ArchR-Project.rds"))

# Get marker genes again ------------------------------------------------------

markersGS <- getMarkerFeatures(
    ArchRProj = projAging4, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Harmony_Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

immune_markersGS <- getMarkerFeatures(
    ArchRProj = projAging4, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Harmony_Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    useGroups = c("C3","C9")
)
immune_markerList <- getMarkers(immune_markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

projAging4 <- addClusters(
    input = projAging4,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Harmony_Clusters_v2",
    resolution = 0.8,
    force = T
)

# Group clusters into cell type ------------------------------------------------------

# C1, C2, C10 - Endothelial 1
# C3, C4 - Smooth Muscle
# C5 - MuSC
# C6 - Tenocyte
# C7, C8, C11 - MSC
# C9 - Schwann
# C12 - Neutrophil
# C13 - Neural
# C14 - Immune

projAging4$celltype = factor(sapply(projAging4$Harmony_Clusters, function(x) { 
                        switch(x,
                            "C1" = "Endothelial",
                            "C2" = "Endothelial",
                            "C8" = "Endothelial",
                            "C4" = "Smooth Muscle",
                            "C5" = "Smooth Muscle",
                            "C6" = "MSC",
                            "C13" = "MSC",
                            "C12" = "MuSC",
                            "C7" = "Schwann Cell",
                            "C10" = "Schwann Cell",
                            "C11" = "Schwann Cell",
                            "C14" = "Tenocyte",
                            "C9" = "B/T Cells",
                            "C3" = "Macrophage") }),
                    levels = c("Endothelial","Smooth Muscle","MuSC","Tenocyte","MSC","Schwann Cell","B/T Cells","Macrophage"))


# Modified from https://github.com/NoemieL/bubble-plot-ArchR/blob/main/Script

# Bubble plot for ArchR data script ------------------------------------------------------

genescores = getMatrixFromProject(projAging4, useMatrix="GeneScoreMatrix")
genescoredata = assay(genescores)
genedata = as.data.frame(rowData(genescores))
celldata = as.data.frame(colData(genescores))
numcells = genescoredata@Dim[2]

geneintegrationscores = getMatrixFromProject(projAging4, useMatrix="GeneIntegrationMatrix")
geneintegrationscoredata = assay(geneintegrationscores)
geneintdata = as.data.frame(rowData(geneintegrationscores))
cellintdata = as.data.frame(colData(geneintegrationscores))
numintcells = geneintegrationscoredata@Dim[2]

custom_DotPlot = function(gene_list, genedata, celldata, genescoredata, clusters, cluster_levels, legend_label) {
    cluster_size = as.data.frame(table(celldata[[clusters]]))
    colnames(cluster_size) = c("Cluster","Cluster_size")
    
    # Get row indices of gene score activity matrix from gene list
    idx = as.numeric(rownames(genedata[genedata$name %in% gene_list, ]))

    genescoredataframe = list()
    for(i in unique(celldata[[clusters]])) {
        tmp = data.frame(avg_exp = Matrix::rowMeans(genescoredata[idx, rownames(celldata)[celldata[[clusters]]==i]]),
                        pct_exp = Matrix::rowSums(genescoredata[idx, rownames(celldata)[celldata[[clusters]]==i]]>0) / cluster_size$Cluster_size[cluster_size$Cluster == i] * 100,
                        gene = genedata[genedata$name %in% gene_list, "name"],
                        Cluster = i)
        genescoredataframe[[i]] = tmp
    }
    genescoredataframe = bind_rows(genescoredataframe, .id="Cluster")
    genescoredataframe = genescoredataframe %>% mutate(gene = factor(gene, levels=gene_list), Cluster = factor(Cluster, levels=cluster_levels), log1p_exp = log(avg_exp + 1))

    plt = ggplot(data = genescoredataframe, mapping = aes(x = gene, y = Cluster)) +
    geom_point(mapping = aes(size = pct_exp, color = log1p_exp)) +
    scale_colour_gradient(limits = c(0, max(genescoredataframe$log1p_exp)), low="lightgrey", high="blue") +
    scale_size(range=c(0,6), limits = c(0, max(genescoredataframe$pct_exp))) +
    cowplot::theme_cowplot() +
    labs(
        x = NULL,
        y = NULL
    ) +
    guides(color = guide_colorbar(title=legend_label), size = guide_legend(title="Percent\nExpressed")) +
    theme(axis.text.x = element_text(size=14, angle=35, hjust=1),
            axis.text.y = element_text(size=14))

    return(plt)
}

gene_list  <- c("Fabp4","Cdh5","Flt1",
                "Il17b","Pdgfrb","Accsl",
                "Pax7","Sdc4","Runx1",
                "Loxl4","Prelp","Thbs4",
                "Gfpt2","Entpd2","Smoc2",
                "Nkain2","Ano3","Ryr3",
                "Dusp2","H2-Aa","Itk",
                "Cd300a","Ptafr","Plek")

plt = custom_DotPlot(gene_list, genedata, celldata, genescoredata, "celltype", levels(projAging4@cellColData$celltype), "Average\nGene\nActivity") + theme(text = element_text(family="Arial"))
ggsave(plt, filename=file.path(getOutputDirectory(projAging4), "Plots", "Bubble-plot-ATAC-celltype.png"), width=10, height=4.5, dpi=300)

plt2 = custom_DotPlot(gene_list, geneintdata, cellintdata, geneintegrationscoredata, "celltype", levels(projAging4@cellColData$celltype), "Average\nIntegrated\nExpression") + theme(text = element_text(family="Arial"))
ggsave(plt2, filename=file.path(getOutputDirectory(projAging4), "Plots", "Bubble-plot-integrated-RNA-celltype.png"), width=10, height=4.5, dpi=300)

# young_cells = projAging4$cellNames[projAging4$Age=="Young"]
# aged_cells = projAging4$cellNames[projAging4$Age=="Aged"]

# p1 = plotEmbedding(projAging4, embedding="UMAP_Harmony", colorBy="cellColData", name="celltype")
# p2 = plotEmbedding(projAging4, embedding="UMAP_Harmony", colorBy="cellColData", name="Harmony_Clusters_v2")
# p3 = plotEmbedding(projAging4, embedding="UMAP_Harmony", colorBy="cellColData", name="Age", highlightCells=aged_cells, pal = c(Aged = ggsci::pal_nejm()(2)[1]))
# p4 = plotEmbedding(projAging4, embedding="UMAP_Harmony", colorBy="cellColData", name="Age", highlightCells=young_cells, pal = c(Young = ggsci::pal_nejm()(2)[2]))

# projAging4 = addImputeWeights(projAging4)
# p5 = plotEmbedding(projAging4, embedding="UMAP_Harmony", colorBy="GeneScoreMatrix", name="Pax7", imputeWeights = getImputeWeights(projAging4))
# p6 = plotEmbedding(projAging4, embedding="UMAP_Harmony", colorBy="GeneScoreMatrix", name="Myod1", imputeWeights = getImputeWeights(projAging4))
# plotPDF(p1,p2,p3,p4,p5,p6, name = "clustering_plots.pdf", ArchRProj = projAging4, addDOC = FALSE, width = 5, height = 5)


# Add pseudo-bulk per cell type ------------------------------------------------------

library(BSgenome.Mmusculus.UCSC.mm10)

projAging4 <- addGroupCoverages(ArchRProj = projAging4, groupBy = "celltype", minCells = 50, maxCells=500, minReplicates = 3, maxReplicates=6)

pathToMacs2 <- findMacs2()

projAging5 <- addReproduciblePeakSet(
    ArchRProj = projAging4, 
    groupBy = "celltype", 
    pathToMacs2 = pathToMacs2
)

projAging5 = addPeakMatrix(projAging5)

saveArchRProject(ArchRProj = projAging5, outputDirectory = "Save-projAging5-01", load = FALSE)

# Look at cell type-specific peaks  ------------------------------------------------------

# cluster IDs need to be characters not factors (https://github.com/GreenleafLab/ArchR/issues/238)
projAging5$celltype = as.character(projAging5$celltype)

markersPeaks <- getMarkerFeatures(
    ArchRProj = projAging5, 
    useMatrix = "PeakMatrix", 
    groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  binarize =  T
)

markersPeaks_nb <- getMarkerFeatures(
    ArchRProj = projAging5, 
    useMatrix = "PeakMatrix", 
    groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  binarize =  F
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = T)
markerList_nb <- getMarkers(markersPeaks_nb, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = T)

saveRDS(markersPeaks, file.path(getOutputDirectory(projAging5), "marker_peaks_all_list.RDS"))
saveRDS(markerList, file.path(getOutputDirectory(projAging5), "marker_peaks_diff_list.RDS"))

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 1",
  transpose = TRUE
)

heatmapPeaks_nb <- plotMarkerHeatmap(
  seMarker = markersPeaks_nb, 
  cutOff = "FDR <= 0.1 & Log2FC >= 1",
  transpose = TRUE
)

plotPDF(draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot"),
        name = "Plot-Marker-Peaks-Heatmap.pdf", ArchRProj = projAging5, addDOC = FALSE, width = 5, height = 5)

# Get QC metrics  ------------------------------------------------------

cellcoldata = getCellColData(projAging4)
cellcoldata %>% as.data.frame() %>% group_by(Age) %>% summarise(median(nFrags))
cellcoldata %>% as.data.frame() %>% group_by(Age) %>% summarise(median(TSSEnrichment))
cellcoldata %>% as.data.frame() %>% group_by(Age) %>% summarise(median(BlacklistRatio))

peakmatrix = getMatrixFromProject(projAging5, "PeakMatrix")
young_cellnames = rownames(colData(peakmatrix)[colData(peakmatrix)$Sample %in% c("Young","Young_v2")])
aged_cellnames = rownames(colData(peakmatrix)[colData(peakmatrix)$Sample=="Aged",])

summary(Matrix::colSums(assay(peakmatrix)[,young_cellnames] > 0))
summary(Matrix::colSums(assay(peakmatrix)[,aged_cellnames] > 0))

