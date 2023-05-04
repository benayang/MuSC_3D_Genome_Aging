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

projAging6 = readRDS(file.path(projdir, "Save-projAging6-01", "Save-ArchR-Project.rds"))
JAL_MuSC = readRDS("/nas/homes/benyang/HiC/13_MultiOme/JAL_RNA/aging_MuSC_integrated.RDS")

# Re-cluster MuSCs by age  ------------------------------------------------------

MuSC_cellnames = projAging6$cellNames[projAging6$celltype=="MuSC"]

MuSC_projAging1 = subsetArchRProject(projAging6, cells = MuSC_cellnames, outputDirectory = file.path(projdir,"MuSC_ArchR","all_MuSC"), force = T)

MuSC_projAging1 <- addIterativeLSI(
    ArchRProj = MuSC_projAging1,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force=T
)

MuSC_projAging1 <- addUMAP(
    ArchRProj = MuSC_projAging1, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = T
)

MuSC_projAging1 <- addHarmony(
    ArchRProj = MuSC_projAging1,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = T
)

MuSC_projAging1 <- addUMAP(
    ArchRProj = MuSC_projAging1, 
    reducedDims = "Harmony", 
    name = "UMAP_Harmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = T
)

p1 = plotEmbedding(ArchRProj = MuSC_projAging1, colorBy = "cellColData", name = "Age", embedding = "UMAP", size = 1, pal = c(Young = ggsci::pal_nejm()(2)[2], Aged = ggsci::pal_nejm()(2)[1]))
p2 = plotEmbedding(ArchRProj = MuSC_projAging1, colorBy = "cellColData", name = "Age", embedding = "UMAP_Harmony", size = 1, pal = c(Young = ggsci::pal_nejm()(2)[2], Aged = ggsci::pal_nejm()(2)[1]))

plotPDF(p1, p2, name = "Plot-MuSC-Reclustered-Age.pdf", ArchRProj = MuSC_projAging1, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = MuSC_projAging1, outputDirectory = "Save-MuSC_projAging1-01", load = FALSE)

# MuSC QC metrics  ------------------------------------------------------

MuSC_TSS = plotTSSEnrichment(MuSC_projAging1, groupBy="Age", pal=ggsci::pal_nejm()(2))

png(file.path(getOutputDirectory(MuSC_projAging1), "Plots", "MuSC_TSS.png"), res=300, width=4, height=4, units='in')
MuSC_TSS + theme(legend.position="top", legend.title = element_blank(), legend.text = element_text(size=12), axis.text = element_text(size=12), axis.title = element_text(size=12))
dev.off()

MuSC_Frag = plotFragmentSizes(MuSC_projAging1, groupBy="Age", pal=ggsci::pal_nejm()(2))

png(file.path(getOutputDirectory(MuSC_projAging1), "Plots", "MuSC_Frag.png"), res=300, width=4, height=4, units='in')
MuSC_Frag + theme(legend.position="top", legend.title = element_blank(), legend.text = element_text(size=12), axis.text = element_text(size=12), axis.title = element_text(size=12))
dev.off()

plotPDF(MuSC_TSS, MuSC_Frag, name = "Plot-MuSC-QC-Age.pdf", ArchRProj = MuSC_projAging1, addDOC = FALSE, width = 5, height = 5)

# Pseudo-bulk and peak calling ------------------------------------------------------

library(BSgenome.Mmusculus.UCSC.mm10)

MuSC_projAging1 = addGroupCoverages(MuSC_projAging1, groupBy="Age", minCells = 50, maxCells=500, minReplicates = 3, maxReplicates=6, force = T)

pathToMacs2 <- findMacs2()

MuSC_projAging1 <- addReproduciblePeakSet(
    ArchRProj = MuSC_projAging1, 
    groupBy = "Age", 
    pathToMacs2 = pathToMacs2
)

MuSC_projAging1 = addPeakMatrix(MuSC_projAging1)

saveArchRProject(ArchRProj = MuSC_projAging1, outputDirectory = "Save-MuSC_projAging1-01", load = FALSE)

# Get MuSC QC metrics  ---------------------------------------------------

cellcoldata = getCellColData(MuSC_projAging1)
cellcoldata %>% as.data.frame() %>% group_by(Age) %>% summarise(median(nFrags))
cellcoldata %>% as.data.frame() %>% group_by(Age) %>% summarise(median(TSSEnrichment))
cellcoldata %>% as.data.frame() %>% group_by(Age) %>% summarise(median(BlacklistRatio))
cellcoldata %>% as.data.frame() %>% group_by(Age) %>% summarise(median(FRIP))

peakmatrix = getMatrixFromProject(MuSC_projAging1, "PeakMatrix")
young_MuSC_cellnames = rownames(colData(peakmatrix)[(colData(peakmatrix)$Sample %in% c("Young","Young_v2")),])
aged_MuSC_cellnames = rownames(colData(peakmatrix)[(colData(peakmatrix)$Sample=="Aged"),])

summary(Matrix::colSums(assay(peakmatrix)[,young_MuSC_cellnames] > 0))
summary(Matrix::colSums(assay(peakmatrix)[,aged_MuSC_cellnames] > 0))

# Isolate young and aged MuSCs  ------------------------------------------------------

aged_MuSC_cellnames = projAging6$cellNames[projAging6$Age=="Aged" & projAging6$celltype=="MuSC"]
young_MuSC_cellnames = projAging6$cellNames[projAging6$Age=="Young" & projAging6$celltype=="MuSC"]

young_MuSC_projAging1 = subsetArchRProject(MuSC_projAging1, cells = young_MuSC_cellnames, outputDirectory = file.path(projdir,"MuSC_ArchR","young_MuSC"), force = T)
aged_MuSC_projAging1 = subsetArchRProject(MuSC_projAging1, cells = aged_MuSC_cellnames, outputDirectory = file.path(projdir,"MuSC_ArchR","aged_MuSC"), force = T)

young_MuSC_projAging1 <- addIterativeLSI(
    ArchRProj = young_MuSC_projAging1,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = T
)

young_MuSC_projAging1 <- addUMAP(
    ArchRProj = young_MuSC_projAging1, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = T
)

aged_MuSC_projAging1 <- addIterativeLSI(
    ArchRProj = aged_MuSC_projAging1,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = T
)

aged_MuSC_projAging1 <- addUMAP(
    ArchRProj = aged_MuSC_projAging1, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = T
)

young_MuSC_projAging1 <- addGeneScoreMatrix(young_MuSC_projAging1, force=T)
aged_MuSC_projAging1 <- addGeneScoreMatrix(aged_MuSC_projAging1, force=T)

young_MuSC_projAging1 <- addImputeWeights(young_MuSC_projAging1)
aged_MuSC_projAging1 <- addImputeWeights(aged_MuSC_projAging1)

young_MuSC_projAging1 = addGroupCoverages(young_MuSC_projAging1, groupBy="Age", minCells = 40, maxCells=500, minReplicates = 3, maxReplicates=6, force = T)
aged_MuSC_projAging1 = addGroupCoverages(aged_MuSC_projAging1, groupBy="Age", minCells = 40, maxCells=500, minReplicates = 3, maxReplicates=6, force = T)

saveArchRProject(ArchRProj = young_MuSC_projAging1, outputDirectory = "Save-young_MuSC_projAging1-01", load = FALSE)
saveArchRProject(ArchRProj = aged_MuSC_projAging1, outputDirectory = "Save-aged_MuSC_projAging1-01", load = FALSE)

# Make cell data set objects from MuSCs ------------------------------------------------------

library(monocle3)
library(cicero)

young_peakmatrix = getMatrixFromProject(young_MuSC_projAging1, "PeakMatrix")
aged_peakmatrix = getMatrixFromProject(aged_MuSC_projAging1, "PeakMatrix")
MuSC_peakmatrix = getMatrixFromProject(MuSC_projAging1, "PeakMatrix")

make_sparse_peakmatrix = function(mat) {
    peakmatrix_sparse = as.data.frame(summary(assay(mat)))
    peakmatrix_sparse = cbind(peakmatrix_sparse, 
                                    as.data.frame(rowRanges(mat))[peakmatrix_sparse$i, ],
                                    cell = rownames(colData(mat))[peakmatrix_sparse$j])
    peakmatrix_sparse = peakmatrix_sparse %>% mutate(range_string = paste(seqnames, start, end, sep="_"))
    return(peakmatrix_sparse)
}

young_peakmatrix_sparse = make_sparse_peakmatrix(young_peakmatrix)
aged_peakmatrix_sparse = make_sparse_peakmatrix(aged_peakmatrix)
MuSC_peakmatrix_sparse = make_sparse_peakmatrix(MuSC_peakmatrix)

young_atac_cds = make_atac_cds(young_peakmatrix_sparse[ ,c("range_string","cell","x")], binarize=T)
aged_atac_cds = make_atac_cds(aged_peakmatrix_sparse[ ,c("range_string","cell","x")], binarize=T)
MuSC_atac_cds = make_atac_cds(MuSC_peakmatrix_sparse[ ,c("range_string","cell","x")], binarize=T)

# young_peakmatrix_data = assay(young_peakmatrix)
# rownames(young_peakmatrix_data) = rowRanges(young_peakmatrix) %>% as.data.frame() %>% mutate(site_name = paste(seqnames, start, end, sep="_")) %>% pull(site_name)
# young_peakmatrix_data@x[young_peakmatrix_data@x > 0] = 1
# young_peakmatrix_site_data = as.data.frame(rowRanges(young_peakmatrix))[,1:3]
# rownames(young_peakmatrix_site_data) = rownames(young_peakmatrix_data)
# colnames(young_peakmatrix_site_data) = c("chr","bp1","bp2")

# young_cds <- new_cell_data_set(young_peakmatrix_data,
#                             cell_metadata = colData(young_peakmatrix),
#                             gene_metadata = young_peakmatrix_site_data)

# aged_peakmatrix_data = assay(aged_peakmatrix)
# rownames(aged_peakmatrix_data) = rowRanges(aged_peakmatrix) %>% as.data.frame() %>% mutate(site_name = paste(seqnames, start, end, sep="_")) %>% pull(site_name)
# aged_peakmatrix_data@x[aged_peakmatrix_data@x > 0] = 1
# aged_peakmatrix_site_data = as.data.frame(rowRanges(aged_peakmatrix))[,1:3]
# rownames(aged_peakmatrix_site_data) = rownames(aged_peakmatrix_data)
# colnames(aged_peakmatrix_site_data) = c("chr","bp1","bp2")

# aged_cds <- new_cell_data_set(aged_peakmatrix_data,
#                             cell_metadata = colData(aged_peakmatrix),
#                             gene_metadata = aged_peakmatrix_site_data)

# Prepare objects for cicero ------------------------------------------------------

set.seed(2017)

prepare_for_cicero = function(cds, proj) {
    cds <- cds[Matrix::rowSums(exprs(cds)) != 0,] 
    cds <- detect_genes(cds)
    cds <- estimate_size_factors(cds)
    reducedDims(cds)[["UMAP"]] <- getEmbedding(proj, "UMAP")
    return(cds)
}

young_atac_cds = prepare_for_cicero(young_atac_cds, young_MuSC_projAging1)
aged_atac_cds = prepare_for_cicero(aged_atac_cds, aged_MuSC_projAging1)
MuSC_atac_cds = prepare_for_cicero(MuSC_atac_cds, MuSC_projAging1)

# Look at effects of varying k parameters ------------------------------------------------------

# use k = 35
young_umap_coords <- reducedDims(young_atac_cds)$UMAP
young_atac_cds.cicero.klist <- lapply(seq(10,60,5), function(x) make_cicero_cds(young_atac_cds, reduced_coordinates = young_umap_coords, 
                                                                                k = x, return_agg_info = T))

young_atac_cds.cicero.klist.df = lapply(young_atac_cds.cicero.klist, function(x) as.data.frame(x[2]))
young_atac_cds.cicero.klist.df = lapply(young_atac_cds.cicero.klist.df, function(x) table(x$cell, x$agg_cell))
young_atac_cds.cicero.klist.df = lapply(young_atac_cds.cicero.klist.df, function(x) rowMeans(x, na.rm=T))
names(young_atac_cds.cicero.klist.df) = seq(10,60,5)
boxplot(young_atac_cds.cicero.klist.df, ylab="Mean Bins Per Cell", xlab="Cells Per Bin")

# use k = 25
aged_umap_coords <- reducedDims(aged_atac_cds)$UMAP
aged_atac_cds.cicero.klist <- lapply(seq(10,60,5), function(x) make_cicero_cds(aged_atac_cds, reduced_coordinates = aged_umap_coords, 
                                                                                k = x, return_agg_info = T))

aged_atac_cds.cicero.klist.df = lapply(aged_atac_cds.cicero.klist, function(x) as.data.frame(x[2]))
aged_atac_cds.cicero.klist.df = lapply(aged_atac_cds.cicero.klist.df, function(x) table(x$cell, x$agg_cell))
aged_atac_cds.cicero.klist.df = lapply(aged_atac_cds.cicero.klist.df, function(x) rowMeans(x, na.rm=T))
names(aged_atac_cds.cicero.klist.df) = seq(10,60,5)
boxplot(aged_atac_cds.cicero.klist.df, ylab="Mean Bins Per Cell", xlab="Cells Per Bin")

# use k = 50
MuSC_umap_coords <- reducedDims(MuSC_atac_cds)$UMAP
MuSC_atac_cds.cicero.klist <- lapply(seq(10,60,5), function(x) make_cicero_cds(MuSC_atac_cds, reduced_coordinates = MuSC_umap_coords, 
                                                                                k = x, return_agg_info = T))

MuSC_atac_cds.cicero.klist.df = lapply(MuSC_atac_cds.cicero.klist, function(x) as.data.frame(x[2]))
MuSC_atac_cds.cicero.klist.df = lapply(MuSC_atac_cds.cicero.klist.df, function(x) table(x$cell, x$agg_cell))
MuSC_atac_cds.cicero.klist.df = lapply(MuSC_atac_cds.cicero.klist.df, function(x) rowMeans(x, na.rm=T))
names(MuSC_atac_cds.cicero.klist.df) = seq(10,60,5)
boxplot(MuSC_atac_cds.cicero.klist.df, ylab="Mean Bins Per Cell", xlab="Cells Per Bin")

png(file.path(getOutputDirectory(aged_MuSC_projAging1),"Plots","aged_MuSC_covered_peaks_per_bin.png"),res=300,units='in',width=5,height=5)
boxplot(lapply(aged_atac_cds.cicero.klist, function(x) colData(x[[1]])$num_genes_expressed), 
        names = seq(10,60,5), xlab = "Cells Per Bin", ylab = "# Covered Peaks Per Bin", main="Aged MuSCs")
dev.off()

png(file.path(getOutputDirectory(young_MuSC_projAging1),"Plots","young_MuSC_covered_peaks_per_bin.png"),res=300,units='in',width=5,height=5)
boxplot(lapply(young_atac_cds.cicero.klist, function(x) colData(x[[1]])$num_genes_expressed), 
        names = seq(10,60,5), xlab = "Cells Per Bin", ylab = "# Covered Peaks Per Bin", main="Young MuSCs")
dev.off()

png(file.path(getOutputDirectory(MuSC_projAging1),"Plots","MuSC_covered_peaks_per_bin.png"),res=300,units='in',width=5,height=5)
boxplot(lapply(MuSC_atac_cds.cicero.klist, function(x) colData(x[[1]])$num_genes_expressed), 
        names = seq(10,60,5), xlab = "Cells Per Bin", ylab = "# Covered Peaks Per Bin", main="All MuSCs")
dev.off()

# png(file.path(projdir,"integrated_data_cicero","young_MuSC_cellaggregates_per_bin.png"),res=300,units='in',width=5,height=5)
# barplot(lapply(unlist(young_atac_cds.cicero.klist), function(x) nrow(colData(x))), names.arg = seq(10,60,5), xpd = F,
#         xlab = "Cells Per Bin", ylab = "# Cell Aggregates", main="Young MuSCs")
# dev.off()
# png(file.path(projdir,"integrated_data_cicero","aged_MuSC_cellaggregates_per_bin.png"),res=300,units='in',width=5,height=5)
# barplot(unlist(lapply(aged_MuSC_cds.cicero.klist,nrow)), names.arg = seq(10,60,5), xpd = F,
#         xlab = "Cells Per Bin", ylab = "# Cell Aggregates", main="Aged MuSCs")
# dev.off()

# Make cicero cds objects ------------------------------------------------------

young_atac_cds.cicero = make_cicero_cds(young_atac_cds, reduced_coordinates = young_umap_coords, k = 35)
aged_atac_cds.cicero = make_cicero_cds(aged_atac_cds, reduced_coordinates = aged_umap_coords, k = 20)
MuSC_atac_cds.cicero = make_cicero_cds(MuSC_atac_cds, reduced_coordinates = MuSC_umap_coords, k = 50)

saveRDS(young_atac_cds.cicero, file.path(getOutputDirectory(young_MuSC_projAging1), "young_MuSC_cicero.RDS"))
saveRDS(aged_atac_cds.cicero, file.path(getOutputDirectory(aged_MuSC_projAging1), "aged_MuSC_cicero.RDS"))
saveRDS(MuSC_atac_cds.cicero, file.path(getOutputDirectory(MuSC_projAging1), "all_MuSC_cicero.RDS"))

# Add peak2gene linkages ------------------------------------------------------

source("/nas/homes/benyang/HiC/13_MultiOme/custom_addPeak2GeneLinks.R")

MuSC_projAging1 <- addPeak2GeneLinks.mod(
    ArchRProj = MuSC_projAging1,
    reducedDims = "IterativeLSI",
    addEmpiricalPval = FALSE,
    maxDist = 5e5,
    group1 = grep("Young",getCellNames(MuSC_projAging1),fixed=T,value=T),
    group2 = grep("Aged",getCellNames(MuSC_projAging1),fixed=T,value=T)
)

p2g = getPeak2GeneLinks(MuSC_projAging1, returnLoops=F)
saveRDS(p2g, "/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/MuSC_ArchR/all_MuSC/Peak2GeneLinks/all_MuSC_peak2gene.RDS")

# p4 <- plotPeak2GeneHeatmap(ArchRProj = MuSC_projAging3, groupBy = "Age", k=20, nPlot=10000, returnMatrices = TRUE)
# p5 <- plotPeak2GeneHeatmap(ArchRProj = MuSC_projAging3, groupBy = "Age", k=20, nPlot=10000, returnMatrices = FALSE)
p6_files <- plotPeak2GeneHeatmap(ArchRProj = MuSC_projAging1, groupBy = "Age", k=5, nPlot=dim(p2g)[1], returnMatrices = TRUE)
p6 <- plotPeak2GeneHeatmap(ArchRProj = MuSC_projAging1, palGroup = c("Young"=ggsci::pal_nejm()(2)[2], "Aged"=ggsci::pal_nejm()(2)[1]), groupBy = "Age", k=5, nPlot=dim(p2g)[1], returnMatrices = FALSE)

plotPDF(p6, name = "Plots-Peak-2-Gene-MuSC-Age_custom.pdf", ArchRProj = MuSC_projAging1, addDOC = FALSE, width = 5, height = 7)

saveRDS(p6_files, "/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/MuSC_ArchR/all_MuSC/Peak2GeneLinks/all_MuSC_peak2gene_clustered.RDS")

p2g_list = lapply(unique(p6_files$ATAC$kmeansId), function(x) rownames(p6_files$ATAC$matrix)[which(p6_files$ATAC$kmeansId==x)])
names(p2g_list) = unique(p6_files$ATAC$kmeansId)

saveArchRProject(ArchRProj = MuSC_projAging1, outputDirectory = "Save-MuSC_projAging1-01", load = FALSE)

## K = 20
# 1 - 
# 2 - Chromosome segregation
# 7 - Cell cycle
# 8 - Mitochondrial translation
# 10 - Mitochondrial translation
# 13 - Mitochondrial translation
# 18 - cellular response to stress, cell cycle, chromosome organization
# 19 - Regulation of necroptotic cell death, regulation of TNFR1 signaling, TNF singaling, heat stress/shock, rRNA processing, TP53, stimuli, transcription

## K = 10
# 1 - NOTCH3 signaling, mitochondrial translation
# 2 - Mitochondrial translation
# 4 - Mitochondrial translation
# 5 - senescence, stress, cell cycle, transcription
# 6 - resolution of AP sites, base excision repair, mitochondrial translation, SUMOylation, cell cycle, chromosome segregation, cell cycle, transcription
# 7 - ribosomal scanning, hedgehog signaling, cell cycle, cell stress
# 8 - TLR4 signaling, IFNa, NFkB, mRNA splicing, integrin, SUMOylation, cell stimuli, transcription, chromatin organization
# 9 - transcription, TP53
# 10 - ribonucleoprotein organization/biogenesis, metabolism

# output genes per kmeans group for GO term enrichment analysis
for(n in names(p2g_list)) {
    write.table(unique(p6_files$Peak2GeneLinks[p2g_list[[n]], 'gene']), file.path(getOutputDirectory(MuSC_projAging1), "Peak2GeneLinks", paste0("MuSC_projAging1_p2g_c",n,".txt")), row.names=F, col.names=F, quote=F)
}

gene_cluster_list = lapply(names(p2g_list), function(x) p6_files$Peak2GeneLinks[p2g_list[[x]], 'gene'])
gene_cluster_df = data.frame(gene = unlist(gene_cluster_list))
gene_cluster_df$Cluster = rep(names(p2g_list), sapply(gene_cluster_list, length))

# get average expression per Age group 
young_MuSC_cellnames = grep("Young",getCellNames(MuSC_projAging1),fixed=T,value=T)
aged_MuSC_cellnames = grep("Aged",getCellNames(MuSC_projAging1),fixed=T,value=T)

exp_mat = getMatrixFromProject(MuSC_projAging1)

get_avg_exp = function(pattern, exp_mat, young_MuSC_cellnames, aged_MuSC_cellnames) {
    pattern_idx = grep(pattern,rowData(exp_mat)$name,fixed=T)
    pattern_avg_exp = data.frame(Young = Matrix::rowMeans(assay(exp_mat)[pattern_idx,young_MuSC_cellnames]),
                                Aged = Matrix::rowMeans(assay(exp_mat)[pattern_idx,aged_MuSC_cellnames]),
                                Young_sem = rowSds(assay(exp_mat)[pattern_idx,young_MuSC_cellnames]) / sqrt(length(young_MuSC_cellnames)),
                                Aged_sem = rowSds(assay(exp_mat)[pattern_idx,aged_MuSC_cellnames]) / sqrt(length(aged_MuSC_cellnames)),
                                row.names = rowData(exp_mat)$name[pattern_idx])
    return(pattern_avg_exp)
}

mrpl_df = get_avg_exp("Mrpl",exp_mat,young_MuSC_cellnames,aged_MuSC_cellnames)
ndufa_df = get_avg_exp("Ndufa",exp_mat,young_MuSC_cellnames,aged_MuSC_cellnames)

# Add gene names to the Peak2Gene heatmap

kdm_idx = grep("Kdm",mcols(metadata(p6_files$Peak2GeneLinks)$geneSet)$name,fixed=T) # index of Kdm genes in the gene set
kmt_idx = grep("Kmt",mcols(metadata(p6_files$Peak2GeneLinks)$geneSet)$name,fixed=T) # index of Kmt genes in the gene set
nduf_idx = grep("Nduf",mcols(metadata(p6_files$Peak2GeneLinks)$geneSet)$name,fixed=T) # index of Kmt genes in the gene set
notch_idx = which(mcols(metadata(p6_files$Peak2GeneLinks)$geneSet)$name %in% c("Notch3","Dll1","Jag2","Aph1a","Aph1b","Aph1c")) # index of Kmt genes in the gene set
# gene_idx = c(kdm_idx, kmt_idx, notch_idx)
# gene_p2g_idx = which(p6_files$Peak2GeneLinks$idxRNA %in% gene_idx) # location of Kdm genes in the peak-to-gene linkages
# gene_names = p6_files$Peak2GeneLinks$gene[gene_p2g_idx]
kdm_p2g_idx = (p6_files$Peak2GeneLinks$idxRNA %in% kdm_idx) + 0
nduf_p2g_idx = (p6_files$Peak2GeneLinks$idxRNA %in% nduf_idx) + 0

png(file.path(getOutputDirectory(MuSC_projAging1),"Peak2GeneLinks","annotated_p2g_heatmap.png"), units='in', res=300, width = 5, height = 7)
p6 + 
    Heatmap(kdm_p2g_idx[order(p6_files$RNA$kmeansId)], name="Kdm", width=unit(0.5,"cm"), col = c("0" = "white", "1" = "purple"), show_heatmap_legend = FALSE) +
    Heatmap(nduf_p2g_idx[order(p6_files$RNA$kmeansId)], name="Nduf", width=unit(0.5,"cm"), col = c("0" = "white", "1" = "green"), show_heatmap_legend = FALSE) +
    rowAnnotation(link = anno_mark(at = which(nduf_p2g_idx[order(p6_files$RNA$kmeansId)]=='1'), 
        labels = rep_len("hello", length(which(nduf_p2g_idx[order(p6_files$RNA$kmeansId)]=='1'))), 
        labels_gp = gpar(fontsize = 5), padding = unit(1, "mm")))
dev.off()

# Make browser track plot only peak2gene links involving target gene
aged_TADs = read.table("/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed", sep='\t')[,1:3]
colnames(aged_TADs) = c("seqnames","start","end")
young_TADs = read.table("/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed", sep='\t')[,1:3]
colnames(young_TADs) = c("seqnames","start","end")

aged_TADs = GRanges(aged_TADs)
young_TADs = GRanges(young_TADs)

young_conns_gi = readRDS("/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/MuSC_ArchR/young_MuSC/young_conns_gi_peakmatrix.RDS")
aged_conns_gi = readRDS("/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/MuSC_ArchR/aged_MuSC/aged_conns_gi_peakmatrix.RDS")

make_conns_plt_df = function(conns_gi) {
    conns_gi_left = resize(GenomicInteractions::anchorOne(conns_gi), 1, 'center')
    conns_gi_right = resize(GenomicInteractions::anchorTwo(conns_gi), 1, 'center')
    conns_gi_plot = GRanges(data.frame(seqnames=seqnames(conns_gi_left), start=start(conns_gi_left), end=end(conns_gi_right)))
    mcols(conns_gi_plot)$value = mcols(conns_gi)$coaccess
    return(conns_gi_plot)
}

young_conns_gi_plot = make_conns_plt_df(young_conns_gi)
aged_conns_gi_plot = make_conns_plt_df(aged_conns_gi)

atac_ranges = getPeak2GeneLinks(MuSC_projAging1)[[1]]
gene_position = start(metadata(p2g)$geneSet[metadata(p2g)$geneSet$name=="Ndufc1",])
atac_ranges_subset = atac_ranges[start(atac_ranges)==gene_position | end(atac_ranges)==gene_position,]

plt_region = GRanges(data.frame(seqnames=unique(as.character(seqnames(atac_ranges_subset))), 
                        start = min(start(atac_ranges_subset))-5e4,
                        end = max(end(atac_ranges_subset))+5e4))
aged_TAD_feature = subsetByOverlaps(aged_TADs, plt_region, ignore.strand = T)
aged_TAD_feature$name=paste0("Aged TAD",seq_len(length(aged_TAD_feature)))
aged_TAD_feature = data.frame(aged_TAD_feature)
young_TAD_feature = subsetByOverlaps(young_TADs, plt_region, ignore.strand = T)
young_TAD_feature$name=paste0("Young TAD",seq_len(length(young_TAD_feature)))
young_TAD_feature = data.frame(young_TAD_feature)
all_TAD_feature = rbind(aged_TAD_feature, young_TAD_feature)
all_TAD_feature$start = unlist(lapply(all_TAD_feature$start, function(x) ifelse(x<start(plt_region), start(plt_region), x)))
all_TAD_feature$end = unlist(lapply(all_TAD_feature$end, function(x) ifelse(x>end(plt_region), end(plt_region), x)))

p <- plotBrowserTrack(
    ArchRProj = MuSC_projAging1, 
    groupBy = "Age", 
    baseSize = 7,
    tileSize = 1e3,
    sizes = c(8, 1.5, 3, 4),
    region = plt_region,
    loops = GRangesList(P2G = atac_ranges_subset, 
                        Young_conns = young_conns_gi_plot[which(young_conns_gi_plot$value>0.1),],
                        Aged_conns = aged_conns_gi_plot[which(aged_conns_gi_plot$value>0.1),]),
    features = GRangesList(Peaks = getPeakSet(MuSC_projAging1), 
                            young_TAD = GRanges(all_TAD_feature[3:4,]), 
                            aged_TAD = GRanges(all_TAD_feature[1:2,]))
    #features = young_TADs,
    #loops = young_conns_gi_plot[which(young_conns_gi_plot$value>0.1),]
    #features = aged_TADs,
    #loops = aged_conns_gi_plot[which(aged_conns_gi_plot$value>0.1),]
)
png(file.path(getOutputDirectory(MuSC_projAging1),"Peak2GeneLinks","Ndufc1_p2g_conns_plot.png"), units='in', res=300, width = 5.5, height = 4.5)
grid::grid.newpage()
grid::grid.draw(p)
dev.off()

# featureO$name = as.character(featureO$name)
# featureO$name[32] = "TAD2"
# featureO$name <- factor(paste0(featureO$name), levels=c("Peaks","TAD","TAD2"))
# featureO$facet = "FeatureTracks"

# featureO[featureO$name != "Peaks", 'start'] = unlist(lapply(featureO[featureO$name != "Peaks", 'start'], 
#     function(x) ifelse(x<start(region), start(region), x) ))
# featureO[featureO$name != "Peaks", 'end'] = unlist(lapply(featureO[featureO$name != "Peaks", 'end'], 
#     function(x) ifelse(x>end(region), end(region), x) ))

# ggplot(data=featureO, aes(color=name)) +
# facet_grid(facet~.) +
# geom_segment(data=featureO, aes(x=start,xend=end,y=name,yend=name,color=name), size=2) +
# scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0))

getGroupBW(
  ArchRProj = MuSC_projAging1,
  groupBy = "Age",
  normMethod = "ReadsInTSS",
  tileSize = 1e2,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)


# mATAC = readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ]
# mRNA = readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ]

# mATAC = assay(mATAC)
# mRNA = assay(mRNA)

# .rowZscores <- function(m = NULL, min = -2, max = 2, limit = FALSE){
#   z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m),`/`)
#   if(limit){
#     z[z > max] <- max
#     z[z < min] <- min
#   }
#   return(z)
# }

# mATAC = .rowZscores(mATAC)
# mRNA = .rowZscores(mRNA)

# rownames(mATAC) <- NULL
# rownames(mRNA) <- NULL
# colnames(mATAC) <- paste0("K_", seq_len(ncol(mATAC)))
# colnames(mRNA) <- paste0("K_", seq_len(ncol(mRNA)))
# rownames(mATAC) <- paste0("P2G_", seq_len(nrow(mATAC)))
# rownames(mRNA) <- paste0("P2G_", seq_len(nrow(mRNA)))
# rownames(p2g) <- paste0("P2G_", seq_len(nrow(p2g)))

# mATAC_kmeans = lapply(seq(50,100,10), function(k) kmeans(mATAC, centers=k, iter.max=20))
# mATAC_kmeans_fit = lapply(mATAC_kmeans, function(k) k$betweenss/k$totss)

# mATAC_kmeans_nbclust = NbClust::NbClust(mATAC,min.nc=5,max.nc=25,method='kmeans',index='silhouette')

# set.seed(123)
# fviz_nbclust(mATAC, kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+
#   labs(subtitle = "Gap statistic method")

bind_rows(lapply(metadata(seATAC)$KNNList, function(i) table(sapply(i, function(x) ifelse(grepl("Young",x,fixed=T),"Young","Aged"))))) %>%
mutate(cluster = 1:n()) %>%
pivot_longer(cols=c(Young,Aged)) %>%
ggplot(aes(x=cluster, y=value)) +
geom_col(aes(fill=name))

saveArchRProject(ArchRProj = MuSC_projAging2, outputDirectory = "Save-MuSC_projAging2-01", load = FALSE)
saveArchRProject(ArchRProj = young_MuSC_projAging1, outputDirectory = "Save-young_MuSC_projAging1-01", load = FALSE)
saveArchRProject(ArchRProj = aged_MuSC_projAging1, outputDirectory = "Save-aged_MuSC_projAging1-01", load = FALSE)

p1 <- plotPeak2GeneHeatmap(ArchRProj = young_MuSC_projAging1, groupBy = "Age", k=10, nPlot=10000)
p2 <- plotPeak2GeneHeatmap(ArchRProj = MuSC_projAging2, groupBy = "Age", k=30, nPlot=10000)
p3 <- plotPeak2GeneHeatmap(ArchRProj = MuSC_projAging2, groupBy = "Sample")

plotPDF(p1, name = "Plots-Peak-2-Gene-Young-MuSC.pdf", ArchRProj = young_MuSC_projAging1, addDOC = FALSE, width = 7, height = 20)
plotPDF(p2, name = "Plots-Peak-2-Gene-MuSC-Age.pdf", ArchRProj = MuSC_projAging1, addDOC = FALSE, width = 3, height = 5)
plotPDF(p3, name = "Plots-Peak-2-Gene-MuSC-Sample.pdf", ArchRProj = MuSC_projAging1, addDOC = FALSE, width = 10, height = 25)

png(file.path(getOutputDirectory(MuSC_projAging1),"Plots","Plots-Peak-2-Gene-MuSC-Age.png"),res=300,units='in',width=10,height=25)
print(p2)
dev.off()

png(file.path(getOutputDirectory(projAging6),"Plots","Plots-Peak-2-Gene-All.png"),res=300,units='in',width=10,height=25)
print(p3)
dev.off()

markersPeaks <- getMarkerFeatures(
    ArchRProj = MuSC_projAging1, 
    useMatrix = "PeakMatrix", 
    groupBy = "Age",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  plotLog2FC = T,
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

markerTest <- getMarkerFeatures(
  ArchRProj = MuSC_projAging1, 
  useMatrix = "PeakMatrix",
  groupBy = "Age",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Young",
  bgdGroups = "Aged",
  binarize = T
)

pma <- plotMarkers(seMarker = markerTest, name = "Young", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pma <- plotMarkers(seMarker = markerTest, name = "Young", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pma

png(file.path(getOutputDirectory(MuSC_projAging1),"Plots","MuSC-Diff-Age-Peaks.png"),res=300,units='in',width=5,height=5)
pma + scale_color_manual(values=c(ggsci::pal_nejm()(2)[1], "lightgrey", ggsci::pal_nejm()(2)[2]))
dev.off()