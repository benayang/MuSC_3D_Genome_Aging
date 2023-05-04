library(ArchR)
library(dplyr)
library(tidyr)
library(parallel)

#ArchR::installExtraPackages()
projdir = '/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis'

addArchRGenome("mm10")
addArchRThreads(threads = 45) 

projAging3 = readRDS(file.path(projdir, "Save-projAging3-01", "Save-ArchR-Project.rds"))

scRNA = readRDS("/nas/homes/benyang/HiC/13_MultiOme/TabulaMurisSenis/aging_RNA_integrated_annotated.RDS")

# Unconstrained and constrained RNA integration ------------------------------------------------------

projAging4 <- addGeneIntegrationMatrix(
    ArchRProj = projAging3, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = scRNA,
    addToArrow = FALSE,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

cM <- as.matrix(confusionMatrix(projAging4$Harmony_Clusters, projAging4$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
preClust_assignments = cbind(preClust, rownames(cM)) #Assignments

unique(unique(projAging4$predictedGroup_Un))

# atac_MuSC <- "C12"
# atac_Tenocyte <- "C14"
# atac_MSC <- "C13"

atac_cellnames = lapply(c("C12","C13","C14"), function(x) projAging4$cellNames[projAging4$Harmony_Clusters == x])
names(atac_cellnames) = c("MuSC","MSC","Tenocyte")
atac_non_cellnames = lapply(c("C12","C13","C14"), function(x) projAging4$cellNames[projAging4$Harmony_Clusters != x])
names(atac_non_cellnames) = c("MuSC","MSC","Tenocyte")

rna_cellnames = lapply(c("MuSC","Tenocyte","MSC"), function(x) colnames(scRNA)[grep(x, scRNA$celltype, fixed=T)])
names(rna_cellnames) = c("MuSC","Tenocyte","MSC")
rna_non_cellnames = lapply(c("MuSC","Tenocyte","MSC"), function(x) colnames(scRNA)[grep(x, scRNA$celltype, fixed=T, invert=T)])
names(rna_non_cellnames) = c("MuSC","Tenocyte","MSC")

groupList <- SimpleList(
    MuSC = SimpleList(
        ATAC = atac_cellnames[['MuSC']],
        RNA = rna_cellnames[['MuSC']]
    ),
    NonMuSC = SimpleList(
        ATAC = atac_non_cellnames[['MuSC']],
        RNA = rna_non_cellnames[['MuSC']]
    )
)

# groupList <- SimpleList(
#     MuSC = SimpleList(
#         ATAC = atac_cellnames[['MuSC']],
#         RNA = rna_cellnames[['MuSC']]
#     ),
#     NonMuSC = SimpleList(
#         ATAC = atac_non_cellnames[['MuSC']],
#         RNA = rna_non_cellnames[['MuSC']]
#     ),
#     Tenocyte = SimpleList(
#         ATAC = atac_cellnames[['Tenocyte']],
#         RNA = rna_cellnames[['Tenocyte']]
#     ),
#     NonTenocyte = SimpleList(
#         ATAC = atac_non_cellnames[['Tenocyte']],
#         RNA = rna_non_cellnames[['Tenocyte']]
#     ),
#     MSC = SimpleList(
#         ATAC = atac_cellnames[['MSC']],
#         RNA = rna_cellnames[['MSC']]
#     ),
#     NonMSC = SimpleList(
#         ATAC = atac_non_cellnames[['MSC']],
#         RNA = rna_non_cellnames[['MSC']]
#     )
# )


projAging4 <- addGeneIntegrationMatrix(
    ArchRProj = projAging4, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = scRNA,
    addToArrow = FALSE, 
    groupList = groupList,
    groupRNA = "celltype",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)

p1 = plotEmbedding(projAging4, embedding="UMAP_Harmony", colorBy="cellColData", name="predictedGroup_Un")
p2 = plotEmbedding(projAging4, embedding="UMAP_Harmony", colorBy="cellColData", name="predictedGroup_Co")
p3 = plotEmbedding(projAging4, embedding="UMAP_Harmony", colorBy="cellColData", name="Harmony_Clusters")
ggAlignPlots(p1, p2, p3, type='h')

# Actually add the RNA integration to the object ------------------------------------------------------

projAging4 <- addGeneIntegrationMatrix(
    ArchRProj = projAging4, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = scRNA,
    addToArrow = T, 
    force = T,
    groupList = groupList,
    groupRNA = "celltype",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)

saveArchRProject(ArchRProj = projAging4, outputDirectory = "Save-projAging4-01", load = FALSE)

pal <- paletteDiscrete(values = scRNA$celltype)
p1 <- plotEmbedding(
    projAging4, 
    embedding = "UMAP_Harmony",
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    pal = pal
)
p2 <- plotEmbedding(
    projAging4, 
    embedding = "UMAP_Harmony",
    colorBy = "cellColData", 
    name = "predictedGroup_Co", 
    pal = pal
)
p3 <- plotEmbedding(
    projAging4, 
    embedding = "UMAP_Harmony",
    colorBy = "GeneScoreMatrix", 
    name = c("Pax7","Myod1"), 
    continuousSet = "horizonExtra",
    imputeWeights = getImputeWeights(projAging4)
)
p4 <- plotEmbedding(
    projAging4, 
    embedding = "UMAP_Harmony",
    colorBy = "GeneIntegrationMatrix", 
    name = c("Pax7","Myod1"), 
    continuousSet = "horizonExtra",
    imputeWeights = getImputeWeights(projAging4)
)
p5 <- plotEmbedding(
    projAging4, 
    embedding = "UMAP_Harmony",
    colorBy = "cellColData", 
    name = "Age",
    highlightCells = grep("Aged", projAging4$cellNames, fixed=T, value=T),
    pal = c(Aged = ggsci::pal_nejm()(2)[1])
)
p6 <- plotEmbedding(
    projAging4, 
    embedding = "UMAP_Harmony",
    colorBy = "cellColData", 
    name = "Age",
    highlightCells = grep("Young", projAging4$cellNames, fixed=T, value=T),
    pal = c(Young = ggsci::pal_nejm()(2)[2])
)

plotPDF(p1,p2,p3,p4,p5,p6, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = projAging4, addDOC = FALSE, width = 5, height = 5)
