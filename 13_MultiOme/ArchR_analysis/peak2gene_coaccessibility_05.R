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

projAging5 = readRDS(file.path(projdir, "Save-projAging5-01", "Save-ArchR-Project.rds"))

# Add co-accessibility  ------------------------------------------------------

projAging6 <- addCoAccessibility(
    ArchRProj = projAging5,
    reducedDims = "Harmony"
)

cA <- getCoAccessibility(
    ArchRProj = projAging6,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
)

p <- plotBrowserTrack(
    ArchRProj = projAging6, 
    groupBy = "celltype", 
    geneSymbol = "Pax7", 
    upstream = 1.5e5,
    downstream = 5e4,
    loops = getCoAccessibility(projAging6)
)
grid::grid.newpage()
grid::grid.draw(p$Pax7)

plotPDF(p$Pax7, name = "Plot-Pax7-Coaccessibility-All-Cells.pdf", ArchRProj = projAging6, addDOC = FALSE, width = 5, height = 5)

# Add peak 2 gene linkage  ------------------------------------------------------

projAging6 <- addPeak2GeneLinks(
    ArchRProj = projAging5,
    reducedDims = "Harmony"
)

p <- plotPeak2GeneHeatmap(ArchRProj = projAging6, groupBy = "celltype", nPlot=10000, k=20)
plotPDF(p, name = "Plots-Peak-2-Gene-celltype.pdf", ArchRProj = projAging6, addDOC = FALSE, width = 10, height = 25)

saveArchRProject(ArchRProj = projAging6, outputDirectory = "Save-projAging6-01", load = FALSE)