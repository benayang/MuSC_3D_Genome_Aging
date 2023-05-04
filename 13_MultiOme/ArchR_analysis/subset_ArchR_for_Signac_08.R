library(ArchR)
library(dplyr)
library(tidyr)
library(parallel)
library(SingleCellExperiment)
library(Signac)
library(Seurat)

#ArchR::installExtraPackages()
projdir = '/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis'

addArchRGenome("mm10")
addArchRThreads(threads = 45) 

projAging6 = readRDS(file.path(projdir, "Save-projAging6-01", "Save-ArchR-Project.rds"))

# subset by age ------------------------------------------------------

peakmatrix = getMatrixFromProject(projAging6, useMatrix="PeakMatrix")
rD = getReducedDims(projAging6)

young_cellnames = grep("Young#", getCellNames(projAging6), value=T, fixed=T)
young_v2_cellnames = grep("Young_v2#", getCellNames(projAging6), value=T, fixed=T)
aged_cellnames = grep("Aged#", getCellNames(projAging6), value=T, fixed=T)

young_se = as(peakmatrix[,young_cellnames], "SingleCellExperiment")
young_v2_se = as(peakmatrix[,young_v2_cellnames], "SingleCellExperiment")
aged_se = as(peakmatrix[,aged_cellnames], "SingleCellExperiment")

se_list = list(young_se, young_v2_se, aged_se)
cellnames_list = list(young_cellnames, young_v2_cellnames, aged_cellnames)
obj_list = list()

for (se in 1:length(se_list)) {
    names(assays(se_list[[se]])) = "counts"
    rownames(se_list[[se]]) = paste0(seqnames(se_list[[se]]), ":", start(se_list[[se]]),"-",end(se_list[[se]]))
    obj_list[[se]] = as.Seurat(se_list[[se]], data=NULL)
    obj_list[[se]][["iterativeLSI"]] = CreateDimReducObject(embeddings = rD[cellnames_list[[se]], ], key = "ILSI_", assay = "ATAC")
}

names(assays(young_se)) = "counts"
names(assays(young_v2_se)) = "counts"
names(assays(aged_se)) = "counts"

rownames(young_se) <- paste0(seqnames(young_se), ":", start(young_se),"-",end(young_se))