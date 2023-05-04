library(dplyr)
library(tidyr)
library(ArchR)
library(plotgardener)
library(GenomicInteractions)

projdir = "/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/MuSC_ArchR"

# USE CICERO CONNECTIONS PROCESSED ON LAPTOP
# ## Get Cicero connections
# load(file.path(projdir, "aged_MuSC", "Aged_peakmatrix_MuSC_cicero_conns_ccans.RData"))
# aged_conns = conns
# aged_ccans = ccans
# load(file.path(projdir, "young_MuSC", "Young_peakmatrix_MuSC_cicero_conns_ccans.RData"))
# young_conns = conns
# young_ccans = ccans

# rm(list=c("conns","ccans"))
# gc()

# # Convert conns to GenomicInteractions object
# conns_to_gi = function(conns) {
#         peak1_GR = Signac::StringToGRanges(conns$Peak1, sep=c("_","_"))
#         peak2_GR = Signac::StringToGRanges(conns$Peak2, sep=c("_","_"))
#         conns_GI = GenomicInteractions(peak1_GR, peak2_GR, mode="strict")
#         mcols(conns_GI)$coaccess = conns$coaccess
#         conns_GI = unique(conns_GI)
#         return(conns_GI)
# }

# young_conns_gi = conns_to_gi(young_conns)
# saveRDS(young_conns_gi, file.path(projdir, "plotgardener", "young_conns_GI_peakmatrix.RDS"))
# aged_conns_gi = conns_to_gi(aged_conns)
# saveRDS(aged_conns_gi, file.path(projdir, "plotgardener", "aged_conns_GI_peakmatrix.RDS"))

## Get MACS2 MuSC peaks
young_MuSC_projAging1 = readRDS("/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/Save-young_MuSC_projAging1-01/Save-ArchR-Project.rds")
aged_MuSC_projAging1 = readRDS("/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/Save-aged_MuSC_projAging1-01/Save-ArchR-Project.rds")

young_MuSC_peaks = getPeakSet(young_MuSC_projAging1)
aged_MuSC_peaks = getPeakSet(aged_MuSC_projAging1)

saveRDS(young_MuSC_peaks, file.path(projdir, "plotgardener", "young_MuSC_peaks.RDS"))
saveRDS(aged_MuSC_peaks, file.path(projdir, "plotgardener", "aged_MuSC_peaks.RDS"))

## Get GeneIntegrationMatrix
young_gene_matrix = getMatrixFromProject(young_MuSC_projAging1, "GeneIntegrationMatrix")
aged_gene_matrix = getMatrixFromProject(aged_MuSC_projAging1, "GeneIntegrationMatrix")

saveRDS(young_gene_matrix, file.path(projdir, "plotgardener", "young_MuSC_genematrix.RDS"))
saveRDS(aged_gene_matrix, file.path(projdir, "plotgardener", "aged_MuSC_genematrix.RDS"))

## Get MuSC objects with LinkPeaks run
MuSC_projAging1 = readRDS("/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/Save-MuSC_projAging1-01/Save-ArchR-Project.rds")
p2g = getPeak2GeneLinks(MuSC_projAging1, corCutOff=0, resolution=1, returnLoops=T)[[1]]
# young_MuSC_links = getPeak2GeneLinks(young_MuSC_projAging1, corCutOff=0, resolution=1, returnLoops=T)[[1]]
# aged_MuSC_links = getPeak2GeneLinks(aged_MuSC_projAging1, corCutOff=0, resolution=1, returnLoops=T)[[1]]

## Convert LinkPeaks to bedpe files 
p2g_to_bedpe = function(links) {
        links_bedpe = links %>% 
                        as.data.frame() %>% 
                        dplyr::select(seqnames, start, end, value) %>%
                        mutate(start1 = start, end1 = start + 1,
                                start2 = end, end2 = end + 1)
        links_bedpe = cbind(links_bedpe[,c("seqnames","start1","end1")], links_bedpe[,c("seqnames","start2","end2","value")])
        colnames(links_bedpe)[1] = "chrom1"
        colnames(links_bedpe)[4] = "chrom2"
        links_bedpe_GI = GenomicInteractions(GRanges(setNames(links_bedpe[,1:3], c("chrom","start","end"))),
                                                GRanges(setNames(links_bedpe[,4:6], c("chrom","start","end"))),
                                                counts = links_bedpe[,7],
                                                mode = "strict")
        return(list(bedpe = links_bedpe, bedpe_GI = links_bedpe_GI))
}
p2g_bedpe = p2g_to_bedpe(p2g)


# boxplot(list(young = mcols(young_links_bedpe_GI)$counts, aged=mcols(aged_links_bedpe_GI)$counts))

saveRDS(p2g_bedpe$bedpe, file.path(projdir, "plotgardener", "MuSC_p2g_bedpe.RDS"))
saveRDS(p2g_bedpe$bedpe_GI, file.path(projdir, "plotgardener", "MuSC_p2g_GI.RDS"))
