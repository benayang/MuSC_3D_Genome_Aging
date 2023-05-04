library(Seurat)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)
library(future)

plan("multicore", workers=45)
options(future.globals.maxSize = 75 * 1024^3)

projdir = '/nas/homes/benyang/HiC/13_MultiOme'

aged_MuSC = readRDS(file.path(projdir, "all_reps_aged_MuSC_subset.RDS"))
young_MuSC = readRDS(file.path(projdir, "all_reps_young_MuSC_subset.RDS"))

# first compute the GC content for each peak
DefaultAssay(aged_MuSC) = "ATAC"
DefaultAssay(young_MuSC) = "ATAC"

main.chroms = standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
young.keep.peaks <- which(as.character(seqnames(granges(young_MuSC))) %in% main.chroms)
young_MuSC[["ATAC"]] <- subset(young_MuSC[["ATAC"]], features = rownames(young_MuSC[["ATAC"]])[young.keep.peaks])
aged.keep.peaks <- which(as.character(seqnames(granges(aged_MuSC))) %in% main.chroms)
aged_MuSC[["ATAC"]] <- subset(aged_MuSC[["ATAC"]], features = rownames(aged_MuSC[["ATAC"]])[aged.keep.peaks])

aged_MuSC = RegionStats(aged_MuSC, genome = BSgenome.Mmusculus.UCSC.mm10)
young_MuSC = RegionStats(young_MuSC, genome = BSgenome.Mmusculus.UCSC.mm10)

stopifnot(identical(rownames(aged_MuSC[["RNA"]]), rownames(young_MuSC[["RNA"]])))
known.genes = rownames(aged_MuSC[["RNA"]])
known.genes = known.genes[!grepl("^Gm", known.genes) & !grepl("Rik", known.genes)]

aged_MuSC = LinkPeaks(aged_MuSC, 
                        peak.assay = "ATAC",
                        peak.slot = "data",
                        expression.assay = "RNA",
                        distance = 5e5,
                        min.cells = 5,
                        n_sample = 200,
                        genes.use = known.genes)

young_MuSC = LinkPeaks(young_MuSC, 
                        peak.assay = "ATAC",
                        peak.slot = "data",
                        expression.assay = "RNA",
                        distance = 5e5,
                        min.cells = 5,
                        n_sample = 200,
                        genes.use = known.genes)

saveRDS(aged_MuSC, file.path(projdir, "all_reps_aged_MuSC_subset_linkPeaks.RDS"))
saveRDS(young_MuSC, file.path(projdir, "all_reps_young_MuSC_subset_linkPeaks.RDS"))