library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

plan("multicore", workers = 30)
options(future.globals.maxSize = 70 * 1024 ^ 3) # for 70 Gb RAM

data = readRDS("/nas/homes/benyang/HiC/13_MultiOme/aging_integrated.RDS")
DefaultAssay(data) = "ATAC"

# create blacklist object
blacklist = read.table("/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed", sep='\t')
colnames(blacklist) = c("chrom","start","end","name")
blacklist = GRanges(blacklist)

# split integrated data into Young and Aged
age_split_data = SplitObject(data, split.by="orig.ident")

young.counts = age_split_data$Young@assays$ATAC@counts
colnames(young.counts) = sapply(colnames(young.counts), function(x) gsub("_1","",x,fixed=T))
aged.counts = age_split_data$Aged@assays$ATAC@counts
colnames(aged.counts) = sapply(colnames(aged.counts), function(x) gsub("_2","",x,fixed=T))

young_assay = CreateChromatinAssay(counts = young.counts, fragments = "/nas/homes/benyang/HiC/13_MultiOme/RawData/Young/10x_analysis_5644-JL/Sample_5644-JL-1/atac_fragments.tsv.gz")
aged_assay = CreateChromatinAssay(counts = aged.counts, fragments = "/nas/homes/benyang/HiC/13_MultiOme/RawData/Aged/10x_analysis_5422-JL/Sample_5422-JL-1/atac_fragments.tsv.gz")

young_md = age_split_data$Young@meta.data
rownames(young_md) = sapply(rownames(young_md), function(x) gsub("_1","",x,fixed=T))
aged_md = age_split_data$Aged@meta.data
rownames(aged_md) = sapply(rownames(aged_md), function(x) gsub("_2","",x,fixed=T))

young_ATAC_obj = CreateSeuratObject(counts = young_assay, assay="ATAC", project="Young", meta.data=young_md)
aged_ATAC_obj = CreateSeuratObject(counts = aged_assay, assay="ATAC", project="Aged", meta.data=aged_md)

# Call peaks to replace 10X peaks
call_peaks = function(data, name, blacklist, outdir) {
    peaks = CallPeaks(
        data, 
        macs2.path="/opt/miniconda3/envs/py38/bin/macs2", 
        name=name,
        group.by="celltype_wsnn",
        combine.peaks = FALSE,
        effective.genome.size=1.87e9, 
        outdir=outdir, 
        cleanup=F, 
        format = "BED",
        extsize = 200,
        shift = -100,
        additional.args = '--nomodel --keep-dup all --qval 0.05 -B --SPMR')

    # Retain only known chromosomes and remove peaks overlapping blacklisted regions
    peaks <- lapply(peaks, function(x) keepStandardChromosomes(x, pruning.mode = "coarse"))
    peaks <- lapply(peaks, function(x) subsetByOverlaps(x = x, ranges = blacklist, invert = TRUE))
    names(peaks) <- sapply(peaks, function(x) unique(x$ident))

    return(peaks)
}

young_peaks = call_peaks(young_ATAC_obj, "Young", blacklist, '/nas/homes/benyang/HiC/13_MultiOme/integrated_data/Young')
saveRDS(young_peaks, "/nas/homes/benyang/HiC/13_MultiOme/integrated_data/young_celltype_MACS2_peaks.RDS")

# MuSC feature matrices
Idents(young_ATAC_obj) = "celltype_wsnn"
young_peak.MuSC_counts = FeatureMatrix(fragments = Fragments(young_ATAC_obj), features = young_peaks[['MuSC']], cells = WhichCells(young_ATAC_obj, ident="MuSC"))
saveRDS(young_peak.MuSC_counts, "/nas/homes/benyang/HiC/13_MultiOme/integrated_data/young_ATAC.MuSC.counts.RDS")

aged_peaks = call_peaks(aged_ATAC_obj, "Aged_celltype", blacklist, '/nas/homes/benyang/HiC/13_MultiOme/integrated_data/Aged')
saveRDS(aged_peaks, "/nas/homes/benyang/HiC/13_MultiOme/integrated_data/aged_celltype_MACS2_peaks.RDS")

# MuSC feature matrices
Idents(aged_ATAC_obj) = "celltype_wsnn"
aged_peak.MuSC_counts = FeatureMatrix(fragments = Fragments(aged_ATAC_obj), features = aged_peaks[['MuSC']], cells = WhichCells(aged_ATAC_obj, ident="MuSC"))
saveRDS(aged_peak.MuSC_counts, "/nas/homes/benyang/HiC/13_MultiOme/integrated_data/aged_ATAC.MuSC.counts.RDS")


