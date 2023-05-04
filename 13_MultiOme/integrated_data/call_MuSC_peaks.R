library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

plan("multicore", workers = 45)
options(future.globals.maxSize = 70 * 1024 ^ 3) # for 70 Gb RAM

data = readRDS("/nas/homes/benyang/HiC/13_MultiOme/all_reps_aging_integrated.RDS")
DefaultAssay(data) = "ATAC"

# create blacklist object
blacklist = read.table("/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed", sep='\t')
colnames(blacklist) = c("chrom","start","end","name")
blacklist = GRanges(blacklist)

# split integrated data into Young and Aged
# Need to add a new metadata column to account for pooled Young
data$Age = sapply(data$orig.ident, function(x) ifelse(x=="Aged", "Aged", "Young"))
age_split_data = SplitObject(data, split.by="Age")

young.counts = age_split_data$Young@assays$ATAC@counts
aged.counts = age_split_data$Aged@assays$ATAC@counts

sapply(colnames(young.counts), function(x) gsub("_1","",x,fixed=T))

Idents(data) = "orig.ident"
young.frag = CreateFragmentObject("/nas/homes/benyang/HiC/13_MultiOme/RawData/Young/10x_analysis_5644-JL/Sample_5644-JL-1/atac_fragments.tsv.gz", cells=sapply(WhichCells(data, idents="Young"), function(x) gsub("_1","",x,fixed=T)), validate.fragments=T)
young_v2.frag = CreateFragmentObject("/nas/datasets/6098-JL/10x_analysis_6098-JL/Sample_6098-JL-1/atac_fragments.tsv.gz", cells=sapply(WhichCells(data, idents="Young_v2"), function(x) gsub("_2","",x,fixed=T)), validate.fragments=T)
aged.frag = CreateFragmentObject("/nas/homes/benyang/HiC/13_MultiOme/RawData/Aged/10x_analysis_5422-JL/Sample_5422-JL-1/atac_fragments.tsv.gz", cells=sapply(WhichCells(data, idents="Aged"), function(x) gsub("_3","",x,fixed=T)), validate.fragments=T)

young_assay = CreateChromatinAssay(counts = young.counts, fragments = list(young.frag, young_v2.frag))
aged_assay = CreateChromatinAssay(counts = aged.counts, fragments = aged.frag)
# young_assay = CreateChromatinAssay(counts = young.counts, fragments = "/nas/homes/benyang/HiC/13_MultiOme/RawData/Young/10x_analysis_5644-JL/Sample_5644-JL-1/atac_fragments.tsv.gz")
# young_v2_assay = CreateChromatinAssay(counts = young.counts, fragments = "/nas/datasets/6098-JL/10x_analysis_6098-JL/Sample_6098-JL-1/atac_fragments.tsv.gz")
# aged_assay = CreateChromatinAssay(counts = aged.counts, fragments = "/nas/homes/benyang/HiC/13_MultiOme/RawData/Aged/10x_analysis_5422-JL/Sample_5422-JL-1/atac_fragments.tsv.gz")

# young_md = age_split_data$Young@meta.data
# rownames(young_md) = sapply(rownames(young_md), function(x) gsub("_1","",x,fixed=T))
# young_v2_md = age_split_data$Young_v2@meta.data
# rownames(young_v2_md) = sapply(rownames(young_v2_md), function(x) gsub("_2","",x,fixed=T))
# aged_md = age_split_data$Aged@meta.data
# rownames(aged_md) = sapply(rownames(aged_md), function(x) gsub("_3","",x,fixed=T))

# young_ATAC_obj = CreateSeuratObject(counts = young_assay, assay="ATAC", project="Young", meta.data=young_md)
# young_v2_ATAC_obj = CreateSeuratObject(counts = young_v2_assay, assay="ATAC", project="Young_v2", meta.data=young_v2_md)
# aged_ATAC_obj = CreateSeuratObject(counts = aged_assay, assay="ATAC", project="Aged", meta.data=aged_md)
young_ATAC_obj = CreateSeuratObject(counts = young_assay, assay="ATAC", project="Young", meta.data=age_split_data$Young@meta.data)
aged_ATAC_obj = CreateSeuratObject(counts = aged_assay, assay="ATAC", project="Aged", meta.data=age_split_data$Aged@meta.data)


# Call peaks to replace 10X peaks
call_peaks = function(data, name, blacklist, outdir) {
    peaks = CallPeaks(
        data, 
        macs2.path="/opt/miniconda3/envs/py38/bin/macs2", 
        name=name,
        combine.peaks = FALSE,
        effective.genome.size=1.87e9, 
        outdir=outdir, 
        cleanup=F, 
        format = "BED",
        extsize = 200,
        shift = -100,
        additional.args = '--nomodel --keep-dup all --qval 0.05 -B --SPMR')

    # Retain only known chromosomes and remove peaks overlapping blacklisted regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist, invert = TRUE)

    return(peaks)
}

young_peaks = call_peaks(young_ATAC_obj, "Young", blacklist, '/nas/homes/benyang/HiC/13_MultiOme/integrated_data/all_reps_Young_MuSC')
saveRDS(young_peaks, "/nas/homes/benyang/HiC/13_MultiOme/integrated_data/all_reps_young_MuSC_MACS2_peaks.RDS")

# MuSC feature matrices
Idents(young_ATAC_obj) = "wsnn_res.0.3"
young_peak.MuSC_counts = FeatureMatrix(fragments = Fragments(young_ATAC_obj), features = young_peaks, cells=WhichCells(young_ATAC_obj, idents="9"))
saveRDS(young_peak.MuSC_counts, "/nas/homes/benyang/HiC/13_MultiOme/integrated_data/all_reps_young_ATAC.MuSC.counts.RDS")

aged_peaks = call_peaks(aged_ATAC_obj, "Aged", blacklist, '/nas/homes/benyang/HiC/13_MultiOme/integrated_data/all_reps_Aged_MuSC')
saveRDS(aged_peaks, "/nas/homes/benyang/HiC/13_MultiOme/integrated_data/all_reps_aged_MuSC_MACS2_peaks.RDS")

# MuSC feature matrices
Idents(aged_ATAC_obj) = "wsnn_res.0.3"
aged_peak.MuSC_counts = FeatureMatrix(fragments = Fragments(aged_ATAC_obj), features = aged_peaks, cells=WhichCells(aged_ATAC_obj, idents="9"))
saveRDS(aged_peak.MuSC_counts, "/nas/homes/benyang/HiC/13_MultiOme/integrated_data/all_reps_aged_ATAC.MuSC.counts.RDS")
