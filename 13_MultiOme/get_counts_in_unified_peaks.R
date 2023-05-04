library(future)
library(Seurat)
library(SeuratObject)
library(Signac)
library(GenomicRanges)

plan("multicore", workers = 45)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM

# Import data
aged_subset_ATAC = readRDS("/nas/homes/benyang/HiC/13_MultiOme/aged_subset_ATAC.RDS")
young_subset_ATAC = readRDS("/nas/homes/benyang/HiC/13_MultiOme/young_subset_ATAC.RDS")
young_v2_subset_ATAC = readRDS("/nas/homes/benyang/HiC/13_MultiOme/young_v2_subset_ATAC.RDS")
outdir = "/nas/homes/benyang/HiC/13_MultiOme/MACS2_peaks_by_age/"

# create blacklist object
blacklist = read.table("/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed", sep='\t')
colnames(blacklist) = c("chrom","start","end","name")
blacklist = GRanges(blacklist)

# Need to reassign fragment location since ATAC subsets were processed on local machine
young = CreateChromatinAssay(counts = GetAssayData(young_subset_ATAC, assay="ATAC", slot="counts"), fragments = "/nas/homes/benyang/HiC/13_MultiOme/RawData/Young/10x_analysis_5644-JL/Sample_5644-JL-1/atac_fragments.tsv.gz")
young_v2 = CreateChromatinAssay(counts = GetAssayData(young_v2_subset_ATAC, assay="ATAC", slot="counts"), fragments = "/nas/datasets/6098-JL/10x_analysis_6098-JL/Sample_6098-JL-1/atac_fragments.tsv.gz")
aged = CreateChromatinAssay(counts = GetAssayData(aged_subset_ATAC, assay="ATAC", slot="counts"), fragments = "/nas/homes/benyang/HiC/13_MultiOme/RawData/Aged/10x_analysis_5422-JL/Sample_5422-JL-1/atac_fragments.tsv.gz")

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

    peak.counts = FeatureMatrix(fragments = Fragments(data), features = peaks, cells = colnames(data))

    return(peak.counts)
}

young_peak.counts = call_peaks(young, "young_subset_ATAC", blacklist, outdir)
saveRDS(young_peak.counts, file.path(outdir, "young_peak_counts.RDS"))

young_v2_peak.counts = call_peaks(young_v2, "young_v2_subset_ATAC", blacklist, outdir)
saveRDS(young_v2_peak.counts, file.path(outdir, "young_v2_peak_counts.RDS"))

aged_peak.counts = call_peaks(aged, "aged_subset_ATAC", blacklist, outdir)
saveRDS(aged_peak.counts, file.path(outdir, "aged_peak_counts.RDS"))

# create a new assay using the MACS2 peak set and add it to the Seurat object
young_MACS2 <- CreateChromatinAssay(
  counts = young_peak.counts,
  fragments = "/nas/homes/benyang/HiC/13_MultiOme/RawData/Young/10x_analysis_5644-JL/Sample_5644-JL-1/atac_fragments.tsv.gz"
)

young_v2_MACS2 <- CreateChromatinAssay(
  counts = young_v2_peak.counts,
  fragments = "/nas/datasets/6098-JL/10x_analysis_6098-JL/Sample_6098-JL-1/atac_fragments.tsv.gz"
)

aged_MACS2 <- CreateChromatinAssay(
  counts = aged_peak.counts,
  fragments = "/nas/homes/benyang/HiC/13_MultiOme/RawData/Aged/10x_analysis_5422-JL/Sample_5422-JL-1/atac_fragments.tsv.gz"
)

# Identify unified set of peaks (Reduced)
aging.unified.peaks = UnifyPeaks(list(young_MACS2, young_v2_MACS2, aged_MACS2), mode="reduce")
aging.unified.peaks$width = width(aging.unified.peaks)

young_ATAC.unified_counts = FeatureMatrix(fragments = Fragments(young_MACS2), features=aging.unified.peaks, cells = colnames(young_MACS2))
young_v2_ATAC.unified_counts = FeatureMatrix(fragments = Fragments(young_v2_MACS2), features=aging.unified.peaks, cells = colnames(young_v2_MACS2))
aged_ATAC.unified_counts = FeatureMatrix(fragments = Fragments(aged_MACS2), features=aging.unified.peaks, cells = colnames(aged_MACS2))

saveRDS(young_ATAC.unified_counts, file.path(outdir, "young_unified_peak_counts.RDS"))
saveRDS(young_v2_ATAC.unified_counts, file.path(outdir, "young_v2_unified_peak_counts.RDS"))
saveRDS(aged_ATAC.unified_counts, file.path(outdir, "aged_unified_peak_counts.RDS"))