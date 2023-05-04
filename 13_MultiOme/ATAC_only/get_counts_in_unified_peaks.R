library(future)
library(Seurat)
library(SeuratObject)
library(Signac)
library(GenomicRanges)

plan("multicore", workers = 45)
options(future.globals.maxSize = 75 * 1024 ^ 3) # for 50 Gb RAM

# Import data
aged_subset_ATAC = readRDS("/nas/homes/benyang/HiC/13_MultiOme/ATAC_only/aged_subset_ATAC.RDS")
young_subset_ATAC = readRDS("/nas/homes/benyang/HiC/13_MultiOme/ATAC_only/young_subset_ATAC.RDS")
young_v2_subset_ATAC = readRDS("/nas/homes/benyang/HiC/13_MultiOme/ATAC_only/young_v2_subset_ATAC.RDS")
outdir = "/nas/homes/benyang/HiC/13_MultiOme/ATAC_only/unified_peaks/"

# create blacklist object
blacklist = read.table("/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed", sep='\t')
colnames(blacklist) = c("chrom","start","end","name")
blacklist = GRanges(blacklist)

update_frag_path = function(obj, new.path) {
  frags = Fragments(obj)[[1]]
  Fragments(obj) = NULL
  frags = UpdatePath(frags, new.path=new.path)
  Fragments(obj) = frags
  return(obj)
}

# Need to reassign fragment location since ATAC subsets were processed on local machine
young_subset_ATAC = update_frag_path(young_subset_ATAC, "/nas/homes/benyang/HiC/13_MultiOme/RawData/Young/10x_analysis_5644-JL/Sample_5644-JL-1/atac_fragments.tsv.gz")
young_v2_subset_ATAC = update_frag_path(young_v2_subset_ATAC, "/nas/datasets/6098-JL/10x_analysis_6098-JL/Sample_6098-JL-1/atac_fragments.tsv.gz")
aged_subset_ATAC = update_frag_path(aged_subset_ATAC, "/nas/homes/benyang/HiC/13_MultiOme/RawData/Aged/10x_analysis_5422-JL/Sample_5422-JL-1/atac_fragments.tsv.gz")

# Identify unified set of peaks (Reduced)
aging.unified.peaks = UnifyPeaks(list(young_subset_ATAC, young_v2_subset_ATAC, aged_subset_ATAC), mode="reduce")
aging.unified.peaks$width = width(aging.unified.peaks)

young_ATAC.unified_counts = FeatureMatrix(fragments = Fragments(young_subset_ATAC), features=aging.unified.peaks, cells = colnames(young_subset_ATAC))
young_v2_ATAC.unified_counts = FeatureMatrix(fragments = Fragments(young_v2_subset_ATAC), features=aging.unified.peaks, cells = colnames(young_v2_subset_ATAC))
aged_ATAC.unified_counts = FeatureMatrix(fragments = Fragments(aged_subset_ATAC), features=aging.unified.peaks, cells = colnames(aged_subset_ATAC))

saveRDS(young_ATAC.unified_counts, file.path(outdir, "young_unified_peak_counts.RDS"))
saveRDS(young_v2_ATAC.unified_counts, file.path(outdir, "young_v2_unified_peak_counts.RDS"))
saveRDS(aged_ATAC.unified_counts, file.path(outdir, "aged_unified_peak_counts.RDS"))