library(EnsDb.Mmusculus.v79)
library(Seurat)
library(Signac)
library(future)
library(dplyr)
library(tidyr)

create_seurat_objects = function(indir, age, annotations, md_path, max.lines) {
  raw.data = Read10X(file.path(indir, "raw_feature_bc_matrix"))
  # extract RNA and ATAC data
  rna_counts <- raw.data$`Gene Expression`
  atac_counts <- raw.data$Peaks
  
  tmp <- CreateSeuratObject(counts = raw.data$`Gene Expression`, project = age)
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^mt-")
  
  # Now add in the ATAC-seq data
  # we'll only use peaks in standard chromosomes
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  
  frag.file <- file.path(indir,"atac_fragments.tsv.gz")
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = frag.file,
    max.lines = max.lines,
    annotation = annotations
  )
  
  tmp[["ATAC"]] <- chrom_assay
  
  DefaultAssay(tmp) = "ATAC"
  tmp = TSSEnrichment(tmp, fast=F)  
  tmp = NucleosomeSignal(tmp)

  metadata = read.table(md_path, sep=',', header=T)
  metadata = metadata %>% mutate(FRiP = atac_peak_region_fragments / atac_fragments)
  tmp$FRiP = metadata[match(colnames(tmp), metadata$barcode), "FRiP"]
  tmp$blacklist_fraction = FractionCountsInRegion(tmp, assay="ATAC", regions=blacklist_mm10)

  return(tmp)
}

plan("multicore", workers = 45)
options(future.globals.maxSize = 75 * 1024 ^ 3) # for 50 Gb RAM

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

aged = create_seurat_objects("/nas/homes/benyang/HiC/13_MultiOme/RawData/Aged/10x_analysis_5422-JL/Sample_5422-JL-1","Aged", annotations,"/nas/homes/benyang/HiC/13_MultiOme/RawData/Aged/10x_analysis_5422-JL/Sample_5422-JL-1/per_barcode_metrics.csv",5e7)
saveRDS(aged, "/nas/homes/benyang/HiC/13_MultiOme/RawData/aged_raw_seurat.RDS")
write.table(aged@meta.data, "/nas/homes/benyang/HiC/13_MultiOme/RawData/barcode_rank_plots/aged_raw_metadata.tsv", sep='\t', quote=F)

young = create_seurat_objects("/nas/homes/benyang/HiC/13_MultiOme/RawData/Young/10x_analysis_5644-JL/Sample_5644-JL-1","Young", annotations,"/nas/homes/benyang/HiC/13_MultiOme/RawData/Young/10x_analysis_5644-JL/Sample_5644-JL-1/per_barcode_metrics.csv",NULL)
saveRDS(young, "/nas/homes/benyang/HiC/13_MultiOme/RawData/young_raw_seurat.RDS")
write.table(young@meta.data, "/nas/homes/benyang/HiC/13_MultiOme/RawData/barcode_rank_plots/young_raw_metadata.tsv", sep='\t', quote=F)

young_v2 = create_seurat_objects("/nas/datasets/6098-JL/10x_analysis_6098-JL/Sample_6098-JL-1","Young_v2", annotations, "/nas/datasets/6098-JL/10x_analysis_6098-JL/Sample_6098-JL-1/per_barcode_metrics.csv", 5e7)
saveRDS(young_v2, "/nas/homes/benyang/HiC/13_MultiOme/RawData/young_v2_raw_seurat.RDS")
write.table(young_v2@meta.data, "/nas/homes/benyang/HiC/13_MultiOme/RawData/barcode_rank_plots/young_v2_raw_metadata.tsv", sep='\t', quote=F)