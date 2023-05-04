library(Seurat)
library(Signac)

# Young
young_MuSC = readRDS("/nas/homes/benyang/HiC/13_MultiOme/young_MuSC.RDS")
young_fragments = Fragments(young_MuSC)[[1]]
young_fragments = UpdatePath(young_fragments, new.path="/nas/homes/benyang/HiC/13_MultiOme/RawData/Young/10x_analysis_5644-JL/Sample_5644-JL-1/atac_fragments.tsv.gz")

Fragments(young_MuSC) = NULL
Fragments(young_MuSC) = young_fragments

CoveragePlot(young_MuSC, region="chr7-46033475-46710819", features="Myod1", assay="ATAC", expression.assay="RNA", extend.upstream=1e4, extend.downstream=1e4, peaks=T, bigwig=list("H3K4me3"="/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw", "ATAC"="/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw"), bigwig.scale="separate")
ggsave(file.path("/nas/homes/benyang/HiC/13_MultiOme/integrated_data/young_Myod1_Coverage_bigwig.png"), dpi=300, width=8, height=5)


# Aged
aged_MuSC = readRDS("/nas/homes/benyang/HiC/13_MultiOme/aged_MuSC.RDS")
aged_fragments = Fragments(aged_MuSC)[[1]]
aged_fragments = UpdatePath(aged_fragments, new.path="/nas/homes/benyang/HiC/13_MultiOme/RawData/Aged/10x_analysis_5422-JL/Sample_5422-JL-1/atac_fragments.tsv.gz")

Fragments(aged_MuSC) = NULL
Fragments(aged_MuSC) = aged_fragments

CoveragePlot(aged_MuSC, region="chr7-46033475-46710819", features="Myod1", assay="ATAC", expression.assay="RNA", extend.upstream=1e4, extend.downstream=1e4, peaks=T, bigwig=list("H3K4me3"="/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw", "ATAC"="/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw"), bigwig.scale="separate")
ggsave(file.path("/nas/homes/benyang/HiC/13_MultiOme/integrated_data/aged_Myod1_Coverage_bigwig.png"), dpi=300, width=8, height=5)
