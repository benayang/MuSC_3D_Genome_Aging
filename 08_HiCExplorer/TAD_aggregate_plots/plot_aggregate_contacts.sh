prefix=/nas/homes/benyang/HiC

# hicAggregateContacts \
# --matrix $prefix/08_HiCExplorer/aged.merged_10kb_KR.cool \
# --outFileName aged_merged_10kb_CTCF_median_intraChr \
# --BED $prefix/HistoneTracks/C2C12_CTCF_ENCFF784ASD.sorted.bed \
# --mode intra-chr \
# --range 300000:1000000 --numberOfBins 30


hicAggregateContacts \
--matrix $prefix/08_HiCExplorer/aged.merged_40kb_KR.cool \
--outFileName aged_merged_40kb_TAD_boundary_intraChr \
--BED $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
--mode intra-chr \
--range 1000000:2000000 --numberOfBins 60