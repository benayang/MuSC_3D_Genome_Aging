prefix=/nas/homes/benyang/HiC

# hicAverageRegions \
# --matrix $prefix/08_HiCExplorer/aged.merged_40kb_KR.cool \
# --regions $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
# --rangeInBins 10 10 \
# --outFileName aged_merged_TAD_domains

hicPlotAverageRegions --dpi 300 --log1p -m aged_merged_TAD_domains.npz -o aged_merged_TAD_domains.png

# hicAverageRegions \
# --matrix $prefix/08_HiCExplorer/young.merged_40kb_KR.cool \
# --regions $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
# --rangeInBins 10 10 \
# --outFileName young_merged_TAD_domains

hicPlotAverageRegions --dpi 300 --log1p -m young_merged_TAD_domains.npz -o young_merged_TAD_domains.png
