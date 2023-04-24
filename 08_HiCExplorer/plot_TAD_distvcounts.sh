prefix=/nas/homes/benyang/HiC/08_HiCExplorer

hicPlotDistVsCounts -m $prefix/young.merged_40kb_KR.cool \
--plotFile young_TAD_dist_vs_counts.png \
--labels Young \
--maxdepth 3000000 \
--domains $prefix/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
--outFileData young_TAD_dist_vs_counts.txt \
--plotsize 5 5

hicPlotDistVsCounts -m $prefix/aged.merged_40kb_KR.cool \
--plotFile aged_TAD_dist_vs_counts.png \
--labels Aged \
--maxdepth 3000000 \
--domains $prefix/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
--outFileData aged_TAD_dist_vs_counts.txt \
--plotsize 5 5
