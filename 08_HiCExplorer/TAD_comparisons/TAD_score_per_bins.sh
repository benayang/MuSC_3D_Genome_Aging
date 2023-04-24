multiBigwigSummary bins \
-b "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bw" \
"/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bw" \
--binSize 40000 \
-bl "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" \
-p 45 \
--labels Young Aged \
--outFileName TAD_scores_per_bin.npz --outRawCounts TAD_scores_per_bin.tab