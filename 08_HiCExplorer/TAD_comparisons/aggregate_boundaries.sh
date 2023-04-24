cat "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed" \
"/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed" | \
cut -f1-4 | sort -k1,1 -k2,2n | uniq > ./aggregate_TAD_boundaries.bed
