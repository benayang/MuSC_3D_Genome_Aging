bedtools unionbedg -i young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bedgraph \
aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bedgraph \
-g /nas/homes/benyang/Genome_References/sizes.mm10 \
-empty -filler 'NA' > union_young_aged_40kb_TAD_score.txt
