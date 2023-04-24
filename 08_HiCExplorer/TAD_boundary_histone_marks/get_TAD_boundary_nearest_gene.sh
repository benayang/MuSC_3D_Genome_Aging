prefix=/nas/homes/benyang/HiC

bedtools sort -i $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed | bedtools closest -d -a $prefix/get_tss/tss_1kb_slop.bed -b stdin > genes_closest_to_Young_TAD_boundaries.txt

bedtools sort -i $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed | bedtools closest -d -a $prefix/get_tss/tss_1kb_slop.bed -b stdin > genes_closest_to_Aged_TAD_boundaries.txt