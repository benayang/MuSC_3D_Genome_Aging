prefix=/nas/homes/benyang/HiC/08_HiCExplorer


bedtools intersect -c -a $prefix/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
-b $prefix/homer/ctcf.scan.out.bed > aged_boundaries_CTCF.bed
bedtools intersect -c -a $prefix/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
-b $prefix/homer/ctcf.scan.out.bed > aged_domains_CTCF.bed

bedtools intersect -c -a $prefix/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
-b $prefix/homer/ctcf.scan.out.bed > young_boundaries_CTCF.bed
bedtools intersect -c -a $prefix/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
-b $prefix/homer/ctcf.scan.out.bed > young_domains_CTCF.bed
