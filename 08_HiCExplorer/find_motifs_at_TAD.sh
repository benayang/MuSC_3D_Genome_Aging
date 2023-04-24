prefix=/nas/homes/benyang/HiC

# TAD boundaries
bedtools intersect -wo -a $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
-b $prefix/08_HiCExplorer/homer/ctcf.scan.out.bed > aged_40kb_TAD_boundary_CTCF.bed

bedtools intersect -wo -a $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
-b $prefix/08_HiCExplorer/homer/ctcf.scan.out.bed > young_40kb_TAD_boundary_CTCF.bed

# Whole TAD including boundaries
cat $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed | \
bedtools slop -i stdin -g /nas/homes/benyang/Genome_References/sizes.mm10 -b 20000 | \
bedtools makewindows -b stdin -w 20000 -i srcwinnum | \
bedtools annotate -i stdin -counts -files $prefix/08_HiCExplorer/homer/ctcf.scan.out.pos.bed $prefix/08_HiCExplorer/homer/ctcf.scan.out.neg.bed -names pos neg > aged_40kb_TAD_CTCF.bed
# bedtools intersect -wao -a stdin \
# -b $prefix/08_HiCExplorer/homer/ctcf.scan.out.bed > aged_40kb_TAD_CTCF.bed

cat $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed | \
bedtools slop -i stdin -g /nas/homes/benyang/Genome_References/sizes.mm10 -b 20000 | \
bedtools makewindows -b stdin -w 20000 -i srcwinnum | \
bedtools annotate -i stdin -counts -files $prefix/08_HiCExplorer/homer/ctcf.scan.out.pos.bed $prefix/08_HiCExplorer/homer/ctcf.scan.out.neg.bed -names pos neg > young_40kb_TAD_CTCF.bed
# bedtools intersect -wao -a stdin \
# -b $prefix/08_HiCExplorer/homer/ctcf.scan.out.bed > young_40kb_TAD_CTCF.bed
