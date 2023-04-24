prefix=/nas/homes/benyang/HiC

for age in young aged
do

cat $prefix/08_HiCExplorer/$age.merged/40kb/$age.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed | \
bedtools slop -b 20000 -i stdin -g /nas/homes/benyang/Genome_References/sizes.mm10 | \
bedtools intersect -u -wa -f 0.8 \
-a $prefix/05_Arrowhead_HICCUPS/assign_genes_to_loops/${age}_merged_loop_domains.bed \
-b stdin | \
sort -k1,1 -k2,2n > ${age}_loops_in_TADs.bed

grep -wvf ${age}_loops_in_TADs.bed $prefix/05_Arrowhead_HICCUPS/assign_genes_to_loops/${age}_merged_loop_domains.bed | \
sort -k1,1 -k2,2n > ${age}_loops_outside_TADs.bed

bedtools intersect -u -wa \
-a $prefix/08_HiCExplorer/$age.merged/40kb/$age.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
-b $prefix/05_Arrowhead_HICCUPS/HOMER/${age}_merged_loop_anchors.bed | \
sort -k1,1 -k2,2n > ${age}_loops_in_TAD_boundaries.bed

grep -wvf ${age}_loops_in_TAD_boundaries.bed $prefix/08_HiCExplorer/$age.merged/40kb/$age.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed | \
sort -k1,1 -k2,2n > ${age}_loops_outside_TAD_boundaries.bed

done