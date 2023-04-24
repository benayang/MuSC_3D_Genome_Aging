prefix=/nas/homes/benyang/HiC

# Get CTCF distribution across loop domains
cat $prefix/05_Arrowhead_HICCUPS/assign_genes_to_loops/aged_merged_loop_domains.bed | \
awk '{OFS="\t"}{print $s,NR}' | \
bedtools makewindows -b stdin -w 5000 -i srcwinnum | \
bedtools annotate -i stdin -counts -files $prefix/08_HiCExplorer/homer/ctcf.scan.out.pos.bed $prefix/08_HiCExplorer/homer/ctcf.scan.out.neg.bed -names pos neg | \
sort -k1,1 -k2,2n > aged_loop_domains_CTCF.bed

cat $prefix/05_Arrowhead_HICCUPS/assign_genes_to_loops/young_merged_loop_domains.bed | \
awk '{OFS="\t"}{print $s,NR}' | \
bedtools makewindows -b stdin -w 5000 -i srcwinnum | \
bedtools annotate -i stdin -counts -files $prefix/08_HiCExplorer/homer/ctcf.scan.out.pos.bed $prefix/08_HiCExplorer/homer/ctcf.scan.out.neg.bed -names pos neg | \
sort -k1,1 -k2,2n > young_loop_domains_CTCF.bed