prefix=/nas/homes/benyang/HiC
tss=$prefix/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed

bedtools slop -i $tss -b 1000 -g /nas/homes/benyang/Genome_References/sizes.mm10 | \
bedtools intersect -wo -a stdin -b $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed | \
sort -k1,1 -k2,2n | uniq | \
bedtools groupby -i stdin -grp 5-13 -c 4 -o collapse | sort -k1,1 -k2,2n > aged_genes_in_domains.txt

bedtools slop -i $tss -b 1000 -g /nas/homes/benyang/Genome_References/sizes.mm10 | \
bedtools intersect -wo -a stdin -b $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed | \
sort -k1,1 -k2,2n | uniq | \
bedtools groupby -i stdin -grp 5-13 -c 4 -o collapse | sort -k1,1 -k2,2n > young_genes_in_domains.txt