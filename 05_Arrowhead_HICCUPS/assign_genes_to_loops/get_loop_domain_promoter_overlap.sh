mainDir="/nas/homes/benyang/HiC"
mm10=/nas/homes/benyang/Genome_References/sizes.mm10
tss="$mainDir/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed"
outDir="$mainDir/05_Arrowhead_HICCUPS/assign_genes_to_loops"

# get gene names with boundaries overlapping promoter
bedtools slop -b 1000 -i $tss -g /nas/homes/benyang/Genome_References/sizes.mm10 | \
bedtools intersect -wa \
-a stdin \
-b "$mainDir/05_Arrowhead_HICCUPS/assign_genes_to_loops/young_merged_loop_domains.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/young.merged.loop_domain.knownGenes.bed"

bedtools slop -b 1000 -i $tss -g /nas/homes/benyang/Genome_References/sizes.mm10 | \
bedtools intersect -wa \
-a stdin \
-b "$mainDir/05_Arrowhead_HICCUPS/assign_genes_to_loops/aged_merged_loop_domains.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/aged.merged.loop_domain.knownGenes.bed"

bedtools slop -b 1000 -i $tss -g /nas/homes/benyang/Genome_References/sizes.mm10 | \
bedtools intersect -wa \
-a stdin \
-b "$mainDir/05_Arrowhead_HICCUPS/assign_genes_to_loops/young_diff_loop_domains.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/young.merged.diff_loop_domain.knownGenes.bed"

bedtools slop -b 1000 -i $tss -g /nas/homes/benyang/Genome_References/sizes.mm10 | \
bedtools intersect -wa \
-a stdin \
-b "$mainDir/05_Arrowhead_HICCUPS/assign_genes_to_loops/aged_diff_loop_domains.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/aged.merged.diff_loop_domain.knownGenes.bed"