mainDir="/nas/homes/benyang/HiC"
mm10=/nas/homes/benyang/Genome_References/sizes.mm10
tss="$mainDir/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed"
outDir="$mainDir/05_Arrowhead_HICCUPS/assign_genes_to_loops"

# get gene names with boundaries overlapping promoter
bedtools slop -b 1000 -i "$tss" -g "$mm10" | \
bedtools intersect -wa \
-a stdin \
-b "$mainDir/05_Arrowhead_HICCUPS/HOMER/young_merged_loop_anchors.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/young.merged.loop_anchor.knownGenes.bed"

bedtools slop -b 1000 -i "$tss" -g "$mm10" | \
bedtools intersect -wa \
-a stdin \
-b "$mainDir/05_Arrowhead_HICCUPS/HOMER/aged_merged_loop_anchors.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/aged.merged.loop_anchor.knownGenes.bed"

bedtools slop -b 1000 -i "$tss" -g "$mm10" | \
bedtools intersect -wa \
-a stdin \
-b "$mainDir/05_Arrowhead_HICCUPS/HOMER/young_merged_diff_loop_anchors.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/young.merged.diff_loop_anchor.knownGenes.bed"

bedtools slop -b 1000 -i "$tss" -g "$mm10" | \
bedtools intersect -wa \
-a stdin \
-b "$mainDir/05_Arrowhead_HICCUPS/HOMER/aged_merged_diff_loop_anchors.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/aged.merged.diff_loop_anchor.knownGenes.bed"

# get gene names with boundaries not overlapping promoter
grep -vwf young.merged.loop_anchor.knownGenes.bed $tss > young.merged.non_loop_anchor.knownGenes.bed
grep -vwf aged.merged.loop_anchor.knownGenes.bed $tss > aged.merged.non_loop_anchor.knownGenes.bed
grep -vwf young.merged.diff_loop_anchor.knownGenes.bed $tss > young.merged.diff_non_loop_anchor.knownGenes.bed
grep -vwf aged.merged.diff_loop_anchor.knownGenes.bed $tss > aged.merged.diff_non_loop_anchor.knownGenes.bed

# bedtools slop -b 1000 -i "$tss" -g "$mm10" | \
# bedtools intersect -wa -v \
# -a stdin \
# -b "$mainDir/05_Arrowhead_HICCUPS/HOMER/young_merged_diff_loop_anchors.bed" | \
# sort -k1,1 -k2,2n | uniq > "$outDir/young.merged.non_loop_anchor.knownGenes.bed"

# bedtools slop -b 1000 -i "$tss" -g "$mm10" | \
# bedtools intersect -wa -v \
# -a stdin \
# -b "$mainDir/05_Arrowhead_HICCUPS/HOMER/aged_merged_diff_loop_anchors.bed" | \
# sort -k1,1 -k2,2n | uniq > "$outDir/aged.merged.non_loop_anchor.knownGenes.bed"