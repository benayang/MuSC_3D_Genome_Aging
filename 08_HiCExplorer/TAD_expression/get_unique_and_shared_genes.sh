mainDir="/nas/homes/benyang/HiC"
mm10=/nas/homes/benyang/Genome_References/sizes.mm10
genebody="$mainDir/get_tss/genebodies.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed"
tss="$mainDir/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed"

# get gene names with boundaries overlapping genebody
# bedtools slop -b 1000 -i "$genebody" -g "$mm10" | \
# bedtools intersect -wa \
# -a stdin \
# -b "$mainDir/08_HiCExplorer/TAD_boundary_strength/unique_young_TAD_boundary.bed" | \
# sort -k1,1 -k2,2n | uniq > unique.young.TAD.knownGenes.bed

# bedtools slop -b 1000 -i "$genebody" -g "$mm10" | \
# bedtools intersect -wa \
# -a stdin \
# -b "$mainDir/08_HiCExplorer/TAD_boundary_strength/unique_aged_TAD_boundary.bed" | \
# sort -k1,1 -k2,2n | uniq > unique.aged.TAD.knownGenes.bed

# get gene names with domains overlapping promoter
#bedtools slop -b 20000 -i shared_TAD_domain.bed -g "$mm10" > shared_TAD_domain_extended.bed
#bedtools slop -b 20000 -i unique_aged_TAD_domain.bed -g "$mm10" > unique_aged_TAD_domain_extended.bed
#bedtools slop -b 20000 -i unique_young_TAD_domain.bed -g "$mm10" > unique_young_TAD_domain_extended.bed

bedtools intersect -wa -f 1 \
-a $genebody \
-b unique_aged_TAD_domain.bed | \
sort -k1,1 -k2,2n | uniq > unique.aged.TAD.domain.knownGenes.bed

bedtools intersect -wa -f 1 \
-a $genebody \
-b unique_young_TAD_domain.bed | \
sort -k1,1 -k2,2n | uniq > unique.young.TAD.domain.knownGenes.bed

bedtools intersect -wa -f 1 \
-a $genebody \
-b shared_TAD_domain.bed | \
sort -k1,1 -k2,2n | uniq > shared.TAD.domain.knownGenes.bed

# bedtools intersect -wa \
# -a "$mainDir/08_HiCExplorer/$age.merged/40kb/$age.merged_$suffix" \
# -b "$mainDir/04_FANC/compartmentExpression/compartmentBed/$f.bed" | \
# sort -k1,1 -k2,2n | uniq > "$outDir/$f.TADDomain.bed"

bedtools slop -i $tss -g $mm10 -b 1000 | \
bedtools intersect -wb \
-a $tss \
-b $mainDir/08_HiCExplorer/TAD_boundary_strength/shared_TAD_boundary.bed | \
sort -k1,1 -k2,2n | uniq > shared.TAD.boundary.knownGenes.bed

bedtools slop -i $tss -g $mm10 -b 1000 | \
bedtools intersect -wb \
-a $tss \
-b $mainDir/08_HiCExplorer/TAD_boundary_strength/unique_young_TAD_boundary.bed | \
sort -k1,1 -k2,2n | uniq > unique.young.TAD.boundary.knownGenes.bed

bedtools slop -i $tss -g $mm10 -b 1000 | \
bedtools intersect -wb \
-a $tss \
-b $mainDir/08_HiCExplorer/TAD_boundary_strength/unique_aged_TAD_boundary.bed | \
sort -k1,1 -k2,2n | uniq > unique.aged.TAD.boundary.knownGenes.bed