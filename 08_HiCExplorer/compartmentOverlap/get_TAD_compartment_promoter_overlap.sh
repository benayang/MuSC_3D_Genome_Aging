mainDir="/nas/homes/benyang/HiC"
mm10="/nas/homes/benyang/Genome_References/sizes.mm10"
genesDir=$mainDir/04_FANC/compartmentExpression
tss="$mainDir/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed"
outDir="$mainDir/08_HiCExplorer/compartmentOverlap/promoterOverlap"

bedtools slop -b 1000 -i "$tss" -g "$mm10" > tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.1kb.bed

# get gene names with boundaries overlapping promoter
for f in "A_to_B" "B_to_A" "StaticA" "StaticB"
do
    for a in young aged
    do

    bedtools intersect -wa \
    -a $genesDir/$f.1kb.tss.bed \
    -b $mainDir/08_HiCExplorer/$a.merged/40kb/$a.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed | \
    sort -k1,1 -k2,2n | uniq > "$outDir/$a.$f.TADBoundary.promoter.knownGenes.bed"

    #grep -vwf "$outDir/$a.$f.TADBoundary.promoter.knownGenes.bed" $genesDir/$f.1kb.tss.bed > "$outDir/$a.$f.non_TADBoundary.promoter.knownGenes.bed"

    bedtools intersect -wa \
    -a $genesDir/$f.1kb.tss.bed \
    -b $mainDir/08_HiCExplorer/$a.merged/40kb/$a.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed | \
    sort -k1,1 -k2,2n | uniq > "$outDir/$a.$f.TADDomain.promoter.knownGenes.bed"

    #grep -vwf "$outDir/$a.$f.TADDomain.promoter.knownGenes.bed" $genesDir/$f.1kb.tss.bed > "$outDir/$a.$f.non_TADDomain.promoter.knownGenes.bed"

    done
done

for f in "A" "B"
do
    for a in young aged
    do

    bedtools intersect -wa \
    -a $genesDir/$a.$f.1kb.tss.bed \
    -b $mainDir/08_HiCExplorer/$a.merged/40kb/$a.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed | \
    sort -k1,1 -k2,2n | uniq > "$outDir/$a.$f.TADBoundary.promoter.knownGenes.bed"

    #grep -vwf "$outDir/$a.$f.TADBoundary.promoter.knownGenes.bed" $genesDir/$a.$f.1kb.tss.bed > "$outDir/$a.$f.non_TADBoundary.promoter.knownGenes.bed"

    bedtools intersect -wa \
    -a $genesDir/$a.$f.1kb.tss.bed \
    -b $mainDir/08_HiCExplorer/$a.merged/40kb/$a.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed | \
    sort -k1,1 -k2,2n | uniq > "$outDir/$a.$f.TADDomain.promoter.knownGenes.bed"

    #grep -vwf "$outDir/$a.$f.TADDomain.promoter.knownGenes.bed" $genesDir/$a.$f.1kb.tss.bed > "$outDir/$a.$f.non_TADDomain.promoter.knownGenes.bed"
    
    done
done


# get gene names with domains overlapping promoter
# for f in young.A.TADDomain young.B.TADDomain aged.A.TADDomain aged.B.TADDomain
# do

# bedtools intersect -wa \
# -a tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.1kb.bed \
# -b "$outDir/$f.bed" > \
# "$outDir/$f.promoter.knownGenes.bed"

# done