prefix=/nas/homes/benyang/HiC

outDir="$prefix/08_HiCExplorer/TAD expression/promoter_distance"
tss=$prefix/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed
boundarysuffix=min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed

cat $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_$boundarysuffix | \
sort -k1,1 -k2,2n | uniq | \
bedtools intersect -wa -a stdin -b "$prefix/04_FANC/compartmentExpression/compartmentBed/aged.A.bed" | \
bedtools closest -d -a $tss -b stdin > "$outDir/aged.promoter.A.boundaries.knownGenes.distance.bed"

cat $prefix/08_HiCExplorer/young.merged/40kb/young.merged_$boundarysuffix | \
sort -k1,1 -k2,2n | uniq | \
bedtools intersect -wa -a stdin -b "$prefix/04_FANC/compartmentExpression/compartmentBed/young.A.bed" | \
bedtools closest -d -a $tss -b stdin > "$outDir/young.promoter.A.boundaries.knownGenes.distance.bed"

cat $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_$boundarysuffix | \
sort -k1,1 -k2,2n | uniq | \
bedtools intersect -wa -a stdin -b "$prefix/04_FANC/compartmentExpression/compartmentBed/aged.B.bed" | \
bedtools closest -d -a $tss -b stdin > "$outDir/aged.promoter.B.boundaries.knownGenes.distance.bed"

cat $prefix/08_HiCExplorer/young.merged/40kb/young.merged_$boundarysuffix | \
sort -k1,1 -k2,2n | uniq | \
bedtools intersect -wa -a stdin -b "$prefix/04_FANC/compartmentExpression/compartmentBed/young.B.bed" | \
bedtools closest -d -a $tss -b stdin > "$outDir/young.promoter.B.boundaries.knownGenes.distance.bed"
