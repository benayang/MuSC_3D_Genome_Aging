mainDir="/nas/homes/benyang/HiC"
mm10=/nas/homes/benyang/Genome_References/sizes.mm10
tss="$mainDir/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed"
boundarysuffix=min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed
suffix=min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed
outDir="$mainDir/08_HiCExplorer/TAD expression"

# get gene names with boundaries overlapping promoter
bedtools slop -b 1000 -i "$tss" -g "$mm10" | \
bedtools intersect -wa \
-a stdin \
-b "$mainDir/08_HiCExplorer/young.merged/40kb/young.merged_$boundarysuffix" | sort -k1,1 -k2,2n | uniq > "$outDir/young.merged.TAD.promoter.knownGenes.bed"

bedtools slop -b 1000 -i "$tss" -g "$mm10" | \
bedtools intersect -wa \
-a stdin \
-b "$mainDir/08_HiCExplorer/aged.merged/40kb/aged.merged_$boundarysuffix" | sort -k1,1 -k2,2n | uniq > "$outDir/aged.merged.TAD.promoter.knownGenes.bed"

# get boundary coordinates for promoter and non-promoter boundaries
bedtools slop -b 1000 -i "$tss" -g "$mm10" | \
bedtools intersect -wa \
-a "$mainDir/08_HiCExplorer/young.merged/40kb/young.merged_$boundarysuffix" \
-b stdin | sort -k1,1 -k2,2n | uniq > "$outDir/young.merged.TAD.promoter.boundaries.knownGenes.bed"

bedtools slop -b 1000 -i "$tss" -g "$mm10" | \
bedtools intersect -wa \
-a "$mainDir/08_HiCExplorer/aged.merged/40kb/aged.merged_$boundarysuffix" \
-b stdin | sort -k1,1 -k2,2n | uniq > "$outDir/aged.merged.TAD.promoter.boundaries.knownGenes.bed"

bedtools intersect -wa -v \
-a "$mainDir/08_HiCExplorer/young.merged/40kb/young.merged_$boundarysuffix" \
-b "$outDir/young.merged.TAD.promoter.boundaries.knownGenes.bed" | sort -k1,1 -k2,2n | uniq > "$outDir/young.merged.TAD.nonpromoter.boundaries.knownGenes.bed"

bedtools intersect -wa -v \
-a "$mainDir/08_HiCExplorer/aged.merged/40kb/aged.merged_$boundarysuffix" \
-b "$outDir/aged.merged.TAD.promoter.boundaries.knownGenes.bed" | sort -k1,1 -k2,2n | uniq > "$outDir/aged.merged.TAD.nonpromoter.boundaries.knownGenes.bed"