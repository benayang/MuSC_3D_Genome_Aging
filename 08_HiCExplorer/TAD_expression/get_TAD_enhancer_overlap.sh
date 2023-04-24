mainDir="/mnt/c/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C"
mm10="$mainDir/sizes.mm10"
atac="$mainDir/12_ATAC"
exons="$mainDir/get_tss/genebodies.gencode.vM25.basic.annotation.filtered.uniq.bed"
enhancers="$mainDir/08_HiCExplorer/TAD expression/F5.mm10.enhancers.sorted.bed"
suffix="40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed"
boundarysuffix="40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed"
outDir="$mainDir/08_HiCExplorer/TAD expression"

bedtools intersect -wb \
-a "$enhancers" \
-b "$mainDir/08_HiCExplorer/young.merged/40kb/young.merged_$boundarysuffix" | \
bedtools intersect -v -wb \
-a "$exons" \
-b stdin > "$outDir/young.merged.TAD.enhancers.bed"

bedtools intersect -wb \
-a "$enhancers" \
-b "$mainDir/08_HiCExplorer/aged.merged/40kb/aged.merged_$boundarysuffix" | \
bedtools intersect -v -wb \
-a "$exons" \
-b stdin > "$outDir/aged.merged.TAD.enhancers.bed"

cat "$outDir/young.merged.TAD.enhancers.bed" | wc -l
cat "$outDir/aged.merged.TAD.enhancers.bed" | wc -l