mainDir="/mnt/c/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C"
mm10="$mainDir/sizes.mm10"
genebodies="$mainDir/get_tss/genebodies.gencode.vM25.basic.annotation.filtered.uniq.bed"
suffix="40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed"
boundarysuffix="40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed"
outDir="$mainDir/08_HiCExplorer/TAD expression"

bedtools intersect -wa \
-a "$genebodies" \
-b "$mainDir/08_HiCExplorer/young.merged/40kb/young.merged_$suffix" > "$outDir/young.merged.TAD.genebodies.bed"

bedtools intersect -wa \
-a "$genebodies" \
-b "$mainDir/08_HiCExplorer/aged.merged/40kb/aged.merged_$suffix" > "$outDir/aged.merged.TAD.genebodies.bed"

