mainDir="/nas/homes/benyang/HiC"
mm10=/nas/homes/benyang/Genome_References/sizes.mm10
tss="$mainDir/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed"
boundarysuffix=min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed
outDir="$mainDir/08_HiCExplorer/kmeans_TAD_boundaries"

# get gene names with boundaries overlapping promoter
for c in `seq 1 1 6`
do

bedtools slop -b 1000 -i "$tss" -g "$mm10" | \
bedtools intersect -wa \
-a stdin \
-b $mainDir/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster$c.bedgraph | sort -k1,1 -k2,2n | uniq > $outDir/young_cluster$c.promoter.knownGenes.bed

bedtools slop -b 1000 -i "$tss" -g "$mm10" | \
bedtools intersect -wa \
-a stdin \
-b $mainDir/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster$c.bedgraph | sort -k1,1 -k2,2n | uniq > $outDir/aged_cluster$c.promoter.knownGenes.bed

done