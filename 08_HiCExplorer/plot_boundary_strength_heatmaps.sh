prefix=/nas/homes/benyang/HiC/08_HiCExplorer

# LC_COLLATE=C sort -k1,1 -k2,2n $prefix/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bedgraph > \
# aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.sorted.bedgraph

# LC_COLLATE=C sort -k1,1 -k2,2n $prefix/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bedgraph > \
# young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.sorted.bedgraph

# bedGraphToBigWig aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.sorted.bedgraph \
# /nas/homes/benyang/Genome_References/sizes.mm10 \
# $prefix/TAD_boundary_strength/aged_TAD_separation_score.bw

# bedGraphToBigWig young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.sorted.bedgraph \
# /nas/homes/benyang/Genome_References/sizes.mm10 \
# $prefix/TAD_boundary_strength/young_TAD_separation_score.bw

######### 40kb flanking region
# computeMatrix reference-point \
# -R $prefix/TAD_boundary_strength/shared_TAD_boundary.bed \
# $prefix/TAD_boundary_strength/unique_young_TAD_boundary.bed \
# $prefix/TAD_boundary_strength/unique_aged_TAD_boundary.bed \
# -S $prefix/TAD_boundary_strength/young_TAD_separation_score.bw \
# $prefix/TAD_boundary_strength/aged_TAD_separation_score.bw \
# --outFileName TAD_boundary_strength.mat.gz \
# --referencePoint center \
# --beforeRegionStartLength 40000 \
# --afterRegionStartLength 40000 \
# -p 55 \
# --binSize 10 \
# --samplesLabel Young Aged

# plotHeatmap -m TAD_boundary_strength.mat.gz \
# --dpi 300 \
# --outFileName TAD_boundary_strength.png

######### 500kb flanking region
computeMatrix reference-point \
-R $prefix/TAD_boundary_strength/shared_TAD_boundary.bed \
$prefix/TAD_boundary_strength/unique_young_TAD_boundary.bed \
$prefix/TAD_boundary_strength/unique_aged_TAD_boundary.bed \
-S $prefix/TAD_boundary_strength/young_TAD_separation_score.bw \
$prefix/TAD_boundary_strength/aged_TAD_separation_score.bw \
--outFileName TAD_boundary_strength_500kb_flank.mat.gz \
--referencePoint center \
--beforeRegionStartLength 500000 \
--afterRegionStartLength 500000 \
-p 55 \
--binSize 10 \
--samplesLabel Young Aged

plotHeatmap -m TAD_boundary_strength_500kb_flank.mat.gz \
--dpi 300 \
--regionsLabel Shared Young Aged \
--outFileName TAD_boundary_strength_500kb_flank.png \
--legendLocation none \
--xAxisLabel "distance (bp)"