prefix=/nas/homes/benyang/HiC

#Young
# computeMatrix reference-point \
# -R $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
# -S $prefix/HistoneTracks/get_signal_bigwigs/H3K27me3.count.rpkm.track.nodup.Young.bw \
# $prefix/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw \
# --outFileName young_H3K27me3_H3K4me3.mat.gz \
# --referencePoint center -b 80000 -a 80000 \
# --binSize 100 \
# -p 55

# plotHeatmap -m young_H3K27me3_H3K4me3.mat.gz \
# --outFileName young_H3K27me3_H3K4me3_kmeans4.png \
# --kmeans 3 \
# --outFileSortedRegions young_H3K27me3_H3K4me3_kmeans4_sortedRegions.txt \
# --outFileNameMatrix young_H3K27me3_H3K4me3_kmeans4_matrix.txt \
# --xAxisLabel "distance (bp)" \
# --samplesLabel H3K27me3 H3K4me3 \
# --dpi 300

computeMatrix scale-regions \
-R $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
-S $prefix/HistoneTracks/get_signal_bigwigs/H3K27me3.count.rpkm.track.nodup.Young.bw \
$prefix/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw \
--outFileName young_H3K27me3_H3K4me3_scale.mat.gz \
-b 500 -a 500 \
--regionBodyLength 1000 \
--binSize 5 \
-p 55

plotHeatmap -m young_H3K27me3_H3K4me3_scale.mat.gz \
--outFileName young_H3K27me3_H3K4me3_scale_kmeans4.png \
--kmeans 3 \
--outFileSortedRegions young_H3K27me3_H3K4me3_scale_kmeans4_sortedRegions.txt \
--outFileNameMatrix young_H3K27me3_H3K4me3_scale_kmeans4_matrix.txt \
--xAxisLabel "distance (bp)" \
--samplesLabel H3K27me3 H3K4me3 \
--dpi 300

# Aged
# computeMatrix reference-point \
# -R $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
# -S $prefix/04_FANC/without_KR_normalization/aged.merged_100kb_noKR.bw \
# --outFileName aged_AB.mat.gz \
# --referencePoint center -b 500000 -a 500000 \
# --binSize 10000 \
# -p 55

# plotHeatmap -m aged_AB.mat.gz \
# --outFileName aged_AB_kmeans4.png \
# --kmeans 4 \
# --outFileSortedRegions aged_AB_kmeans4_sortedRegions.txt \
# --outFileNameMatrix aged_AB_kmeans4_matrix.txt \
# --xAxisLabel "distance (bp)" \
# --samplesLabel Aged \
# --dpi 300