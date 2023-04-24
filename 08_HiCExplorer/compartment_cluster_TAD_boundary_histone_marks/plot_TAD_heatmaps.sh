prefix=/nas/homes/benyang/HiC

# computeMatrix reference-point \
# -R $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
# $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
# -S $prefix/04_FANC/without_KR_normalization/aged.merged_100kb_noKR.bw $prefix/04_FANC/without_KR_normalization/young.merged_100kb_noKR.bw \
# --outFileName aged_young_AB.mat.gz \
# --referencePoint center -b 500000 -a 500000 \
# --binSize 250 \
# -p 50

# plotHeatmap -m aged_young_AB.mat.gz \
# --outFileName aged_young_AB.png \
# --xAxisLabel "distance (bp)" \
# --samplesLabel Aged Young \
# --regionsLabel Aged Young \
# --dpi 300

# plotHeatmap -m aged_young_AB.mat.gz \
# --outFileName aged_young_AB_kmeans4.png \
# --kmeans 4 \
# --outFileSortedRegions aged_young_AB_kmeans4_sortedRegions.txt \
# --outFileNameMatrix aged_young_AB_kmeans4_matrix.txt \
# --xAxisLabel "distance (bp)" \
# --samplesLabel Aged Young \
# --dpi 300

############# Get A/B clusters
# Young
# computeMatrix reference-point \
# -R $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
# -S $prefix/04_FANC/without_KR_normalization/young.merged_100kb_noKR.bw \
# --outFileName young_AB.mat.gz \
# --referencePoint center -b 500000 -a 500000 \
# --binSize 10000 \
# -p 55

# plotHeatmap -m young_AB.mat.gz \
# --outFileName young_AB_kmeans4.png \
# --kmeans 4 \
# --outFileSortedRegions young_AB_kmeans4_sortedRegions.txt \
# --outFileNameMatrix young_AB_kmeans4_matrix.txt \
# --xAxisLabel "distance (bp)" \
# --samplesLabel Young \
# --dpi 300
# # Aged
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

############# Plot histone markers for each cluster
for feature in atac H3K4me3 H4K20me1 H3K27me3
do

    computeMatrix reference-point \
    -R $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/young_AB_kmeans4_cluster1.txt \
    $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/young_AB_kmeans4_cluster2.txt \
    $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/young_AB_kmeans4_cluster3.txt \
    $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/young_AB_kmeans4_cluster4.txt \
    -S $prefix/HistoneTracks/get_signal_bigwigs/$feature.count.rpkm.track.nodup.Young.bw \
    --outFileName young.$feature.bin10kb.mat.gz \
    --referencePoint center -b 500000 -a 500000 \
    --binSize 10000 \
    -p 55

    plotHeatmap -m young.$feature.bin10kb.mat.gz \
    --outFileName young.$feature.bin10kb.kmeans4.png \
    --xAxisLabel "distance (bp)" \
    --samplesLabel $feature \
    --regionsLabel "Cluster 1" "Cluster 2" "Cluster 3" "Cluster 4" \
    --dpi 300

    computeMatrix reference-point \
    -R $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/aged_AB_kmeans4_cluster1.txt \
    $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/aged_AB_kmeans4_cluster3.txt \
    $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/aged_AB_kmeans4_cluster2.txt \
    $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/aged_AB_kmeans4_cluster4.txt \
    -S $prefix/HistoneTracks/get_signal_bigwigs/$feature.count.rpkm.track.nodup.Aged.bw \
    --outFileName aged.$feature.bin10kb.mat.gz \
    --referencePoint center -b 500000 -a 500000 \
    --binSize 10000 \
    -p 55

    plotHeatmap -m aged.$feature.bin10kb.mat.gz \
    --outFileName aged.$feature.bin10kb.kmeans4.png \
    --xAxisLabel "distance (bp)" \
    --samplesLabel $feature \
    --regionsLabel "Cluster 1" "Cluster 2" "Cluster 3" "Cluster 4" \
    --dpi 300

done

######### Get matrices for A/B compartments separately

# computeMatrix reference-point \
# -R $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/young_AB_kmeans4_cluster1.txt \
# $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/young_AB_kmeans4_cluster2.txt \
# $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/young_AB_kmeans4_cluster3.txt \
# $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/young_AB_kmeans4_cluster4.txt \
# -S $prefix/04_FANC/without_KR_normalization/young.merged_100kb_noKR.bw \
# --outFileName young_AB.kmeans4.mat.gz \
# --referencePoint center -b 500000 -a 500000 \
# --binSize 10000 \
# -p 55

# computeMatrix reference-point \
# -R $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/aged_AB_kmeans4_cluster1.txt \
# $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/aged_AB_kmeans4_cluster3.txt \
# $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/aged_AB_kmeans4_cluster2.txt \
# $prefix/08_HiCExplorer/compartment_cluster_TAD_boundary_histone_marks/kmeans_clusters/aged_AB_kmeans4_cluster4.txt \
# -S $prefix/04_FANC/without_KR_normalization/aged.merged_100kb_noKR.bw \
# --outFileName aged_AB.kmeans4.mat.gz \
# --referencePoint center -b 500000 -a 500000 \
# --binSize 10000 \
# -p 55


# plotHeatmap -m ${age}_AB.mat.gz \
# --outFileName ${age}_AB_kmeans4.png \
# --xAxisLabel "distance (bp)" \
# --samplesLabel Aged \
# --regionsLabel "Cluster 1" "Cluster 2" "Cluster 3" "Cluster 4" \
# --dpi 300